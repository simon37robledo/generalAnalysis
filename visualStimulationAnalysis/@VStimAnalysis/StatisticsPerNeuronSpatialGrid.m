function results = StatisticsPerNeuronSpatialGrid(obj, params)
% StatisticsPerNeuronSpatialGrid - Spatial grid analysis of moving ball responses.
%
% Divides the screen into a GridSize × GridSize grid and analyses the response
% at each grid cell as the ball centre crosses it. Responses are compared across
% directions (main factor) and further split by other stimulus factors (offset,
% size, speed, luminosity) for downstream analysis.
%
% Pipeline:
%   1. Detect ball crossings per trial per grid cell (computeBallGridCrossings)
%   2. Extract GridAnalysisWindow ms response at each crossing
%   3. Per-direction permutation test with pooled null across directions
%   4. Best-direction observed stat, bias-corrected z-score
%   5. Per-cell z-scores split by direction × each non-direction factor
%
% Only applies to linearlyMovingBall. Other stimuli use the standard function.
%
% Reference:
%   Nichols & Holmes (2002) Human Brain Mapping 15:1-25

arguments (Input)
    obj
    params.nBoot              = 10000   % number of permutation iterations
    params.randomSeed         = 42      % fixed seed for reproducibility
    params.GridSize           = 9       % 9×9 = 81 grid cells
    params.GridAnalysisWindow = 200     % ms analysis window starting at ball crossing
    params.MinTrialsPerCell   = 3       % minimum trials per cell×direction×factor level
    params.ApplyFDR           = false   % Benjamini-Hochberg FDR correction reduces the number of false positives.
    params.Lock               = {}      % {name,value,...} keep only trials with these factor values; {} = all trials
end

% No caching here: the only caller, StatisticsPerNeuron, checks its params
% against the cache and saves this result, so a stale file can never be
% returned without that check.

% -------------------------------------------------------------------------
% Fix random seed for reproducible permutation results
% -------------------------------------------------------------------------
rng(params.randomSeed);

% -------------------------------------------------------------------------
% Load spike-sorted somatic units
% -------------------------------------------------------------------------
p     = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder,0,1);  % kilosort/phy output
label = string(p.label');                                             % unit quality labels
goodU = p.ic(:, label == 'good');                                     % somatic units only

if isempty(goodU)
    warning('%s has no somatic neurons.', obj.dataObj.recordingName);
    results = [];
    return
end

responseParams = obj.ResponseWindow;
nSpeeds        = numel(unique(obj.VST.speed));   % number of distinct speed conditions
winSize        = params.GridAnalysisWindow;      % analysis window in ms

% =========================================================================
% Main loop over speed conditions (Speed1, Speed2)
% =========================================================================
for s = 1:nSpeeds

    fieldName  = sprintf('Speed%d', s);
    C          = responseParams.(fieldName).C;           % stimulus category matrix
    trialTimes = C(:,1)';                                % trial onset times in ms
    stimDur    = responseParams.(fieldName).stimDur;     % full stimulus duration in ms

    % -------------------------------------------------------------------------
    % Frame-to-time conversion per trial
    % obj.VST.nFrames is [nSpeeds × nOffsets × nDirections] — frame count differs
    % across trials because speed changes trajectory duration.
    % For this speed condition (indexed by s), look up per-(offset, direction) frame
    % count, then map each trial's (offset, direction) to its correct frame count.
    % -------------------------------------------------------------------------
    nFramesFull = obj.VST.nFrames;


    % Build per-trial frame count using C(:,2)=direction and C(:,3)=offset
    uDirsAll        = unique(C(:,2));
    uOffsetsAll     = unique(C(:,3));
    nTrials         = size(C,1);
    nFramesPerTrial = zeros(nTrials, 1);

    if ndims(nFramesFull) == 3
        nFramesThisSpeed = reshape(nFramesFull(s, :, :),size(nFramesFull,2),size(nFramesFull,3));   % [nOffsets × nDirections]
    else
        nFramesThisSpeed = nFramesFull;                      % single-speed fallback
    end


    for t = 1:nTrials
        dIdx = find(uDirsAll    == C(t, 2));   % direction index
        oIdx = find(uOffsetsAll == C(t, 3));   % offset index
        nFramesPerTrial(t) = nFramesThisSpeed(oIdx, dIdx);
    end

    % Per-trial ms-per-frame conversion: [nTrials × 1]
    msPerFramePerTrial = stimDur ./ nFramesPerTrial;

    % -------------------------------------------------------------------------
    % Detect ball crossings per trial per grid cell
    % Returns: crossingFrame [nTrials × nCells], dwellFrames [nTrials × nCells],
    %          validGridPerDir [nCells × nDirs], nTrialsPerCellDir [nCells × nDirs]
    % -------------------------------------------------------------------------
    [crossingFrame, dwellFrames, validGridPerDir, nTrialsPerCellDir] = ...
        computeBallGridCrossings(obj, s, params);

    % -------------------------------------------------------------------------
    % Lock filter: keep only trials matching the locked factor values, using the
    % same helper as StatisticsPerNeuron so both tests select identical trials
    % -------------------------------------------------------------------------
    if ~isempty(params.Lock)
        assert(size(crossingFrame, 1) == size(C, 1), ...
            'Grid crossings (%d rows) and C (%d rows) are not row-aligned.', ...
            size(crossingFrame, 1), size(C, 1));                    % alignment invariant before subsetting
        keepMask = applyCategoryFocus(C, responseParams.colNames{1}(5:end), "all", params.Lock);  % locked trials
        C                  = C(keepMask, :);                         % conditions of kept trials
        trialTimes         = trialTimes(keepMask);                   % onsets of kept trials
        crossingFrame      = crossingFrame(keepMask, :);             % crossings of kept trials
        dwellFrames        = dwellFrames(keepMask, :);               % dwell frames of kept trials
        msPerFramePerTrial = msPerFramePerTrial(keepMask);           % frame duration of kept trials
        nTrials            = size(C, 1);                             % kept trial count

        % Validity and counts recomputed from kept trials only, as in
        % computeBallGridCrossings (directions ordered by unique(C(:,2)))
        uDirsKept         = unique(C(:, 2));                         % directions present after the lock
        validGridPerDir   = false(size(crossingFrame, 2), numel(uDirsKept));  % [nCells x nDirs]
        nTrialsPerCellDir = zeros(size(crossingFrame, 2), numel(uDirsKept));  % [nCells x nDirs]
        for d = 1:numel(uDirsKept)
            inDir = C(:, 2) == uDirsKept(d);                         % kept trials of this direction
            validGridPerDir(:, d)   = any(~isnan(crossingFrame(inDir, :)), 1)';  % cell crossed at least once
            nTrialsPerCellDir(:, d) = sum(~isnan(crossingFrame(inDir, :)), 1)';  % crossings per cell
        end
        fprintf('Lock kept %d trials for %s.\n', nTrials, fieldName);  % report the lock effect
    end

    % Dwell time in ms per trial per cell (per-trial frame rate)
    dwellTimeMs = dwellFrames .* msPerFramePerTrial;   % broadcast [nTrials × 1] across cells

    % Warn for cells with low mean dwell time
    meanDwellPerCell = mean(dwellTimeMs, 1, 'omitnan');     % [1 × nCells]
    lowDwellCells    = find(meanDwellPerCell < winSize/2 & meanDwellPerCell > 0);
    if ~isempty(lowDwellCells)
        fprintf(['Warning: %d grid cells have mean dwell time < %.0fms ' ...
                 '(half of analysis window). Interpret results at these ' ...
                 'cells cautiously.\n'], numel(lowDwellCells), winSize/2);
    end

    % -------------------------------------------------------------------------
    % Set up dimensions
    % -------------------------------------------------------------------------
    directions = C(:,2);                    % direction label per trial
    uDirs      = unique(directions);        % unique direction values
    nDirs      = numel(uDirs);              % number of directions
    nNeurons   = size(goodU, 2);            % number of somatic units
    nCells     = params.GridSize^2;         % total grid cells (e.g. 81)

    % -------------------------------------------------------------------------
    % Build full-duration response burst matrix ONCE for all trials
    % MrFull: [nTrials × nNeurons × stimDur] spike counts per ms bin
    % Then index into this matrix per grid cell using crossing times
    % -------------------------------------------------------------------------
    MrFull = BuildBurstMatrix(goodU, round(p.t), round(trialTimes), round(stimDur));

    % -------------------------------------------------------------------------
    % Baseline burst matrix: winSize ms from start of ITI preceding each trial
    % Mb: [nTrials × nNeurons × winSize]
    % baselines: [nTrials × nNeurons] — mean spikes/ms per trial per neuron
    % One baseline per trial, shared across all grid cells crossed in that trial
    % -------------------------------------------------------------------------
    MbMat = BuildBurstMatrix(goodU, round(p.t), ...
        round(trialTimes - winSize), ...
        winSize);
    baselines = mean(MbMat, 3);   % [nTrials × nNeurons]

    % Pooled baseline SD across all trials — stable z-score normalisation
    sdBase = std(baselines, 0, 1);   % [1 × nNeurons]

    % -------------------------------------------------------------------------
    % Extract per-cell response by indexing into MrFull
    % gridResponse: [nTrials × nCells × nNeurons]
    %   spikes/ms averaged over winSize starting at crossing time
    %   NaN where ball did not cross that cell on that trial
    % -------------------------------------------------------------------------
    gridResponse      = nan(nTrials, nCells, nNeurons);
    crossingMsInTrial = (crossingFrame - 1) .* msPerFramePerTrial;   % [nTrials × nCells]
    nTruncated        = 0;

    for t = 1:nTrials
        for c = 1:nCells
            if isnan(crossingFrame(t,c))
                continue   % ball did not cross this cell on this trial
            end

            startBin = round(crossingMsInTrial(t,c)) + 1;   % 1-indexed start bin in MrFull
            endBin   = startBin + winSize - 1;               % inclusive end bin

            if endBin > size(MrFull, 3)
                endBin     = size(MrFull, 3);                % truncate at stimulus end
                nTruncated = nTruncated + 1;                  % count for diagnostic
            end

            % Mean spikes/ms over the window for all neurons in this trial
            gridResponse(t, c, :) = mean(MrFull(t, :, startBin:endBin), 3);
        end
    end

    % Report truncation diagnostic
    if nTruncated > 0
        totalCrossings = sum(~isnan(crossingFrame(:)));
        fprintf(['Info: %d of %d trial-cell crossings (%.1f%%) had truncated ' ...
                 'analysis windows due to reaching stimulus end.\n'], ...
                 nTruncated, totalCrossings, 100*nTruncated/totalCrossings);
    end

    % -------------------------------------------------------------------------
    % Per-trial per-cell Diff: response minus trial's baseline
    % Broadcasting: baselines [nTrials × nNeurons] → [nTrials × 1 × nNeurons]
    % gridDiff: [nTrials × nCells × nNeurons], NaN where ball did not cross
    % -------------------------------------------------------------------------
    gridDiff = gridResponse - reshape(baselines, nTrials, 1, nNeurons);

    % =========================================================================
    % Per-direction observed stat and null distribution
    %
    % For each direction:
    %   Per-cell mean Diff across trials of that direction (ignoring NaN)
    %   Mask invalid cells (validGridPerDir)
    %   Max across valid cells = observed stat for this direction
    %   Sign-flip trials within direction → null distribution for this direction
    %
    % Pooled null = concatenation of per-direction null distributions
    % Observed overall stat = max across directions per neuron
    % =========================================================================
    obsStatPerDir = zeros(nDirs, nNeurons);                       % [nDirs × nNeurons]
    % Joint null over directions: ONE draw per bootstrap = max over cells AND
    % directions. Initialised at -Inf and updated by a running max across the
    % direction loop, so the null search space matches the observed statistic
    % (also max over cells AND directions) and FWER is controlled over both.
    nullMax = -Inf(params.nBoot, nNeurons);                       % [nBoot × nNeurons]

    for d = 1:nDirs
        trialsD = find(directions == uDirs(d));   % trial indices for this direction
        nTd     = numel(trialsD);

        % Extract Diff for this direction: [nTd × nCells × nNeurons]
        DiffD = gridDiff(trialsD, :, :);

        % Per-cell mean across trials for this direction (omits NaN from uncrossed cells)
        meanDiffD = squeeze(mean(DiffD, 1, 'omitnan'));   % [nCells × nNeurons]
        if nNeurons == 1
            meanDiffD = reshape(meanDiffD, nCells, 1);     % guard singleton collapse
        end

        % Mask cells invalid for this direction
        meanDiffDMasked                              = meanDiffD;
        meanDiffDMasked(~validGridPerDir(:,d), :)    = -Inf;

        % Observed max per neuron for this direction
        obsStatPerDir(d, :) = max(meanDiffDMasked, [], 1);

        % Sign-flip permutations within this direction
        % signs: [nTd × nBoot] — each column is one permutation
        signs = 2 * randi(2, nTd, params.nBoot) - 3;

        for b = 1:params.nBoot
            % Apply signs: broadcast [nTd × 1 × 1] across [nTd × nCells × nNeurons]
            DiffPerm = DiffD .* reshape(signs(:,b), nTd, 1, 1);

            % Per-cell mean under H0
            meanPerm = squeeze(mean(DiffPerm, 1, 'omitnan'));   % [nCells × nNeurons]
            if nNeurons == 1
                meanPerm = reshape(meanPerm, nCells, 1);
            end

            meanPerm(~validGridPerDir(:,d), :) = -Inf;

            % Join this direction's permuted cell-max into the running max across
            % directions for bootstrap b. After the direction loop completes,
            % nullMax(b,:) is the max over cells AND directions — the correct
            % max-statistic null draw matching the observed obsStat.
            nullMax(b, :) = max(nullMax(b, :), max(meanPerm, [], 1));
        end
    end

    % -------------------------------------------------------------------------
    % Overall observed stat and p-value
    % obsStat: max across cells AND directions per neuron
    % pVal: proportion of the joint (cells × directions) max null >= observed
    [obsStat, prefDirection] = max(obsStatPerDir, [], 1);   % [1 × nNeurons]
    pVal                     = mean(nullMax >= obsStat, 1);  % [1 × nNeurons]


    if  params.ApplyFDR 
        [pVal, ~, ~,~] = fdr_BH(pVal, 0.05);

    end

    % -------------------------------------------------------------------------
    % Bias-corrected z-score
    % z = (observed - expected null max) / pooled baseline SD
    % Subtracting mean(nullStatsAll) removes winner's curse inflation
    % -------------------------------------------------------------------------
    nullMean        = mean(nullMax, 1);               % [1 × nNeurons]
    z               = (obsStat - nullMean) ./ sdBase;  % [1 × nNeurons]
    z_mean          = obsStat;    
    z(sdBase == 0)  = 0;   % silent baseline — set to 0

    % -------------------------------------------------------------------------
    % Preferred grid cell per neuron at preferred direction
    % Identified from observed mean Diff across trials of preferred direction
    % -------------------------------------------------------------------------
    prefGridCell = zeros(1, nNeurons);
    for u = 1:nNeurons
        d         = prefDirection(u);
        trialsD   = find(directions == uDirs(d));
        meanDiffD = squeeze(mean(gridDiff(trialsD, :, u), 1, 'omitnan'));  % [nCells × 1]
        meanDiffD(~validGridPerDir(:,d)) = -Inf;
        [~, prefGridCell(u)] = max(meanDiffD);
    end

    % =========================================================================
    % Per-cell z-scores split by direction × each non-direction factor
    % C columns: 1=stimOn, 2=direction, 3=offset, 4=size, 5=speed, 6=luminosity
    % Output struct: one field per factor, each [nCells × nDirs × nLevels × nNeurons]
    % Direction is already a dimension of each array — not a separate field
    % =========================================================================
    factorCols  = 3:size(C,2);                        % non-direction factor columns
    factorNames = {'offset', 'size', 'speed', 'luminosity'};
    factorNames = factorNames(1:numel(factorCols));   % trim to available columns

    ZScorePerGrid = struct();

    for fIdx = 1:numel(factorCols)
        col     = factorCols(fIdx);
        uLevels = unique(C(:, col));                  % unique values for this factor
        nLevels = numel(uLevels);
        fName   = factorNames{fIdx};

        % Pre-allocate [nCells × nDirs × nLevels × nNeurons], NaN default
        zArr = nan(nCells, nDirs, nLevels, nNeurons);

        for d = 1:nDirs
            for lev = 1:nLevels
                % Trials matching this direction AND this factor level
                mask     = (directions == uDirs(d)) & (C(:,col) == uLevels(lev));
                trialsDL = find(mask);

                if numel(trialsDL) < params.MinTrialsPerCell
                    continue   % too few trials — leave NaN
                end

                DiffDL           = gridDiff(trialsDL, :, :);            % [nTdl × nCells × nNeurons]
                nTrialsPerCellDL = sum(~isnan(DiffDL(:,:,1)), 1);       % [1 × nCells]

                meanDL = squeeze(mean(DiffDL, 1, 'omitnan'));            % [nCells × nNeurons]
                if nNeurons == 1
                    meanDL = reshape(meanDL, nCells, 1);
                end

                % Normalise by pooled baseline SD — same for all cells
                zDL = meanDL ./ reshape(sdBase, 1, nNeurons);            % [nCells × nNeurons]

                % Mask cells with too few crossings for this (direction, factor level)
                zDL(nTrialsPerCellDL < params.MinTrialsPerCell, :) = NaN;

                % Store in 4D array
                zArr(:, d, lev, :) = reshape(zDL, nCells, 1, 1, nNeurons);
            end
        end

        ZScorePerGrid.(fName) = zArr;   % [nCells × nDirs × nLevels × nNeurons]
    end

    % =========================================================================
    % Split-half cross-validated peak delta rate (spikes/s)
    % Cell x direction is chosen on the odd repetitions of each stimulus
    % condition and measured on the even ones (and vice versa), so the value is
    % free of the winner's-curse inflation of z_mean (a plain max over cells x
    % directions). Splitting within each condition (direction x offset x size
    % ...) puts every crossed cell's trials into both halves.
    % =========================================================================
    [~, ~, condID] = unique(C(:, 2:end), 'rows');      % condition index per trial (all stimulus factors)
    repRank = zeros(nTrials, 1);                       % repetition rank of each trial within its condition
    for k = 1:max(condID)
        rows          = find(condID == k);             % trials of this condition, in trial order
        repRank(rows) = (1:numel(rows))';              % 1, 2, 3, ... within the condition
    end
    nRepMin = min(accumarray(condID, 1));              % fewest repetitions of any condition
    if nRepMin < 2                                     % a single repetition cannot be split
        warning('%s %s: a condition has %d repetition(s); split-half delta rate is NaN for its cells.', ...
            obj.dataObj.recordingName, fieldName, nRepMin);
    end
    halfID      = 2 - mod(repRank, 2);                 % 1 = odd repetition, 2 = even repetition
    [~, ~, dirID] = unique(directions);                % direction index per trial (selection group)
    DeltaRateCV = splitHalfDelta(permute(gridDiff, [1 3 2]), dirID, halfID);  % [1 × nNeurons] spikes/s

    % =========================================================================
    % Store results for this speed condition
    % =========================================================================
    S.(fieldName).pvalsResponse         = pVal;                % [1 × nNeurons]
    S.(fieldName).ZScoreU               = z;                    % [1 × nNeurons] bias-corrected z
    S.(fieldName).prefDirection         = prefDirection;        % [1 × nNeurons]
    S.(fieldName).prefGridCell          = prefGridCell;         % [1 × nNeurons]
    S.(fieldName).ObsResponse           = gridResponse;         % [nTrials × nCells × nNeurons]
    S.(fieldName).ObsBaseline           = baselines;            % [nTrials × nNeurons]
    S.(fieldName).ObsDiff               = gridDiff;             % [nTrials × nCells × nNeurons]
    S.(fieldName).validGrid             = validGridPerDir;      % [nCells × nDirs]
    S.(fieldName).dwellTimeMs           = dwellTimeMs;           % [nTrials × nCells]
    S.(fieldName).meanDwellPerCell      = meanDwellPerCell;      % [1 × nCells]
    S.(fieldName).nTrialsPerCellDir     = nTrialsPerCellDir;     % [nCells × nDirs]
    S.(fieldName).ZScorePerGrid         = ZScorePerGrid;         % struct of 4D arrays
    S.(fieldName).gridSize              = params.GridSize;       % scalar, grid dimensions
    S.(fieldName).z_mean               = z_mean*1000;
    S.(fieldName).DeltaRateCV           = DeltaRateCV;           % [1 × nNeurons] split-half peak delta rate, spikes/s

    S.params = params;   % store parameters for reproducibility

end % end speed loop

% --- Return (StatisticsPerNeuron tags and saves the result) ---
results = S;

end % end main function