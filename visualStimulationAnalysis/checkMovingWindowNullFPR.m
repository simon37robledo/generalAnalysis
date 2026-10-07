function [fprTable, deltaTable, info] = checkMovingWindowNullFPR(expNum, params)
% checkMovingWindowNullFPR  Null calibration of the responsiveness tests.
%
%   [fprTable, deltaTable, info] = checkMovingWindowNullFPR(97)
%
%   Simulates neurons with NO stimulus effect (homogeneous Bernoulli spiking at
%   the same rate before and during the stimulus) on the real trial layout of one
%   recording, runs each responsiveness test exactly as the production code does,
%   and reports the false-positive rate at alpha. A valid test gives ~alpha.
%
%   Tests (production lines mirrored are cited next to each step):
%     MB_slide     sliding-window MB test, per-neuron null      (StatisticsPerNeuron v2)
%     MB_slide_v1  same, null pooled across neurons            (pre-fix v1 behaviour)
%     MB_grid      spatial-grid MB test                         (StatisticsPerNeuronSpatialGrid)
%     RG           fixed-window RG test, per-neuron null        (StatisticsPerNeuron v3)
%     RG_v1        same, null pooled across neurons            (pre-fix v1 behaviour)
%
%   Each simulated recording ("chunk") mixes neurons at every firing rate, so a
%   null pooled across neurons is exposed rather than hidden.
%
%   deltaTable reports the split-half cross-validated peak delta rate (spikes/s)
%   under the null; an unbiased estimator averages ~0.
%
%   Reads the recording's cached ResponseWindow only; writes nothing to disk.

arguments
    expNum (1,1) double = 97                    % recording supplying the trial layout
    params.ratesHz double = [0.5 1 2 5 10 20]   % simulated firing rates (Hz)
    params.nPerRate double = 240                % simulated neurons per rate
    params.perRateInChunk double = 5            % neurons per rate inside one simulated recording
    params.nBoot double = 2000                  % sign-flip permutations per test
    params.alpha double = 0.05                  % responsiveness threshold (p < alpha)
    params.lock cell = {"luminosities", 255}    % Figure 1 Lock for both stimuli
    params.slideWin double = 250                % Figure 1 durationWindow = BaseRespWindow (ms)
    params.gridWin double = 300                 % GridAnalysisWindow used by StatisticsPerNeuron:154
    params.gridSize double = 9                  % 9 x 9 grid, as in StatisticsPerNeuron:154
    params.rgWin double = 300                   % RG response (from onset) and baseline window: BaseRespWindow forced to 300 in grid mode
    params.seed double = 1                      % RNG seed for reproducibility
    params.selfTest logical = true              % check vectorised grid null against production loop
end

rng(params.seed);                                                  % reproducible simulation

% =========================================================================
% 1. Real trial layouts (cached ResponseWindow, bare call = load only)
% =========================================================================
NP   = loadNPclassFromTable(expNum);                               % recording object
vsMB = linearlyMovingBallAnalysis(NP);                             % moving-ball analysis object
rwMB = vsMB.ResponseWindow;                                        % cached MB response-window struct
assert(isfield(rwMB, 'Speed1'), 'Recording %d has no MB Speed1.', expNum);   % MB layout required

CMB       = rwMB.Speed1.C;                                         % [nTrials x (1+nParams)] MB conditions
stimDurMB = round(rwMB.Speed1.stimDur);                            % MB response length, round(effectiveStimDur) as in BuildBurstMatrix
[keepMB, catMB] = applyCategoryFocus(CMB, rwMB.colNames{1}(5:end), "all", params.lock);  % same focus call as StatisticsPerNeuron
keptMB          = find(keepMB);                                    % trials surviving the lock
[catMBs, ordMB] = sort(catMB);                                     % contiguous category blocks (StatisticsPerNeuron focus block)
rowMB           = keptMB(ordMB);                                   % real trial index of each simulated MB trial
cntMB           = accumarray(catMBs, 1);                           % trials per MB category
assert(all(cntMB == cntMB(1)), 'MB focus is unbalanced: %s', mat2str(cntMB'));  % balance required by the vectorised null
tcMB   = cntMB(1);                                                 % trials per MB category
nCatMB = numel(cntMB);                                             % number of MB categories
nTrMB  = numel(rowMB);                                             % number of simulated MB trials
assert(tcMB >= 2, 'Split-half needs >= 2 trials per MB category.');  % odd/even halves must both exist

vsRG = rectGridAnalysis(NP);                                       % rect-grid analysis object
rwRG = vsRG.ResponseWindow;                                        % cached RG response-window struct
CRG       = rwRG.C;                                                % RG conditions
stimDurRG = round(rwRG.stimDur);                                   % full RG stimulus length (reported only)
assert(params.rgWin <= stimDurRG, 'rgWin exceeds the RG stimulus duration.');  % fixed window must fit in the stimulus
[keepRG, catRG] = applyCategoryFocus(CRG, rwRG.colNames{1}(5:end), "all", params.lock);  % same focus as production
assert(numel(catRG) == sum(keepRG), 'RG catID is not one entry per kept trial.');  % alignment invariant
catRGs = sort(catRG);                                              % category of each kept trial, block-sorted
cntRG  = accumarray(catRGs, 1);                                    % trials per RG category
assert(all(cntRG == cntRG(1)), 'RG focus is unbalanced: %s', mat2str(cntRG'));  % balance required
tcRG   = cntRG(1);                                                 % trials per RG category
nCatRG = numel(cntRG);                                             % number of RG categories
nTrRG  = numel(catRGs);                                            % number of simulated RG trials

% --- Grid crossing times on the MB trials (StatisticsPerNeuronSpatialGrid:84-118, 166) ---
gp.GridSize = params.gridSize;                                     % only field computeBallGridCrossings reads
[crossF, ~, validG] = computeBallGridCrossings(vsMB, 1, gp);      % [nTrials x nCells] crossing frames, [nCells x nDirs] validity
assert(size(crossF, 1) == size(CMB, 1), 'Grid crossings and C are not row-aligned.');  % alignment invariant
nFr = vsMB.VST.nFrames;                                            % frames per (speed, offset, direction)
if ndims(nFr) == 3                                                 % multi-speed layout
    nFr = reshape(nFr(1, :, :), size(nFr, 2), size(nFr, 3));       % Speed1 slice [nOffsets x nDirs]
end
uDirAll = unique(CMB(:, 2));                                       % direction levels (C column 2)
uOffAll = unique(CMB(:, 3));                                       % offset levels (C column 3)
nFrTrial = zeros(size(CMB, 1), 1);                                 % frames of each trial's trajectory
for t = 1:size(CMB, 1)                                             % per-trial frame-count lookup (SpatialGrid:100-104)
    nFrTrial(t) = nFr(uOffAll == CMB(t, 3), uDirAll == CMB(t, 2)); % frames for this trial's offset x direction
end
msPerFrame = rwMB.Speed1.stimDur ./ nFrTrial;                      % ms per frame per trial (SpatialGrid:107)
crossMs    = (crossF(rowMB, :) - 1) .* msPerFrame(rowMB);          % crossing time in ms, simulated-trial order (SpatialGrid:166)
[~, ~, dirIdx] = unique(CMB(rowMB, 2));                            % direction index of each simulated MB trial
nCells = size(crossF, 2);                                          % number of grid cells

% --- Odd/even halves by repetition rank within each category (split-half) ---
halfMB = 2 - mod(repmat((1:tcMB)', nCatMB, 1), 2);                 % 1 = odd repetition, 2 = even repetition
halfRG = 2 - mod(repmat((1:tcRG)', nCatRG, 1), 2);                 % same for RG

info = struct('expNum', expNum, 'nCatMB', nCatMB, 'tcMB', tcMB, 'stimDurMB', stimDurMB, ...
    'nCatRG', nCatRG, 'tcRG', tcRG, 'stimDurRG', stimDurRG, 'nCells', nCells);  % layout summary
fprintf('Layout exp %d: MB %d cats x %d trials, %d ms | RG %d cats x %d trials, %d ms | grid %d cells\n', ...
    expNum, nCatMB, tcMB, stimDurMB, nCatRG, tcRG, stimDurRG, nCells);  % echo the layout used

% =========================================================================
% 2. Simulate recordings with mixed firing rates and run every test
% =========================================================================
nRates    = numel(params.ratesHz);                                 % number of simulated rates
nU        = nRates * params.perRateInChunk;                        % neurons per simulated recording
nChunks   = ceil(params.nPerRate / params.perRateInChunk);         % recordings needed for nPerRate per rate
baseLen   = max(params.slideWin, params.gridWin);                  % simulated pre-stimulus window (ms)
rateOfU   = repelem(params.ratesHz, params.perRateInChunk);        % firing rate of each neuron in a chunk (Hz)
pLam      = reshape(rateOfU / 1000, 1, nU);                        % spike probability per 1 ms bin

tests     = ["MB_slide", "MB_slide_v1", "MB_grid", "RG", "RG_v1"]; % test labels
P         = nan(numel(tests), nChunks * nU);                       % p-values [test x simulated neuron]
layouts   = ["MB_slide", "MB_grid", "RG"];                         % split-half layouts
Dlt       = nan(numel(layouts), nChunks * nU);                     % split-half delta (spikes/s)
rateAll   = repmat(rateOfU, 1, nChunks);                           % rate label of every simulated neuron

for k = 1:nChunks                                                  % one simulated recording per iteration
    cols = (k-1)*nU + (1:nU);                                      % output columns of this chunk

    % ---- MB spikes: [trials x neurons x ms], Bernoulli per 1 ms bin ----
    Mr = single(rand(nTrMB, nU, stimDurMB, 'single') < pLam);      % response spikes over the full stimulus
    Mb = single(rand(nTrMB, nU, baseLen,  'single') < pLam);       % pre-stimulus spikes (last bin = just before onset)
    MbS = Mb(:, :, end-params.slideWin+1:end);                     % sliding-test baseline: [onset-250, onset) (StatisticsPerNeuron:336)

    % ---- (a) MB sliding window (StatisticsPerNeuron:377-381) ----
    mrMov  = movmean(Mr,  params.slideWin, 3, 'Endpoints', 'discard');  % response window means at every position
    respMW = max(mrMov, [], 3);                                    % per-trial best response window (line 379)
    baseMW = max(movmean(MbS, params.slideWin, 3, 'Endpoints', 'discard'), [], 3);  % per-trial best baseline window (line 380)
    [P(1, cols), P(2, cols)] = maxCatSignFlip(double(respMW - baseMW), tcMB, nCatMB, params.nBoot);  % v2 and v1 p-values
    Dlt(1, cols) = splitHalfDelta(mrMov - mean(MbS, 3), catMBs, halfMB);  % split-half delta, 250 ms windows
    clear mrMov                                                    % free the largest array before the grid step

    % ---- (b) MB spatial grid (StatisticsPerNeuronSpatialGrid:151-201) ----
    gridDiff = buildGridDiff(Mr, mean(Mb(:, :, end-params.gridWin+1:end), 3), crossMs, params.gridWin);  % [trials x cells x neurons]
    signsGrid = cell(1, max(dirIdx));                              % one sign matrix per direction, drawn as in SpatialGrid:244
    for d = 1:max(dirIdx)                                          % draw in direction order
        signsGrid{d} = 2 * randi(2, sum(dirIdx == d), params.nBoot) - 3;  % +-1 per trial and permutation
    end
    P(3, cols) = gridSignFlip(gridDiff, dirIdx, validG, signsGrid);% grid p-values
    Dlt(2, cols) = splitHalfDelta(permute(gridDiff, [1 3 2]), dirIdx, halfMB);  % split-half over cell x direction

    if params.selfTest && k == 1                                   % once: vectorised grid null == production loop
        nTest = min(50, params.nBoot);                             % short test run
        sub   = cellfun(@(s) s(:, 1:nTest), signsGrid, 'UniformOutput', false);  % same signs for both versions
        tic; [~, nLoop, oLoop] = gridSignFlipLoop(gridDiff, dirIdx, validG, sub); tLoop = toc;  % production-style loop, timed
        [~, nVec, oVec] = gridSignFlip(gridDiff, dirIdx, validG, sub);  % vectorised version
        % Compare the null and observed values (not p): summation order differs,
        % so exact ties in nullMax >= obs could flip by floating-point rounding.
        assert(max(abs(nLoop(:) - nVec(:))) < 1e-9 && max(abs(oLoop - oVec)) < 1e-9, ...
            'Vectorised grid null differs from the production loop.');  % hard check
        info.gridLoopSecPer10kBootPerNeuron = tLoop / nTest * 1e4 / nU;  % runtime extrapolation for the real loop
        fprintf('Self-test passed: vectorised grid null matches the production loop.\n');  % report
    end
    clear Mr Mb MbS gridDiff                                       % free MB arrays

    % ---- (c) RG fixed-window mean (StatisticsPerNeuron, MovingWindowPVal=false) ----
    MrR = single(rand(nTrRG, nU, params.rgWin, 'single') < pLam);  % RG response spikes, fixed window from onset
    MbR = single(rand(nTrRG, nU, params.rgWin, 'single') < pLam);  % RG baseline spikes, same length
    [P(4, cols), P(5, cols)] = maxCatSignFlip(double(mean(MrR, 3) - mean(MbR, 3)), tcRG, nCatRG, params.nBoot);  % v2 and v1
    Dlt(3, cols) = splitHalfDelta(mean(MrR, 3) - mean(MbR, 3), catRGs, halfRG);  % split-half over positions (one fixed window)
    clear MrR MbR                                                  % free RG arrays

    fprintf('  simulated recording %d/%d done\n', k, nChunks);     % progress
end

% =========================================================================
% 3. Summaries
% =========================================================================
fprTable = table();                                                % false-positive rate per test x rate
for ti = 1:numel(tests)
    for r = params.ratesHz
        m = rateAll == r;                                          % neurons simulated at this rate
        nSig = sum(P(ti, m) < params.alpha);                       % false positives at alpha
        [phat, ci] = binofit(nSig, sum(m));                        % Clopper-Pearson 95% CI
        fprTable = [fprTable; table(tests(ti), r, sum(m), nSig, phat, ci(1), ci(2), ...
            'VariableNames', {'test','rateHz','n','nSig','FPR','ciLo','ciHi'})];  %#ok<AGROW>
    end
end

deltaTable = table();                                              % split-half delta under the null
for li = 1:numel(layouts)
    for r = params.ratesHz
        v = Dlt(li, rateAll == r);                                 % delta of neurons at this rate
        nBad = sum(isnan(v));                                      % neurons with an unmeasurable half
        if nBad > 0                                                % detectable error condition -> explicit warning
            warning('%s at %g Hz: %d neurons had NaN split-half delta.', layouts(li), r, nBad);
        end
        v = v(~isnan(v));                                          % keep measurable neurons
        deltaTable = [deltaTable; table(layouts(li), r, numel(v), mean(v), std(v)/sqrt(numel(v)), ...
            'VariableNames', {'layout','rateHz','n','meanDeltaHz','seDeltaHz'})];  %#ok<AGROW>
    end
end

disp(fprTable); disp(deltaTable);                                  % print both summaries
end

% =========================================================================
% Local functions
% =========================================================================
function [pV2, pV1] = maxCatSignFlip(D, tc, nCats, nBoot)
% Max-category sign-flip test of StatisticsPerNeuron:503-589 (maxCategory=true,
% FilterEmptyResponses=false so every category is valid). D: [nTrials x nNeurons],
% rows in contiguous category blocks. Returns the v2 (per-neuron) p-value and the
% v1 (pooled through the old reshape) p-value from the SAME permutations.
nU       = size(D, 2);                                             % neurons
signs    = 2 * randi(2, size(D, 1), nBoot) - 3;                    % line 503
signsR   = reshape(signs, tc, nCats, nBoot);                       % [trialsCat x nCats x nBoot]
DR       = reshape(D, tc, nCats, nU);                              % [trialsCat x nCats x nNeurons]
catMeans = reshape(mean(DR, 1), nCats, nU);                        % observed category means
obs      = max(catMeans, [], 1);                                   % ObsStat, line 562
cmAll    = pagemtimes(permute(DR, [3 1 2]), permute(signsR, [1 3 2])) / tc;  % [nNeurons x nBoot x nCats]
mx       = max(permute(cmAll, [3 1 2]), [], 1);                    % [1 x nNeurons x nBoot]
nullV2   = reshape(mx, nU, nBoot).';                               % fixed reshape, line 575
nullV1   = reshape(mx, nBoot, nU);                                 % pre-fix reshape (pooled across neurons)
pV2      = mean(nullV2 >= obs, 1);                                 % line 589
pV1      = mean(nullV1 >= obs, 1);                                 % same formula on the pooled null
end

function gridDiff = buildGridDiff(Mr, baseline, crossMs, win)
% Per-trial per-cell response minus baseline (StatisticsPerNeuronSpatialGrid:165-201).
% Mr: [nTrials x nNeurons x T]; baseline: [nTrials x nNeurons] mean rate.
[nTr, nU, T] = size(Mr);                                           % dimensions
nCells   = size(crossMs, 2);                                       % grid cells
gridDiff = nan(nTr, nCells, nU);                                   % NaN where the ball never crossed
for t = 1:nTr
    for c = 1:nCells
        if isnan(crossMs(t, c)), continue; end                     % no crossing on this trial (line 171)
        s = round(crossMs(t, c)) + 1;                              % start bin (line 175)
        e = min(s + win - 1, T);                                   % end bin, truncated at stimulus end (lines 176-181)
        gridDiff(t, c, :) = reshape(mean(Mr(t, :, s:e), 3), 1, 1, nU) ...
            - reshape(baseline(t, :), 1, 1, nU);                   % window mean minus trial baseline (lines 184, 201)
    end
end
end

function [p, nullMax, obs] = gridSignFlip(gridDiff, dirIdx, validG, signsByDir)
% Vectorised form of the joint cells x directions max-statistic null
% (StatisticsPerNeuronSpatialGrid:215-271). Sign-flipped omitnan means are
% computed as (signs' * zero-filled data) ./ crossing counts.
[~, nCells, nU] = size(gridDiff);                                  % dimensions
nBoot   = size(signsByDir{1}, 2);                                  % permutations
obsDir  = -Inf(max(dirIdx), nU);                                   % observed max per direction
nullMax = -Inf(nBoot, nU);                                         % joint running-max null
for d = 1:max(dirIdx)
    Dd    = gridDiff(dirIdx == d, :, :);                           % this direction's trials
    ok    = ~isnan(Dd(:, :, 1));                                   % crossing pattern (same for every neuron)
    nC    = sum(ok, 1);                                            % crossings per cell [1 x nCells]
    D0    = Dd; D0(isnan(D0)) = 0;                                 % zero-fill uncrossed entries
    m     = reshape(sum(D0, 1), nCells, nU) ./ nC';                % observed omitnan cell means (line 230)
    m(~validG(:, d), :) = -Inf;                                    % mask invalid cells (line 237)
    obsDir(d, :) = max(m, [], 1);                                  % observed max for this direction (line 240)
    S     = signsByDir{d}.' * reshape(D0, size(D0, 1), []);        % [nBoot x nCells*nNeurons] sign-flipped sums
    mp    = reshape(S, nBoot, nCells, nU) ./ nC;                   % sign-flipped omitnan means (line 251)
    mp(:, ~validG(:, d), :) = -Inf;                                % mask invalid cells (line 256)
    nullMax = max(nullMax, reshape(max(mp, [], 2), nBoot, nU));    % running max over cells and directions (line 262)
end
obs = max(obsDir, [], 1);                                          % observed max over cells and directions (line 270)
p   = mean(nullMax >= obs, 1);                                     % line 271
end

function [p, nullMax, obs] = gridSignFlipLoop(gridDiff, dirIdx, validG, signsByDir)
% Literal production loop (StatisticsPerNeuronSpatialGrid:222-271), used only to
% verify gridSignFlip and to time the real implementation.
[~, nCells, nU] = size(gridDiff);                                  % dimensions
nBoot   = size(signsByDir{1}, 2);                                  % permutations
obsPerDir = zeros(max(dirIdx), nU);                                % observed max per direction
nullMax = -Inf(nBoot, nU);                                         % joint null
for d = 1:max(dirIdx)
    DiffD = gridDiff(dirIdx == d, :, :);                           % line 227
    nTd   = size(DiffD, 1);                                        % trials in this direction
    mD    = reshape(mean(DiffD, 1, 'omitnan'), nCells, nU);        % line 230
    mD(~validG(:, d), :) = -Inf;                                   % line 237
    obsPerDir(d, :) = max(mD, [], 1);                              % line 240
    for b = 1:nBoot
        DiffPerm = DiffD .* reshape(signsByDir{d}(:, b), nTd, 1, 1);  % line 248
        meanPerm = reshape(mean(DiffPerm, 1, 'omitnan'), nCells, nU);  % line 251
        meanPerm(~validG(:, d), :) = -Inf;                         % line 256
        nullMax(b, :) = max(nullMax(b, :), max(meanPerm, [], 1));  % line 262
    end
end
obs = max(obsPerDir, [], 1);                                       % line 270
p   = mean(nullMax >= obs, 1);                                     % line 271
end

function dlt = splitHalfDelta(X, grp, half)
% Split-half cross-validated peak delta. X: [nTrials x nNeurons x nOptions]
% (options = window positions or grid cells; NaN = not measurable), grp: selection
% group per trial (category or direction), half: 1/2 per trial. Select the
% (group, option) maximising the mean on one half, measure it on the other half,
% swap, average; returns spikes/s [1 x nNeurons].
[nTr, nU, nOpt] = size(X);                                         % dimensions
nGrp = max(grp);                                                   % selection groups
ok   = ~isnan(X);                                                  % measurable entries
X(~ok) = 0;                                                        % zero-fill for summation
M = cell(1, 2);                                                    % group x option means per half
for h = 1:2
    G    = double(sparse(grp(half == h), find(half == h), 1, nGrp, nTr));  % [nGrp x nTrials] membership in this half
    sums = full(G * reshape(double(X),  nTr, []));                 % sums  [nGrp x nNeurons*nOpt]
    cnts = full(G * reshape(double(ok), nTr, []));                 % counts [nGrp x nNeurons*nOpt]
    M{h} = reshape(permute(reshape(sums ./ cnts, nGrp, nU, nOpt), [1 3 2]), nGrp*nOpt, nU);  % [nGrp*nOpt x nNeurons]
end
est = zeros(2, nU);                                                % held-out estimate per direction of the split
for h = 1:2
    sel = M{h}; sel(isnan(sel)) = -Inf;                            % unmeasurable options cannot be selected
    [~, idx] = max(sel, [], 1);                                    % selected (group, option) per neuron
    est(h, :) = M{3-h}(idx + (0:nU-1) * nGrp * nOpt);              % measured on the other half
end
dlt = mean(est, 1) * 1000;                                         % average of both splits, spikes/ms -> spikes/s
end
