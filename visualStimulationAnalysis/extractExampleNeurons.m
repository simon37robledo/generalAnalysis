function T = extractExampleNeurons(cacheFile, params)
% extractExampleNeurons  Pick MIDDLE-of-distribution example neurons per stimulus.
%
%   From the per-neuron table AllExpAnalysis caches (TableStimComp, which holds
%   only responsive neurons), select the few neurons whose response metric sits
%   closest to the MEDIAN of each stimulus's distribution — i.e. representative,
%   not extreme, examples for a figure. For each pick it resolves the Phy cluster
%   ID and the globally-unique experiment ID (by loading that neuron's NP), so you
%   can plot the rasters yourself and then hand-build a markUnits cell.
%
%   IMPORTANT — responsiveness gate: TableStimComp stores the OR-mask neuron set
%   (responsive to at least ONE compared stimulus) replicated across every stimulus
%   column. So a stimulus's block contains neurons responsive only to the OTHER
%   stimulus, whose z here is ~0. Before taking the median, each stimulus is gated to
%   neurons with Z-score >= params.minZ (default 1.64, ~one-tailed p=0.05 for the
%   right-tailed enhancement test) so the "middle" is a genuinely responsive example.
%
%   The returned table exposes, per neuron:
%     animal, realInsertion (original insertion, resets across animals),
%     Exp (experiment ID, globally unique — the loadNPclassFromTable key),
%     PhyID (Phy cluster ID), plus NeurID and the metric value for reference.
%
%   EXAMPLE
%     T = extractExampleNeurons( ...
%           "W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_1-97_Combined_MB-RG.mat");
%     % then plot T.PhyID rasters yourself and build markUnits from the winners.

arguments
    cacheFile                                       % path to AllExpAnalysis *.mat cache, OR a preloaded struct, OR a TableStimComp table
    params.metric      string  = "Z-score"          % ranking column: "Z-score" (default) or "SpkR"
    params.nPer        double  = 5                   % candidate neurons to return per stimulus (around the median)
    params.percentile  double  = 50                  % "middle" of the distribution = 50th percentile (exposed for flexibility)
    params.minZ        double  = 1.64                 % responsiveness gate on Z-score for THIS stimulus (~one-tailed p=0.05, right-tailed enhancement). See note below.
    params.stimuli     string  = ["MB","RG"]          % stimulus categories to sample
    params.expList     double  = []                  % experiment integers; only needed if a bare table/struct without expList is passed
    params.plot        logical = false               % also draw each candidate's raster via that experiment's plotRaster
    params.PaperFig    logical = false               % when true, save the rasters into a "neuron candidates_stimA_stimB" subfolder
    params.outRoot     string  = "W:\Large_scale_mapping_NP\Paper_figs" % root that holds the save subfolder
end

% =========================================================================
% 1. Load TableStimComp and the experiment list
% =========================================================================
if istable(cacheFile)                                % caller passed the table directly
    T0      = cacheFile;                             % use it as-is
    expList = params.expList;                        % expList must then come from the argument
elseif isstruct(cacheFile)                           % caller passed a preloaded cache struct
    assert(isfield(cacheFile,'TableStimComp'), 'Struct lacks TableStimComp.'); % must contain the per-neuron table
    T0      = cacheFile.TableStimComp;               % pull the per-neuron table
    expList = pickExpList(cacheFile, params.expList);% prefer struct.expList, else the argument
else                                                 % caller passed a path to the .mat cache
    S = load(cacheFile);                             % load the cached AllExpAnalysis output
    assert(isfield(S,'TableStimComp'), 'Cache lacks TableStimComp — not an AllExpAnalysis cache?'); % validate contents
    T0      = S.TableStimComp;                       % per-neuron table (responsive neurons only)
    expList = pickExpList(S, params.expList);        % experiment integers saved alongside the table
end

% =========================================================================
% 2. Validate columns and inputs (fail loud, per project standard)
% =========================================================================
metric   = params.metric;                            % ranking column name (may contain a hyphen)
needCols = ["animal","realInsertion","insertion","stimulus","NeurID", metric]; % columns this function relies on
missing  = needCols(~ismember(needCols, string(T0.Properties.VariableNames))); % any required column absent?
assert(isempty(missing), 'TableStimComp is missing required column(s): %s', strjoin(missing, ', ')); % clear error listing them
assert(~isempty(expList), ...                        % Phy/Exp resolution cannot proceed without the experiment list
    'expList is empty — pass a cache with S.expList or supply params.expList so Phy IDs can be resolved.');

% =========================================================================
% 3. Select the middle-of-distribution neurons for each stimulus
% =========================================================================
picks = table();                                     % accumulator for chosen rows across all stimuli
for stim = params.stimuli                            % one stimulus at a time (e.g. "MB", then "RG")
    rows = T0(T0.stimulus == stim, :);               % all responsive neurons for this stimulus
    if isempty(rows)                                 % stimulus not present (e.g. category/level-mode composite labels)
        warning('No rows for stimulus "%s" — skipping.', stim); %#ok<*WNTAG> % tell the user and move on
        continue                                     % nothing to sample for this stimulus
    end

    % Responsiveness gate: TableStimComp replicates the OR-mask neuron set across
    % EVERY stimulus column, so a stimulus's block includes neurons responsive only
    % to the OTHER stimulus (their z here can be ~0). Restrict to neurons genuinely
    % responsive to THIS stimulus before taking the median, else the "middle" lands
    % on a non-responsive unit (e.g. RG median z ~0.2 without this gate).
    zcol  = rows.('Z-score');                        % responsiveness measure for THIS stimulus (always Z-score, even if ranking by SpkR)
    m     = rows.(char(metric));                     % ranking metric per neuron (char() handles the hyphen in "Z-score")
    keep  = ~isnan(m) & ~isnan(zcol) & zcol >= params.minZ; % usable AND responsive to this stimulus
    nDrop = sum((zcol < params.minZ) & ~isnan(zcol));% how many were gated out as non-responsive to this stimulus
    rows  = rows(keep, :);                           % keep only responsive, usable rows
    m     = m(keep);                                 % keep matching ranking values
    if nDrop > 0                                     % report the gating so the pool size is transparent
        fprintf('  %s: gated out %d neuron(s) with Z-score < %.2f (not responsive to %s).\n', stim, nDrop, params.minZ, stim);
    end
    if isempty(m)                                    % nothing passed the gate
        warning('No neuron passes the Z-score >= %.2f responsiveness gate for "%s" — skipping.', params.minZ, stim); % report and skip
        continue                                     % cannot rank an empty set
    end

    target = prctile(m, params.percentile);          % the "middle" value of this stimulus's distribution
    [~, ord] = sort(abs(m - target), 'ascend');      % order neurons by closeness to that middle value
    n      = min(params.nPer, numel(ord));           % cannot return more than exist
    if n < params.nPer                               % fewer neurons than requested
        warning('Only %d neuron(s) available for "%s" (requested %d).', n, stim, params.nPer); % note the shortfall
    end
    sel    = ord(1:n);                               % indices of the n neurons nearest the median

    % --- build the per-pick block; str2double(string()) NEVER double() on categoricals ---
    blk = table();                                                    % this stimulus's chosen block
    blk.stimulus      = repmat(string(stim), n, 1);                   % raw stimulus label (e.g. "MB")
    blk.animal        = string(rows.animal(sel));                     % animal ID (categorical -> string)
    blk.realInsertion = str2double(string(rows.realInsertion(sel)));  % original insertion number (resets across animals)
    blk.insertionSeq  = str2double(string(rows.insertion(sel)));      % table's sequential insertion index (reference only; != Exp)
    blk.NeurID        = str2double(string(rows.NeurID(sel)));         % good-unit index (needed later for markUnits)
    blk.(matlab.lang.makeValidName(metric)) = m(sel);                % metric value (column named e.g. "Zscore" or "SpkR")
    blk.distToMedian  = abs(m(sel) - target);                        % |metric - median|; 0 = most representative

    picks = [picks; blk];                                            %#ok<AGROW> % append this stimulus's block
end

assert(~isempty(picks), 'No candidate neurons selected for any requested stimulus.'); % nothing to resolve -> stop clearly

% =========================================================================
% 4. Resolve experiment ID (Exp) and Phy cluster ID (PhyID) per pick
% =========================================================================
picks.Exp   = nan(height(picks), 1);                 % globally-unique experiment ID (loadNPclassFromTable key)
picks.PhyID = nan(height(picks), 1);                 % Phy cluster ID

% Target (animal, realInsertion) pairs we must locate among the experiments.
pairKey   = @(a,ri) sprintf('%s|%g', a, ri);         % canonical key for one (animal, realInsertion) pair
wantKeys  = arrayfun(@(k) pairKey(picks.animal(k), picks.realInsertion(k)), ...
                     (1:height(picks))', 'UniformOutput', false);    % key per pick row
uniqueKeys = unique(wantKeys);                        % distinct experiments we need to load
resolved   = containers.Map('KeyType','char','ValueType','any');     % key -> struct(exp, phy) once found

for e = expList(:)'                                  % scan the experiment list
    if resolved.Count == numel(uniqueKeys), break; end % stop early once every needed experiment is found
    try                                              % loading may fail for a stale/missing recording
        NP = loadNPclassFromTable(e);                % Neuropixels object for this experiment
    catch ME                                         % could not load this one
        warning('Exp %g: could not load NP (%s) — skipping.', e, ME.message); % report and continue scanning
        continue
    end

    % --- parse animal / realInsertion exactly as AllExpAnalysis / findFalseNegAndPos do ---
    a = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once')); % Pogona animal token
    if a == "", a = string(regexp(NP.recordingName, 'SA\d+', 'match', 'once')); end % else Laudakia token
    tok    = strsplit(NP.recordingName, '_');         % split the recording name on underscores
    realIns = str2double(tok{end});                   % insertion = trailing recording-name token
    key    = pairKey(a, realIns);                     % this experiment's (animal, realInsertion) key

    if ~ismember(key, uniqueKeys) || isKey(resolved, key) % not needed, or already resolved
        continue                                      % skip experiments that hold none of our picks
    end

    % --- good-unit Phy IDs (order matches NeurID); shared across stimuli of this recording ---
    stimHere = picks.stimulus(find(strcmp(wantKeys, key), 1)); % a stimulus present for this experiment's picks
    try                                                        % object build / sorting read may fail
        obj     = buildStimObj(NP, stimHere);                 % stimulus-specific analysis object
        ps      = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder); % Phy sorting table
        phy_IDg = ps.phy_ID(string(ps.label') == 'good');     % Phy IDs of good units, in NeurID order
        resolved(key) = struct('exp', e, 'phy', phy_IDg(:), 'NP', NP); % cache resolution + NP handle (reused for plotting)
    catch ME                                                   % resolution failed for this experiment
        warning('Exp %g (%s): could not resolve Phy IDs (%s).', e, a, ME.message); % report; picks stay NaN
    end
end

% --- fill Exp / PhyID for every pick from the resolved map ---
for k = 1:height(picks)                              % each chosen neuron
    key = wantKeys{k};                               % its (animal, realInsertion) key
    if ~isKey(resolved, key)                         % its experiment was never resolved
        warning('Pick %s ins %g NeurID %g: experiment not found in expList — Exp/PhyID left NaN.', ...
            picks.animal(k), picks.realInsertion(k), picks.NeurID(k)); % flag the gap
        continue                                     % leave Exp/PhyID as NaN
    end
    r = resolved(key);                               % struct(exp, phy) for this experiment
    assert(picks.NeurID(k) >= 1 && picks.NeurID(k) <= numel(r.phy), ... % ALIGNMENT INVARIANT: NeurID must index the good-unit list
        'NeurID %g out of range (exp %g has %d good units) — good-unit ordering mismatch.', ...
        picks.NeurID(k), r.exp, numel(r.phy));
    picks.Exp(k)   = r.exp;                          % globally-unique experiment ID
    picks.PhyID(k) = r.phy(picks.NeurID(k));         % Phy cluster ID (reverse of plotRaster's ismember)
end

% =========================================================================
% 5. Order, report, return
% =========================================================================
T = sortrows(picks, {'stimulus','distToMedian'});    % group by stimulus, most-representative first

fprintf('\nextractExampleNeurons: %d candidate(s) around the %gth percentile of %s.\n', ...
    height(T), params.percentile, metric);           % header line
for stim = params.stimuli                            % per-stimulus summary
    ts = T(T.stimulus == stim, :);                   % this stimulus's picks
    if isempty(ts), continue; end                    % skip stimuli with no picks
    fprintf('  %s: median %s = %.3f | %d neuron(s): %s\n', stim, metric, ...
        median(ts.(matlab.lang.makeValidName(metric))), height(ts), ...  % median metric across picks
        strjoin(compose('%s-ins%g-exp%g-phy%g', ts.animal, ts.realInsertion, ts.Exp, ts.PhyID), ', ')); % one token per pick
end

% =========================================================================
% 6. (optional) Plot each candidate's raster via its experiment's plotRaster
% =========================================================================
if params.plot                                       % only when plotting was requested
    % Save target: subfolder "neuron candidates_stimA_stimB" under outRoot (paper figs only).
    saveDir = fullfile(params.outRoot, "neuron candidates_" + strjoin(params.stimuli, "_")); % requested subfolder name
    if params.PaperFig && ~exist(saveDir, 'dir'), mkdir(saveDir); end % create only when we will save

    okRows = T(~isnan(T.Exp), :);                     % only picks whose experiment resolved
    gExp = findgroups(okRows.Exp, okRows.stimulus);  % group by experiment x stimulus -> one plotRaster call each
    for g = 1:max(gExp)                               % each (experiment, stimulus) group
        sub  = okRows(gExp == g, :);                  % the picks in this group
        stim = sub.stimulus(1);                       % stimulus for this group
        exp  = sub.Exp(1);                            % experiment for this group
        key  = pairKey(char(sub.animal(1)), sub.realInsertion(1)); % (animal, realInsertion) key -> cached NP
        NPh  = resolved(key).NP;                       % reuse the NP loaded during resolution
        phyIds = sub.PhyID(:)';                        % all Phy IDs for this experiment+stimulus (vector -> one call)

        before = findall(0, 'Type', 'figure');        % figures open before the call
        try
            obj = buildStimObj(NPh, stim);            % this experiment's analysis object
            plotPicks(obj, stim, phyIds);             % draw the raster(s); saving is handled here, not by plotRaster
        catch ME
            warning('Exp %g %s: plotRaster failed (%s) — skipping.', exp, stim, ME.message); % report and continue
            continue
        end
        if ~params.PaperFig, continue; end            % nothing to save unless paper figures requested
        newFigs = setdiff(findall(0, 'Type', 'figure'), before); % figure(s) the call created
        for fi = 1:numel(newFigs)                      % save each (plotRaster opens one per neuron)
            fname = sprintf('%s_%s_ins%g_exp%g.jpg', stim, sub.animal(1), sub.realInsertion(1), exp); % base name
            if numel(newFigs) > 1, fname = strrep(fname, '.jpg', sprintf('_%d.jpg', fi)); end % disambiguate multiples
            exportgraphics(newFigs(fi), fullfile(saveDir, fname), 'Resolution', 300); % write into the subfolder
        end
    end
    if params.PaperFig, fprintf('Saved candidate rasters to:\n  %s\n', saveDir); end % tell the user where they went
end

end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function obj = buildStimObj(NP, stim)
% buildStimObj  Analysis object for a stimulus abbreviation (mirrors plotFalseNegPosNeurons).
    switch stim                                          % dispatch on the stimulus label
        case "RG",   obj = rectGridAnalysis(NP);         % rectified grid (RG/SB)
        case "MB",   obj = linearlyMovingBallAnalysis(NP);% moving ball
        case "MBR",  obj = linearlyMovingBarAnalysis(NP); % moving bar
        case {"SDGs","SDGm"}, obj = StaticDriftingGratingAnalysis(NP); % grating
        case "NI",   obj = imageAnalysis(NP);            % natural image
        case "NV",   obj = movieAnalysis(NP);            % natural movie
        case "FFF",  obj = fullFieldFlashAnalysis(NP);   % full-field flash
        otherwise,   error('Unknown stimulus "%s".', stim); % unrecognised label
    end
end

function plotPicks(obj, stim, phyIds)
% plotPicks  Call the experiment's plotRaster for a set of Phy IDs, per stimulus.
%   Argument sets follow the user's canonical MB / RG example calls, EXCEPT
%   PaperFig=false and overwrite=false: plotRaster's printFig self-saves on EITHER
%   PaperFig=true (to Paper_figs) OR overwrite=true (to visualStimPlotsFolder), and
%   PaperFig does not change the figure's appearance — only where it saves. Forcing
%   both off makes this function the sole writer, so rasters land only in the subfolder.
    switch stim                                          % stimulus-specific plotRaster arguments
        case {"MB","MBR"}                                % moving ball / bar
            obj.plotRaster('exNeuronsPhyID', phyIds, 'overwrite', false, ...
                'MergeNtrials', 3, 'PaperFig', false, 'OneLuminosity', 'white'); % MB example call
        case "RG"                                        % rectified grid (RG/SB)
            obj.plotRaster('MergeNtrials', 1, 'overwrite', false, 'AllResponsiveNeurons', false, ...
                'exNeuronsPhyID', phyIds, 'selectedLum', 255, 'oneTrial', true, 'PaperFig', false); % RG example call
        otherwise                                        % other stimuli: minimal generic call
            obj.plotRaster('exNeuronsPhyID', phyIds, 'overwrite', false, ...
                'AllResponsiveNeurons', false, 'PaperFig', false); % generic fallback
    end
end

function el = pickExpList(S, argExpList)
% pickExpList  Prefer the cache's saved expList, else the caller-supplied one.
    if isfield(S,'expList') && ~isempty(S.expList)   % cache stored the experiment list
        el = S.expList;                              % use it
    else                                             % cache had none
        el = argExpList;                             % fall back to the params argument
    end
end
