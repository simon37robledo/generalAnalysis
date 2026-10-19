function T = directionTuningUnits(vsObj, stim, params)
% directionTuningUnits  Per-unit direction tuning for Figure 5 (MB, SDGm, MBR).
%
%   T = directionTuningUnits(vsObj, "MB")      % linearlyMovingBallAnalysis
%   T = directionTuningUnits(vsObj, "SDGm")    % StaticDriftingGratingAnalysis
%   T = directionTuningUnits(vsObj, "MBR")     % linearlyMovingBarAnalysis
%
%   One row per good unit. Responses are baseline-subtracted rates (spikes/s)
%   in the same windows as the Figure 1-4 responsiveness tests:
%     MB   - 300 ms window from the ball crossing the unit's RF grid cell
%            (cached spatial-grid statistics, ObsDiff). The cell is chosen on
%            the response AVERAGED OVER DIRECTIONS, so the choice does not
%            favour any direction; speed = the Figure 1 speed (lowest p).
%     SDGm - first SDGWindows(2) ms after motion onset (one period of the
%            slowest temporal frequency), vs the blank pre-trial baseline of
%            the same length capped by the ITI; temporal frequencies pooled.
%     MBR  - whole sweep vs a 300 ms pre-trial baseline (a bar crosses every
%            RF once per sweep in every direction).
%   Directions are in the moving-ball convention (phi: 0 = up, pi/2 = left,
%   pi = down, 3pi/2 = right; VS_linearlyMovingBall/Bar). Grating angles are
%   mapped to their drift direction (VS_StaticDriftingGrating: procedural
%   sine grating, phase increasing -> drift to -x, rotated clockwise):
%   0/360 deg -> left, 90 -> up, 180 -> right, 270 -> down.
%
%   Tuning: mean response per direction, negatives clamped to 0, then
%   DSI = |sum R e^{i theta}| / sum R and OSI at 2 theta (Ringach et al. 2002;
%   Mazurek et al. 2014). pDSI/pOSI: permutation test shuffling direction labels
%   across the unit's trials (same statistic in the null), which also absorbs
%   the small-sample upward bias of DSI.
%
%   Columns: phyID, unitIdx, pResp, responsive, present, R (cell: 1 x nDir, spikes/s),
%   nPerDir (cell), DSI, OSI, prefDir (rad), pDSI, pOSI, cell (MB), speed (MB), dirs (cell).

arguments
    vsObj                                                          % stimulus analysis object
    stim (1,1) string {mustBeMember(stim, ["MB", "SDGm", "MBR"])}  % stimulus family
    params.nPerm (1,1) double = 2000                               % direction-label permutations
    params.seed  (1,1) double = 1                                  % RNG seed (reproducible p-values)
    params.alpha (1,1) double = 0.05                               % responsiveness threshold
    params.minTrialsPerDir (1,1) double = 3                        % MB: trials per direction at the chosen cell
end

rng(params.seed);                                                  % reproducible permutations
Stats = vsObj.StatisticsPerNeuron;                                 % bare call: figure-design cache (errors if stale)
p     = vsObj.dataObj.convertPhySorting2tIc(vsObj.spikeSortingFolder);   % spike sorting
label = string(p.label');                                          % unit labels
goodU = p.ic(:, label == 'good');                                  % good units
phyID = p.phy_ID(label == 'good');                                 % their phy IDs
nU    = size(goodU, 2);                                            % number of good units
% Grating angle a (deg) -> drift direction in MB phi: the unrotated grating drifts
% left (phi = pi/2) and DrawTexture rotates it clockwise by a, while phi grows
% counter-clockwise, so phi = pi/2 - a. Checks: 0/360 -> left, 90 -> up (0),
% 180 -> right (3pi/2), 270 -> down (pi); 45-degree sets follow the same rule.
sdgDrift = @(a) mod(pi/2 - deg2rad(a), 2*pi);                      % angle (deg) -> drift direction (rad)
assert(abs(sdgDrift(90)) < 1e-9 && abs(sdgDrift(180) - 3*pi/2) < 1e-9 && abs(sdgDrift(360) - pi/2) < 1e-9, ...
    'directionTuningUnits: grating direction mapping changed.');   % guard the convention

switch stim
    case "MB"
        spF  = fieldnames(Stats); spF = spF(startsWith(string(spF), "Speed"));   % Speed1, Speed2, ...
        allP = cell2mat(cellfun(@(f) Stats.(f).pvalsResponse(:), spF', 'UniformOutput', false));  % [nU x nSpeeds]
        [pResp, bestSp] = min(allP, [], 2);                        % Figure 1 rule: lowest-p speed per unit
        present = statPresence(Stats, "MB", "");                   % present in every speed's block
    case "SDGm"
        pResp   = Stats.Moving.pvalsResponse(:);                   % drifting-phase p-values
        present = logical(Stats.Moving.present(:));                % present during the grating block
    case "MBR"
        pResp   = Stats.Speed1.pvalsResponse(:);                   % single bar speed
        present = statPresence(Stats, "MBR", "");                  % present during the bar block
end
assert(numel(pResp) == nU && numel(present) == nU, ...
    'directionTuningUnits: statistics cover %d units, spike sorting has %d good units.', numel(pResp), nU);  % unit alignment

% ---- Per-trial responses X and direction labels ----------------------------
if stim == "MB"
    Xu = cell(nU, 1); Lu = cell(nU, 1); cellU = nan(nU, 1);       % per-unit trials, labels, chosen cell
    for u = 1:nU
        g    = Stats.(spF{bestSp(u)});                             % grid results at this unit's speed
        D    = g.ObsDiff(:, :, u) * 1000;                          % [trials x cells] response - baseline, spikes/s
        dirs = g.trialDirections(:);                               % direction of each trial (MB phi)
        uD   = unique(dirs);                                       % directions in this block
        nPer = zeros(numel(uD), size(D, 2));                       % crossings per direction x cell
        mPer = nan(numel(uD), size(D, 2));                         % mean response per direction x cell
        for k = 1:numel(uD)
            Dk = D(dirs == uD(k), :);                              % this direction's trials
            nPer(k, :) = sum(~isnan(Dk), 1);                       % trials that crossed each cell
            mPer(k, :) = mean(Dk, 1, 'omitnan');                   % mean response per cell
        end
        okCell = all(nPer >= params.minTrialsPerDir, 1);           % cells crossed in every direction
        score  = mean(mPer, 1); score(~okCell) = -Inf;             % direction-AVERAGED response (no direction favoured)
        [best, c] = max(score);                                    % RF cell
        if ~isfinite(best), continue, end                          % no cell crossed in every direction
        keep = ~isnan(D(:, c));                                    % trials that crossed the RF cell
        Xu{u} = D(keep, c); Lu{u} = dirs(keep); cellU(u) = c;      % responses and labels at the RF cell
    end
else
    rw = vsObj.ResponseWindow;                                     % condition matrix and timing
    cn = rw.colNames{1}(5:end);                                    % condition column names
    if stim == "SDGm"
        C  = rw.C;                                                 % [onset angle tf sf]
        ci = find(strcmpi(cn, 'angles'), 1) + 1;                   % angle column
        on = C(:, 1) + vsObj.VST.static_time * 1000;               % motion onsets (ms)
        W  = Stats.params.BaseRespWindow;                          % default window
        if isfield(Stats.params, 'SDGWindows') && ~isempty(Stats.params.SDGWindows)
            W = Stats.params.SDGWindows(2);                        % drifting window of the figure design
        end
        W  = min(W, floor(rw.Moving.stimDur));                     % never past the drifting phase
        bl = min(W, rw.stimInter - 100);                           % blank baseline length (ms), as the stats
        baseOn = C(:, 1) - bl;                                     % baseline before the static onset
        lab = sdgDrift(C(:, ci));                                  % drift direction per trial (MB phi)
        assert(~any(isnan(lab)), 'directionTuningUnits: non-numeric grating angles.');   % every trial labelled
    else
        C  = rw.Speed1.C;                                          % [onset direction ...]
        ci = find(strcmpi(cn, 'directions'), 1) + 1;               % direction column
        on = C(:, 1);                                              % sweep onsets (ms)
        W  = floor(rw.Speed1.stimDur);                             % whole sweep (ms)
        bl = min(300, rw.stimInter - 100);                         % pre-trial baseline (ms)
        baseOn = on - bl;                                          % baseline start
        lab = mod(C(:, ci), 2*pi);                                 % already MB phi (radians)
    end
    Mr = BuildBurstMatrix(goodU, round(p.t), round(on(:)'), W);    % [trials x units x ms] response window
    Mb = BuildBurstMatrix(goodU, round(p.t), round(baseOn(:)'), bl);   % [trials x units x ms] baseline
    X  = (mean(Mr, 3) - mean(Mb, 3)) * 1000;                       % [trials x units] response - baseline, spikes/s
    Xu = num2cell(X, 1)'; Lu = repmat({lab}, nU, 1);               % same trials for every unit
end

% ---- Tuning, selectivity and permutation p-values ---------------------------
allDirs = uniquetol(cell2mat(Lu(~cellfun(@isempty, Lu))), 1e-4);   % directions present
nD = numel(allDirs);                                               % usually 4
R = nan(nU, nD); nPerDir = zeros(nU, nD);                          % tuning curves, trial counts
DSI = nan(nU, 1); OSI = nan(nU, 1); prefDir = nan(nU, 1); pDSI = nan(nU, 1); pOSI = nan(nU, 1);
nullDSImean = nan(nU, 1);                                          % mean DSI under shuffled direction labels
for u = 1:nU
    x = Xu{u}; l = Lu{u};                                          % this unit's trials and labels
    if isempty(x), continue, end                                   % MB unit without a usable RF cell
    [~, li] = ismembertol(l, allDirs, 1e-4);                       % direction index per trial
    assert(all(li > 0), 'directionTuningUnits: unknown direction label.');   % every trial labelled
    [R(u, :), nPerDir(u, :)] = dirMeans(x, li, nD);                % mean response per direction
    [DSI(u), OSI(u), prefDir(u)] = tuneIdx(R(u, :), allDirs);      % observed selectivity
    nullD = zeros(params.nPerm, 1); nullO = zeros(params.nPerm, 1);   % permutation nulls
    for b = 1:params.nPerm
        Rb = dirMeans(x, li(randperm(numel(li))), nD);             % shuffled direction labels
        [nullD(b), nullO(b)] = tuneIdx(Rb, allDirs);               % same statistic under H0
    end
    pDSI(u) = (sum(nullD >= DSI(u)) + 1) / (params.nPerm + 1);     % permutation p (DSI)
    nullDSImean(u) = mean(nullD);                                  % chance-level DSI for this unit's trial counts
    pOSI(u) = (sum(nullO >= OSI(u)) + 1) / (params.nPerm + 1);     % permutation p (OSI)
end

T = table(phyID(:), (1:nU)', pResp, pResp < params.alpha, present, num2cell(R, 2), num2cell(nPerDir, 2), DSI, OSI, prefDir, pDSI, pOSI, ...
    'VariableNames', {'phyID','unitIdx','pResp','responsive','present','R','nPerDir','DSI','OSI','prefDir','pDSI','pOSI'});
T.stim = repmat(stim, nU, 1);                                      % stimulus family
T.cell = nan(nU, 1); T.speed = nan(nU, 1);                         % MB-only columns
if stim == "MB", T.cell = cellU; T.speed = bestSp; end             % RF cell and Figure 1 speed
T.dirs = repmat({allDirs(:)'}, nU, 1);                             % directions of each tuning curve (MB phi); cells: 4 or 8 directions per recording
T.nullDSImean = nullDSImean;                                       % DSI expected by chance (small-sample bias)
T.DSIcorr = DSI - nullDSImean;                                     % bias-corrected DSI: 0 = chance
end

% =========================================================================
function [m, n] = dirMeans(x, li, nD)
% Mean response and trial count per direction index.
m = accumarray(li(:), x(:), [nD 1], @mean, NaN)';                  % [1 x nD]
n = accumarray(li(:), 1, [nD 1])';                                 % [1 x nD]
end

function [dsi, osi, pref] = tuneIdx(R, th)
% Vector-sum selectivity on a tuning curve with negatives clamped to 0.
R = max(R, 0); s = sum(R, 'omitnan');                              % rectified responses
if s == 0, dsi = 0; osi = 0; pref = NaN; return, end               % no response in any direction
v1 = sum(R(:) .* exp(1i * th(:)), 'omitnan');                      % resultant at 1 theta
v2 = sum(R(:) .* exp(2i * th(:)), 'omitnan');                      % resultant at 2 theta
dsi = abs(v1) / s; osi = abs(v2) / s; pref = mod(angle(v1), 2*pi); % DSI, OSI, preferred direction
end
