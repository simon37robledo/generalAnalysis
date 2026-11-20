function R = eyeResponseCorrelation(vsObj, Tr, neurIDs, stim, params)
% EYERESPONSECORRELATION  Is a unit's trial-to-trial response related to the eye movement in
% that trial?
%
%   R = eyeResponseCorrelation(vsObj, Tr, neurIDs, stim)
%     vsObj    stimulus object (linearlyMovingBallAnalysis or rectGridAnalysis)
%     Tr       per-trial eye metrics of the same block (eyeTrialMetrics, presentation order)
%     neurIDs  unit indices among the recording's good units (TableStimComp NeurID)
%     stim     "MB" or "SB"
%
%   Per trial: response = spike rate in the response window minus the rate in the 300 ms before
%   onset. Window: MB = the whole trial (onset -> offset; the ball crosses the RF somewhere in
%   it), SB = 0-300 ms after onset (the SB test window). Eye measures over the same window:
%   max displacement (movement) and |mean position - block median| (where the eye was).
%   Both response and eye measures are centred within each stimulus condition first (MB:
%   direction x offset x size x speed x luminosity; SB: position x tiling x luminosity), so ball
%   position relative to the RF cannot create a correlation. Then Spearman rho per unit over the
%   trials with >= 80% valid eye frames. Slow drift over the session (in firing and in eye
%   position) would correlate the two without any trial-by-trial link, so both are also detrended
%   with a running median over detrendTrials trials in presentation order (primary; 0 = raw).
%   If MB responses were driven by eye movements, rho(response, movement) would be > 0.

arguments
    vsObj                                                                          % stimulus object
    Tr table                                                                       % eyeTrialMetrics output
    neurIDs double                                                                 % units to test
    stim (1,1) string {mustBeMember(stim, ["MB" "SB"])}                            % which stimulus
    params.baseMs double = 300                                                     % pre-onset baseline [ms]
    params.sbWindowMs double = 300                                                 % SB response window [ms]
    params.binMs double = 10                                                       % spike-count bin [ms]
    params.minTrials double = 20                                                   % fewest usable trials for a correlation
    params.detrendTrials double = 51                                               % running-median window over trials removed from both (0 = none)
end

p = vsObj.dataObj.convertPhySorting2tIc(vsObj.spikeSortingFolder, 0, 1, 1);        % spike times and units (as ResponseWindow)
label = string(p.label');                                                          % unit labels
goodU = p.ic(:, label == 'good');                                                  % good units, the order NeurID indexes
NR = vsObj.ResponseWindow;                                                         % cached response window (same unit set as the stats)
if isfield(NR, 'Speed1'), NR = NR.Speed1; end                                     % MB stores one struct per speed (same units)
assert(size(NR.NeuronVals, 1) == size(goodU, 2), 'eyeResponseCorrelation:units', ...
    '%s: %d good units now, %d in ResponseWindow (labels changed?).', stim, size(goodU, 2), size(NR.NeuronVals, 1));   % unit order invariant
assert(all(neurIDs >= 1 & neurIDs <= size(goodU, 2)), 'eyeResponseCorrelation:neurID', 'NeurID outside 1..%d.', size(goodU, 2));

S = vsObj.getSyncedDiodeTriggers;                                                  % onsets/offsets [ms], presentation order
on = S.stimOnFlipTimes(:); off = S.stimOffFlipTimes(:);                            % columns
assert(numel(on) == height(Tr) && all(abs(on - Tr.onMs) < 1e-6), 'eyeResponseCorrelation:trials', ...
    '%s: trial onsets differ from the eye table.', stim);                         % same trials, same order
if stim == "MB"
    winMs = off - on;                                                              % whole trial
    eyeDisp = Tr.stim_dispPx; keep = Tr.stim_kept;                                 % eye movement over the trial
    eyeOff = hypot(Tr.stim_meanXPx, Tr.stim_meanYPx);                              % eye position offset in the trial
    V = vsObj.VST;                                                                 % trial parameters
    lum = 255 * ones(numel(on), 1); if isfield(V, 'Luminosities'), lum = V.Luminosities(:); end   % luminosity (absent = white)
    cond = findgroups(V.directions(:), V.offsets(:), V.ballSizes(:), V.speeds(:), lum);   % condition per trial
else
    winMs = repmat(params.sbWindowMs, numel(on), 1);                               % SB test window
    eyeDisp = Tr.post_dispPx; keep = Tr.post_kept;                                 % eye movement in 0-300 ms
    eyeOff = nan(numel(on), 1);                                                    % mean position: computed over the whole trial only
    if ismember('stim_meanXPx', Tr.Properties.VariableNames), eyeOff = hypot(Tr.stim_meanXPx, Tr.stim_meanYPx); end
    V = vsObj.VST;                                                                 % trial parameters
    tiling = ones(numel(on), 1); if isfield(V, 'tilingRatios'), tiling = V.tilingRatios(:); end   % tiling (absent = 1)
    cond = findgroups(V.pos(:), tiling, V.luminosities(:));                        % condition per trial
end

b = params.binMs;                                                                  % bin [ms]
nWin = round(max(winMs) / b);                                                      % bins covering the longest window
Mr = BuildBurstMatrix(goodU(:, neurIDs), round(p.t / b), round(on' / b), nWin);    % [trials x units x bins] spike counts (onsets as a ROW: a column gives 1 trial)
Mb = BuildBurstMatrix(goodU(:, neurIDs), round(p.t / b), round((on' - params.baseMs) / b), round(params.baseMs / b));   % baseline counts
assert(size(Mr, 1) == numel(on) && size(Mb, 1) == numel(on), 'eyeResponseCorrelation:burst', 'Spike matrix has %d trials, expected %d.', size(Mr, 1), numel(on));   % trial count invariant
nBinT = round(winMs / b);                                                          % bins inside each trial's own window
inWin = (1:nWin) <= nBinT;                                                         % [trials x bins] mask
cnt = squeeze(sum(Mr .* permute(inWin, [1 3 2]), 3));                             % spikes in each trial's window [trials x units]
if isvector(cnt) && numel(neurIDs) == 1, cnt = cnt(:); end                         % one unit: keep a column
base = squeeze(sum(Mb, 3)); if isvector(base) && numel(neurIDs) == 1, base = base(:); end   % baseline spikes
resp = cnt ./ (winMs / 1000) - base / (params.baseMs / 1000);                      % response [spk/s], baseline-subtracted

ok = keep & isfinite(eyeDisp);                                                     % trials with usable eye data
nU = numel(neurIDs);                                                               % units
R = table(neurIDs(:), repmat(stim, nU, 1), nan(nU, 1), nan(nU, 1), repmat(nnz(ok), nU, 1), ...
    'VariableNames', {'NeurID', 'stim', 'rhoMovement', 'rhoPosition', 'nTrials'});    % output
if nnz(ok) < params.minTrials, return, end                                         % too few trials: NaN (counted by the caller)
condOk = findgroups(cond(ok));                                                     % conditions of usable trials, renumbered 1..k
centre = @(x) x - subsref(accumarray(condOk, x, [], @mean), struct('type', '()', 'subs', {{condOk}}));   % x minus its condition mean
detr = @(x) x - (params.detrendTrials > 1) * movmedian(x, max(params.detrendTrials, 1));   % remove slow drift over trials (presentation order)
e1 = detr(centre(eyeDisp(ok)));                                                    % movement, condition-centred, detrended
e2 = eyeOff(ok);                                                                   % position offset of usable trials
if all(isfinite(e2)), e2 = detr(centre(e2)); end                                   % condition-centred, detrended
for u = 1:nU
    r = detr(centre(resp(ok, u)));                                                 % this unit's responses, condition-centred, detrended
    if std(r) == 0, continue, end                                                  % no spikes / constant: undefined
    R.rhoMovement(u) = corr(r, e1, 'Type', 'Spearman');                            % response vs eye movement
    if all(isfinite(e2)), R.rhoPosition(u) = corr(r, e2, 'Type', 'Spearman'); end  % response vs eye position offset
end
end
