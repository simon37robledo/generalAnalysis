function R = stillTrialDelta(vsObj, Tr, neurIDs, stim, params)
% STILLTRIALDELTA  Split-half cross-validated peak response of each unit, from all trials and from
% the trials in which the eye moved least.
%
%   R = stillTrialDelta(vsObj, Tr, neurIDs, stim)
%     vsObj, Tr  stimulus object and its per-trial eye metrics (eyeTrialMetrics, presentation order)
%     neurIDs    unit indices among the recording's good units (TableStimComp NeurID)
%     stim       "MB" or "SB"
%
%   Response of a trial in window w = spike rate in w minus the rate in the 300 ms before onset.
%   Windows: SB = 0-300 ms after onset (the SB test window); MB = consecutive 300 ms segments
%   covering the trial (the ball crosses the RF in one of them; the stats use 300 ms windows).
%   Split half (as DeltaRateCV): repetitions alternate within each condition; the best
%   condition x window is chosen on one half and measured on the other, and vice versa; the two
%   measurements are averaged. Measured three times on the held-out half:
%     dAll    all held-out trials;
%     dStill  held-out trials whose eye displacement is <= the median of their condition
%             (computed over all usable trials of that condition);
%     dMove   held-out trials above that median.
%   Eye displacement: MB = maximum over the trial (stim_dispPx); SB = over 0-300 ms (post_dispPx).
%   Selection never sees the measured trials, so dStill vs dAll is free of selection bias.
%   This is not the Figure 1 statistic (grid-cell crossing windows); it asks whether the response
%   size depends on excluding the trials with the most eye movement.

arguments
    vsObj                                                                          % stimulus object
    Tr table                                                                       % eyeTrialMetrics output
    neurIDs double                                                                 % units
    stim (1,1) string {mustBeMember(stim, ["MB" "SB"])}                            % stimulus
    params.winMs double = 300                                                      % response window / segment length [ms]
    params.baseMs double = 300                                                     % pre-onset baseline [ms]
    params.binMs double = 10                                                       % spike-count bin [ms]
end

p = vsObj.dataObj.convertPhySorting2tIc(vsObj.spikeSortingFolder, 0, 1, 1);        % spikes (as ResponseWindow)
goodU = p.ic(:, string(p.label') == 'good');                                       % good units (NeurID order)
assert(all(neurIDs >= 1 & neurIDs <= size(goodU, 2)), 'stillTrialDelta:neurID', 'NeurID outside 1..%d.', size(goodU, 2));
S = vsObj.getSyncedDiodeTriggers; on = S.stimOnFlipTimes(:); off = S.stimOffFlipTimes(:);   % trial times [ms]
assert(numel(on) == height(Tr) && all(abs(on - Tr.onMs) < 1e-6), 'stillTrialDelta:trials', '%s: trials differ from the eye table.', stim);
V = vsObj.VST;                                                                     % trial parameters
if stim == "MB"
    lum = 255 * ones(numel(on), 1); if isfield(V, 'Luminosities'), lum = V.Luminosities(:); end
    cond = findgroups(V.directions(:), V.offsets(:), V.ballSizes(:), V.speeds(:), lum);   % condition per trial
    nSeg = floor(min(off - on) / params.winMs);                                    % whole segments inside every trial
    eye = Tr.stim_dispPx; keep = Tr.stim_kept;                                     % eye movement over the trial
else
    tiling = ones(numel(on), 1); if isfield(V, 'tilingRatios'), tiling = V.tilingRatios(:); end
    cond = findgroups(V.pos(:), tiling, V.luminosities(:));                        % condition per trial
    nSeg = 1;                                                                      % one window
    eye = Tr.post_dispPx; keep = Tr.post_kept;                                     % eye movement in 0-300 ms
end
assert(nSeg >= 1, 'stillTrialDelta:short', '%s: trials shorter than one window.', stim);

b = params.binMs; nb = round(params.winMs / b);                                    % bins per window
M = BuildBurstMatrix(goodU(:, neurIDs), round(p.t / b), round(on' / b), nSeg * nb);   % [trials x units x bins] (onsets as a row)
Mb = BuildBurstMatrix(goodU(:, neurIDs), round(p.t / b), round((on' - params.baseMs) / b), round(params.baseMs / b));   % baseline
assert(size(M, 1) == numel(on) && size(Mb, 1) == numel(on), 'stillTrialDelta:burst', 'Spike matrix trial count mismatch.');
nT = numel(on); nU = numel(neurIDs);                                               % trials, units
seg = reshape(M, nT, nU, nb, nSeg);                                                % [trials x units x bins x segments] (bins fastest within a segment)
rate = squeeze(sum(seg, 3)) / (params.winMs / 1000);                               % [trials x units x segments] spk/s
rate = reshape(rate, nT, nU, nSeg);                                                % keep 3-D for one unit / one segment
baseRate = reshape(sum(Mb, 3), nT, nU) / (params.baseMs / 1000);                   % [trials x units] spk/s
resp = rate - baseRate;                                                            % baseline-subtracted response per segment

% Repetition index within condition (presentation order) -> halves; still / moving split per condition
rep = zeros(nT, 1); still = false(nT, 1);                                          % repetition number, still flag
for c = unique(cond)'
    k = find(cond == c);                                                           % this condition's trials, in time order
    rep(k) = 1:numel(k);                                                           % repetition number
    ok = k(keep(k) & isfinite(eye(k)));                                            % usable eye data
    still(ok) = eye(ok) <= median(eye(ok));                                        % below-median movement within the condition
end
usable = keep & isfinite(eye);                                                     % trials with eye data
half = 2 - mod(rep, 2);                                                            % 1 = odd repetitions, 2 = even

nC = max(cond);                                                                    % conditions
R = table(neurIDs(:), repmat(stim, nU, 1), nan(nU, 1), nan(nU, 1), nan(nU, 1), ...
    'VariableNames', {'NeurID', 'stim', 'dAll', 'dStill', 'dMove'});               % output
for u = 1:nU
    out = nan(2, 3);                                                               % [selection half x {all, still, move}]
    for h = 1:2
        sel = half == h;                                                           % selection half
        mea = half ~= h & usable;                                                  % measurement half (usable eye data)
        mSel = nan(nC, nSeg);                                                      % mean response per condition x segment, selection half
        for c = 1:nC, mSel(c, :) = mean(reshape(resp(sel & cond == c, u, :), [], nSeg), 1, 'omitnan'); end
        [~, best] = max(mSel(:));                                                  % best condition x segment
        [bc, bs] = ind2sub(size(mSel), best);                                      % its condition and segment
        r = resp(:, u, bs);                                                        % responses in the chosen segment
        out(h, 1) = mean(r(mea & cond == bc), 'omitnan');                          % all held-out trials of the best condition
        out(h, 2) = mean(r(mea & cond == bc & still), 'omitnan');                  % still ones
        out(h, 3) = mean(r(mea & cond == bc & ~still), 'omitnan');                 % moving ones
    end
    R{u, {'dAll', 'dStill', 'dMove'}} = mean(out, 1, 'omitnan');                   % average over the two selection halves
end
end
