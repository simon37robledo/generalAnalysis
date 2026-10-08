function MR = eventAlignedMR(spikeCol, tMs, trialOnsets, eventTimes, preBase, postLen)
% eventAlignedMR  Per-trial spikes: pre-trial baseline followed by the period after an event.
%
%   MR = eventAlignedMR(spikeCol, tMs, trialOnsets, eventTimes, preBase, postLen)
%
%   Row k is [onset_k - preBase, onset_k) concatenated with
%   [event_k, event_k + postLen), so column preBase+1 is the event and every
%   column before it is true (pre-stimulus) baseline. Used when the response is
%   measured relative to an event inside the trial (ball crossing a grid cell,
%   stimulus offset), so the baseline is never taken from the stimulus period.
%
%   INPUTS
%     spikeCol     goodU(:, u) column of the spike-sorting ic matrix
%     tMs          spike times p_sort.t (ms)
%     trialOnsets  [nTrials] stimulus onsets (ms)
%     eventTimes   [nTrials] event times (ms), same trial order
%     preBase      baseline length (ms)
%     postLen      length after the event (ms)
%
%   OUTPUT
%     MR           [nTrials x (preBase + postLen)] spike counts per 1 ms bin
%
%   See also: gridAlignedMR, plotRaster_MultiExp, plotPSTH_MultiExp

nTr = numel(trialOnsets);                                          % trials
assert(numel(eventTimes) == nTr, 'eventAlignedMR: one event per trial required.');      % alignment invariant
assert(all(eventTimes(:) >= trialOnsets(:)), 'eventAlignedMR: an event precedes its trial onset.');  % baseline would overlap
Mb = BuildBurstMatrix(spikeCol, round(tMs), round(trialOnsets(:)' - preBase), round(preBase));  % pre-trial baseline
Me = BuildBurstMatrix(spikeCol, round(tMs), round(eventTimes(:)'), round(postLen));             % from the event
MR = [reshape(Mb, nTr, []), reshape(Me, nTr, [])];                 % [nTrials x (preBase + postLen)]
end
