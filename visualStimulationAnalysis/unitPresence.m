function [present, ratio] = unitPresence(tMs, icGood, blockStart, blockEnd, nBins, minRatio)
% unitPresence  Per-block presence ratio of each unit (BombCell / Allen metric, per block).
%
%   [present, ratio] = unitPresence(tMs, icGood, blockStart, blockEnd)
%   [present, ratio] = unitPresence(tMs, icGood, blockStart, blockEnd, nBins, minRatio)
%
%   Units drift in and out within a session, so a unit sorted as "good" over the
%   whole recording may be absent while a given stimulus runs. The block
%   [blockStart, blockEnd] is split into nBins equal bins; the presence ratio is
%   the fraction of bins in which the unit fires at least once (as BombCell's
%   presenceRatio, computed within the block instead of over the whole recording).
%   A unit is present if ratio >= minRatio.
%
%   INPUTS
%     tMs        spike times (ms), p.t of convertPhySorting2tIc
%     icGood     [4 x nUnits] ic columns of the good units; rows 3 and 4 are the
%                first and last index of each unit's spikes in tMs
%     blockStart, blockEnd   block limits (ms): first onset, last onset + stimDur
%     nBins      bins per block (default 10: blocks can be only a few minutes long)
%     minRatio   minimum presence ratio (default 0.7, BombCell's minPresenceRatio)
%
%   OUTPUTS
%     present    [1 x nUnits] logical
%     ratio      [1 x nUnits] presence ratio in [0, 1]

if nargin < 5 || isempty(nBins),    nBins    = 10;  end           % default number of bins
if nargin < 6 || isempty(minRatio), minRatio = 0.7; end           % default BombCell threshold
assert(blockEnd > blockStart, 'unitPresence: block end must follow block start.');  % valid block
assert(size(icGood, 1) >= 4, 'unitPresence: icGood must have rows 3-4 = spike index range.');  % ic layout

edges = linspace(blockStart, blockEnd, nBins + 1);                % bin edges spanning the block (ms)
nU    = size(icGood, 2);                                           % units
ratio = zeros(1, nU);                                              % presence ratio per unit
for k = 1:nU
    st = tMs(icGood(3, k):icGood(4, k));                           % this unit's spike times
    st = st(st >= blockStart & st <= blockEnd);                    % keep spikes inside the block
    ratio(k) = mean(histcounts(st, edges) > 0);                    % fraction of bins with >= 1 spike
end
present = ratio >= minRatio;                                       % present throughout the block
end
