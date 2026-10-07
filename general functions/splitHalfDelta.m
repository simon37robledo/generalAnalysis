function dlt = splitHalfDelta(X, grp, half)
% splitHalfDelta  Split-half cross-validated peak response (spikes/s).
%
%   dlt = splitHalfDelta(X, grp, half)
%
%   For each neuron, picks the (group, option) with the largest mean on one
%   half of the trials and measures that same (group, option) on the other,
%   disjoint half; the two directions of the split are averaged. Selecting and
%   measuring on independent trials removes the winner's-curse inflation of a
%   plain max, so a neuron with no response averages ~0 (validated under the
%   null in checkMovingWindowNullFPR).
%
%   INPUTS
%     X     [nTrials x nNeurons x nOptions]  response minus baseline per trial,
%           in spikes/ms; options are e.g. grid cells (NaN = not measurable,
%           such as a cell the ball did not cross on that trial)
%     grp   [nTrials x 1]  selection group per trial (category or direction), 1..nGrp
%     half  [nTrials x 1]  1 or 2: which half of the split the trial belongs to
%
%   OUTPUT
%     dlt   [1 x nNeurons]  held-out mean at the selected (group, option), spikes/s.
%           NaN if the selected option has no trials in the other half.

[nTr, nU, nOpt] = size(X);                                         % dimensions
assert(numel(grp) == nTr && numel(half) == nTr, ...
    'splitHalfDelta: grp and half must have one entry per trial (row of X).');  % alignment invariant
assert(all(ismember(half, [1 2])), 'splitHalfDelta: half must contain only 1 and 2.');  % valid split labels
assert(any(half == 1) && any(half == 2), 'splitHalfDelta: both halves need trials.');   % split must exist
grp  = grp(:);                                                     % column vector
half = half(:);                                                    % column vector
nGrp = max(grp);                                                   % selection groups
ok   = ~isnan(X);                                                  % measurable entries
X(~ok) = 0;                                                        % zero-fill for summation
M = cell(1, 2);                                                    % group x option means per half
for h = 1:2
    G    = double(sparse(grp(half == h), find(half == h), 1, nGrp, nTr));  % [nGrp x nTrials] membership in this half
    sums = full(G * reshape(double(X),  nTr, []));                 % sums   [nGrp x nNeurons*nOpt]
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
