function verifySplitHalfStore(sub, tag)
% verifySplitHalfStore  Hard check that a cache's stored split reproduces DeltaRateCV.
%
%   verifySplitHalfStore(sub, tag)
%
%   Recomputes splitHalfDelta from the stored trial-level differences and the
%   stored groups/halves, and errors unless both the delta rate and the
%   selected (group, option) per half match the cache exactly. Run it before
%   plotting split-half trials (splitHalfMR) from a cache, so a misaligned
%   trial order can never be displayed silently.
%
%   INPUTS
%     sub   stats sub-struct (as passed to splitHalfMR)
%     tag   text for the error message (e.g. 'rec 74 RG')

if isfield(sub, 'prefGridCell')                                    % MB spatial grid
    [~, ~, dirID] = unique(sub.trialDirections);                   % direction index per trial (selection group)
    [dlt, sel] = splitHalfDelta(permute(sub.ObsDiff, [1 3 2]), dirID, sub.cvHalf);  % recompute from the cache
    selStored = sub.cvSel;                                         % stored selection
else                                                               % onset-tested stimuli
    ci = sub.cvInfo;                                               % stored split
    if isempty(sub.ObsDiffSeg)
        X = sub.ObsDiff;                                           % single window
    else
        X = sub.ObsDiffSeg;                                        % category x segment
    end
    nCats = size(sub.validCats, 1);                                % categories
    trialsCat = size(X, 1) / nCats;                                % repetitions per category
    assert(trialsCat == round(trialsCat), 'verifySplitHalfStore: %s rows not divisible by categories.', tag);
    invalid = repmat(~repelem(sub.validCats, trialsCat, 1), 1, 1, size(X, 3));  % invalid (category, unit), every segment
    X(invalid) = NaN;                                              % same exclusion as StatisticsPerNeuron
    [dlt, sel] = splitHalfDelta(X, ci.cvGrp, ci.cvHalf);           % recompute from the cache
    selStored = ci.cvSel;                                          % stored selection
end
ok = ~isnan(sub.DeltaRateCV);                                      % units with a defined delta rate
assert(isequaln(dlt, sub.DeltaRateCV), 'verifySplitHalfStore:mismatch', ...
    '%s: recomputed DeltaRateCV differs from the cache (max |diff| %.3g).', tag, ...
    max(abs(dlt(ok) - sub.DeltaRateCV(ok))));                      % trial order or grouping misaligned
assert(isequal(sel, selStored), 'verifySplitHalfStore:mismatch', ...
    '%s: recomputed split-half selection differs from the cache.', tag);  % selection must match
end
