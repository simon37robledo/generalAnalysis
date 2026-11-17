function [pVal, bootDiff] = bootstrapTwoSampleDifference(tbl, pair, nBoot, metric)
% bootstrapTwoSampleDifference  Unpaired hierarchical bootstrap of mean(A) - mean(B).
%
%   [pVal, bootDiff] = bootstrapTwoSampleDifference(tbl, pair, nBoot, metric)
%
%   For AllExpAnalysis presence="perStimulus" (2026-11-18): each stimulus keeps the
%   units present during its own block, so A and B are different (overlapping) unit
%   sets and cannot be paired. The two groups share animals and insertions, so they
%   are resampled JOINTLY: animals with replacement, then insertions within each drawn
%   animal with replacement, then units of group A and of group B within each drawn
%   insertion, each with replacement and at its own count (as hierBootMatchFreq does
%   for one group). Statistic: pooled mean of A minus pooled mean of B per replicate.
%
%   INPUTS
%     tbl    table with columns animal, insertion, stimulus and the metric column
%            (insertion labels must be unique across animals: keyed here on the
%            (animal, insertion) pair).
%     pair   {A, B} stimulus labels.
%     nBoot  number of bootstrap replicates.
%     metric column name, e.g. 'Z-score' or 'DeltaRate'.
%   OUTPUTS
%     pVal     two-tailed p: 2 x the smaller tail of bootDiff around 0, capped at 1
%              (same convention as bootstrapPairDifference).
%     bootDiff 1 x nBoot bootstrap distribution of mean(A) - mean(B).

isA = tbl.stimulus == pair{1};                                % rows of group A
isB = tbl.stimulus == pair{2};                                % rows of group B
assert(any(isA) && any(isB), 'bootstrapTwoSampleDifference:emptyGroup', ...
    'Group %s or %s has no rows.', char(pair{1}), char(pair{2}));   % both groups must exist
val = tbl.(metric);                                           % metric values
assert(~any(isnan(val(isA | isB))), 'bootstrapTwoSampleDifference:nan', ...
    '%s has NaN values in the compared groups.', metric);     % fail loudly, never drop silently

[gAI, aniKey, insKey] = findgroups(string(tbl.animal), string(tbl.insertion));   % one id per (animal, insertion)
[aniList, ~, aniOfIns] = unique(aniKey);                      % animal of every (animal, insertion) group
nAni = numel(aniList);                                        % number of animals
insOfAni = arrayfun(@(a) find(aniOfIns == a), 1:nAni, 'UniformOutput', false);   % insertion groups per animal
nIns = numel(insKey);                                         % number of (animal, insertion) groups
poolA = cell(nIns, 1); poolB = cell(nIns, 1);                 % unit values per insertion and group
for g = 1:nIns
    poolA{g} = val(isA & gAI == g);                           % group A units of this insertion
    poolB{g} = val(isB & gAI == g);                           % group B units of this insertion
end

bootDiff = nan(1, nBoot);                                     % bootstrap distribution
for r = 1:nBoot
    sA = 0; nA = 0; sB = 0; nB = 0;                           % running sums and counts
    aDraw = randi(nAni, nAni, 1);                             % animals with replacement
    for a = aDraw'
        ins = insOfAni{a};                                    % insertions of this animal
        iDraw = ins(randi(numel(ins), numel(ins), 1));        % insertions with replacement
        for g = iDraw(:)'
            pa = poolA{g}; pb = poolB{g};                     % this insertion's units per group
            if ~isempty(pa)
                x = pa(randi(numel(pa), numel(pa), 1));       % group A units with replacement
                sA = sA + sum(x); nA = nA + numel(x);
            end
            if ~isempty(pb)
                y = pb(randi(numel(pb), numel(pb), 1));       % group B units with replacement
                sB = sB + sum(y); nB = nB + numel(y);
            end
        end
    end
    if nA > 0 && nB > 0
        bootDiff(r) = sA / nA - sB / nB;                      % pooled mean difference
    end
end
nEmpty = sum(isnan(bootDiff));                                % replicates with an empty group
if nEmpty > 0
    warning('bootstrapTwoSampleDifference:emptyReplicates', ...
        '%d of %d replicates drew no unit of one group; they are excluded from p.', nEmpty, nBoot);
end
b = bootDiff(~isnan(bootDiff));                               % valid replicates
pVal = min(1, 2 * min(mean(b <= 0), mean(b >= 0)));           % two-tailed, as bootstrapPairDifference
end
