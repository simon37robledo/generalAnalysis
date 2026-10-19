function bt = hierBootStatistic(statFun, nrep, lvl1, lvl2)
% hierBootStatistic  Hierarchical bootstrap of an arbitrary statistic.
%
%   bt = hierBootStatistic(statFun, nrep, lvl1, lvl2)
%
%   Same 3-level resampling as hierBootMatchFreq (Saravanan et al. 2020):
%   animals (lvl1) with replacement, then insertions (lvl2) within each drawn
%   animal, then observations within each drawn insertion, matching the
%   observed counts. Instead of a mean it returns statFun(idx) for every
%   resample, where idx are the drawn ROW INDICES, so statistics that need
%   several paired columns (correlations, within-insertion ranks) stay paired.
%   A drawn insertion keeps its identity, so the resampled rows of one
%   insertion can still be grouped by lvl2(idx).
%
%   statFun - handle: idx (column of row indices) -> scalar
%   nrep    - number of resamples
%   lvl1    - animal ID per row
%   lvl2    - insertion ID per row, unique ACROSS animals (as hierBootMatchFreq)
%
%   See also: hierBootMatchFreq, hierBootCorrSpearman (analyzeStripeNeurons)

lvl1 = lvl1(:); lvl2 = lvl2(:);                                    % column vectors
assert(numel(lvl1) == numel(lvl2), 'hierBootStatistic: lvl1 and lvl2 differ in length.');   % one label per row
for g = unique(lvl2)'                                              % insertion IDs must not span animals
    assert(numel(unique(lvl1(lvl2 == g))) == 1, ...
        'hierBootStatistic: insertion %g belongs to more than one animal (IDs must be unique across animals).', g);
end
grp1 = unique(lvl1);                                               % animals
pools = cell(numel(grp1), 1);                                      % per animal: cell of row-index sets per insertion
for k = 1:numel(grp1)
    g2 = unique(lvl2(lvl1 == grp1(k)));                            % insertions of this animal
    pools{k} = arrayfun(@(g) find(lvl2 == g), g2, 'UniformOutput', false);   % row indices per insertion
end
bt = nan(nrep, 1);                                                 % statistic per resample
for j = 1:nrep
    idx = [];                                                      % drawn rows
    for k = randi(numel(grp1), numel(grp1), 1)'                    % animals with replacement
        P = pools{k};                                              % this animal's insertions
        for m = randi(numel(P), numel(P), 1)'                      % insertions with replacement
            pool = P{m};                                           % rows of this insertion
            idx = [idx; pool(randi(numel(pool), numel(pool), 1))]; %#ok<AGROW> rows with replacement
        end
    end
    bt(j) = statFun(idx);                                          % statistic on the resample
end
end
