function R = depthCorrTest(depth, y, ins, animal, params)
% depthCorrTest  Within-insertion rank correlation between unit depth and a measure.
%
%   R = depthCorrTest(depth, y, ins, animal)
%
%   Statistic (rhoW): depth and y are ranked WITHIN each insertion, the ranks are
%   centred per insertion, and the centred ranks are correlated across all units
%   (insertions with < 3 units are dropped). Only within-insertion association
%   counts, so insertions that differ in depth coverage and in response level
%   cannot create a correlation.
%   p   : depths permuted within each insertion (two-tailed, nPerm permutations).
%   CI  : hierarchical bootstrap (animal -> insertion -> unit) of the same statistic.
%   y may be logical (e.g. responsive yes/no): ranks handle ties.
%
%   Output struct: rho, ciLo, ciHi, p, n, nIns, nAnimals.

arguments
    depth  double                                                  % unit depth (um)
    y      double                                                  % measure per unit
    ins    double                                                  % insertion ID per unit (unique across animals)
    animal                                                         % animal ID per unit (string/categorical/numeric)
    params.nPerm (1,1) double = 10000                              % within-insertion permutations
    params.nBoot (1,1) double = 10000                              % hierarchical bootstrap resamples
    params.minPerIns (1,1) double = 3                              % units per insertion to be used
end

depth = depth(:); y = y(:); ins = ins(:);                          % column vectors
an = double(categorical(string(animal(:))));                       % numeric animal IDs
ok = ~isnan(depth) & ~isnan(y);                                    % complete rows
[g, ~] = findgroups(ins);                                          % insertion groups
nPer = accumarray(g, ok);                                          % complete units per insertion
ok = ok & nPer(g) >= params.minPerIns;                             % insertions large enough to rank
depth = depth(ok); y = y(ok); ins = ins(ok); an = an(ok);          % analysed rows
assert(numel(unique(ins)) >= 2, 'depthCorrTest: fewer than 2 usable insertions.');   % nothing to test otherwise

R.rho = rhoW(depth, y, ins);                                       % observed within-insertion rank correlation
null = zeros(params.nPerm, 1);
[gi, ~] = findgroups(ins);
for b = 1:params.nPerm
    dP = depth;
    for k = 1:max(gi)
        ii = find(gi == k); dP(ii) = depth(ii(randperm(numel(ii))));   % permute depths within the insertion
    end
    null(b) = rhoW(dP, y, ins);                                    % statistic under H0
end
R.p = min(1, (sum(abs(null) >= abs(R.rho)) + 1) / (params.nPerm + 1));   % two-tailed permutation p
bt = hierBootStatistic(@(idx) rhoW(depth(idx), y(idx), ins(idx)), params.nBoot, an, ins);   % hierarchical bootstrap
R.ciLo = prctile(bt, 2.5); R.ciHi = prctile(bt, 97.5);             % 95% CI
R.n = numel(y); R.nIns = numel(unique(ins)); R.nAnimals = numel(unique(an));   % sample sizes
end

% =========================================================================
function r = rhoW(d, y, ins)
% Correlation of within-insertion centred ranks (ties averaged).
rd = nan(size(d)); ry = nan(size(y));                              % centred ranks
for g = unique(ins)'
    ii = ins == g;
    if sum(ii) < 2, continue, end                                  % a single unit carries no within-insertion information
    a = tiedrank(d(ii)); b = tiedrank(y(ii));                      % ranks inside this insertion
    rd(ii) = a - mean(a); ry(ii) = b - mean(b);                    % centred
end
k = ~isnan(rd) & ~isnan(ry);
den = sqrt(sum(rd(k).^2) * sum(ry(k).^2));                         % normaliser
if den == 0, r = NaN; else, r = sum(rd(k) .* ry(k)) / den; end     % pooled correlation of centred ranks
end
