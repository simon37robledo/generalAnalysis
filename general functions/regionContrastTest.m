function R = regionContrastTest(y, isTarget, rec, animal, params)
% regionContrastTest  Within-recording contrast of one region against the other labelled regions.
%
%   R = regionContrastTest(y, isTarget, rec, animal)
%
%   Statistic: in every recording that has >= minPerSide units inside the target region AND
%   >= minPerSide units outside it (other labelled regions), d = mean(y, target) - mean(y, other);
%   the statistic is the average of d over recordings weighted by 1/(1/nIn + 1/nOut). Only
%   within-recording differences count, so recordings that differ in depth coverage, state or
%   response level cannot create a region effect (same logic as depthCorrTest).
%   p   : region labels shuffled among the units of each recording (two-tailed, nPerm shuffles).
%         This is a within-recording test: it does not generalise across animals. pBoot (two-tailed,
%         from the same bootstrap as the CI) does; with k animals it cannot go below ~2/k^k.
%   CI  : hierarchical bootstrap (animal -> recording -> unit) of the same statistic
%         (hierBootStatistic); resamples with no usable recording are counted, not hidden.
%   y may be logical (a fraction is then the contrast of proportions).
%
%   INPUTS  y, isTarget (logical), rec (recording ID, unique across animals), animal: one entry
%           per unit; units without a label ("unlabelled") must be removed BEFORE the call.
%   OUTPUT  struct: diff, ciLo, ciHi, p, nIn, nOut, nRec (usable recordings), nAnimals (usable
%           animals), nBootNaN. diff is NaN (with a reason in R.note) when fewer than
%           minRec recordings or minAnimals animals are usable.

arguments
    y (:,1) double                                                % measure per unit
    isTarget (:,1) logical                                        % unit lies in the target region
    rec (:,1) double                                              % recording per unit
    animal                                                        % animal per unit (string/categorical/numeric)
    params.nPerm (1,1) double = 10000                             % within-recording shuffles
    params.nBoot (1,1) double = 5000                              % hierarchical bootstrap resamples (0 = no CI)
    params.minPerSide (1,1) double = 3                            % units needed inside AND outside the target in a recording
    params.minRec (1,1) double = 3                                % usable recordings needed
    params.minAnimals (1,1) double = 2                            % usable animals needed
end

assert(numel(y) == numel(isTarget) && numel(y) == numel(rec) && numel(y) == numel(animal), ...
    'regionContrastTest: y, isTarget, rec and animal must have one entry per unit.');   % alignment invariant
assert(~any(isnan(y)), 'regionContrastTest: y has NaN values (drop them before the call).');   % never skip silently
an = double(categorical(string(animal(:))));                      % numeric animal IDs
R = struct('diff', NaN, 'ciLo', NaN, 'ciHi', NaN, 'p', NaN, 'pBoot', NaN, 'nIn', sum(isTarget), 'nOut', sum(~isTarget), ...
    'nRec', 0, 'nAnimals', 0, 'nBootNaN', NaN, 'note', "");       % output (NaN until the checks pass)

[g, recList] = findgroups(rec);                                   % recording groups
nG = numel(recList);                                              % number of recordings
nInG = accumarray(g, isTarget);                                   % target units per recording
nOutG = accumarray(g, ~isTarget);                                 % other units per recording
usable = nInG >= params.minPerSide & nOutG >= params.minPerSide;  % recordings that can contribute (counts do not change under shuffling)
animalOfG = accumarray(g, an, [], @(v) v(1));                     % animal of each recording
assert(all(accumarray(g, an, [], @(v) numel(unique(v))) == 1), ...
    'regionContrastTest: a recording belongs to more than one animal (rec IDs must be unique across animals).');   % key check
R.nRec = sum(usable); R.nAnimals = numel(unique(animalOfG(usable)));   % usable recordings and animals
if R.nRec < params.minRec || R.nAnimals < params.minAnimals       % too little data: report why, no test
    R.note = sprintf('%d usable recordings / %d animals', R.nRec, R.nAnimals);
    return
end

rows = arrayfun(@(k) find(g == k), (1:nG)', 'UniformOutput', false);   % unit rows of every recording
rows = rows(usable);                                              % keep usable recordings only
R.diff = contrastStat(y, isTarget, rows);                         % observed statistic

nullD = nan(params.nPerm, 1);                                     % shuffled statistics
for k = 1:params.nPerm
    lab = isTarget;                                               % labels to shuffle
    for m = 1:numel(rows)
        r = rows{m};                                              % units of this recording
        lab(r) = isTarget(r(randperm(numel(r))));                 % shuffle the region label among them
    end
    nullD(k) = contrastStat(y, lab, rows);                        % statistic under no region effect
end
R.p = (1 + sum(abs(nullD) >= abs(R.diff) - 1e-12)) / (params.nPerm + 1);   % two-tailed permutation p

if params.nBoot > 0                                               % hierarchical bootstrap CI
    sf = @(idx) bootContrast(y(idx), isTarget(idx), rec(idx), params.minPerSide);   % statistic on drawn rows
    bt = hierBootStatistic(sf, params.nBoot, an, rec);            % animal -> recording -> unit resampling
    R.nBootNaN = sum(isnan(bt));                                  % resamples without a usable recording
    if R.nBootNaN > 0.05 * params.nBoot                           % explicit warning: CI may be unreliable
        warning('regionContrastTest:bootNaN', '%.1f%% of bootstrap resamples had no usable recording.', 100 * R.nBootNaN / params.nBoot);
    end
    btOk = bt(~isnan(bt));                                        % defined resamples
    ci = prctile(btOk, [2.5 97.5]);                               % 95% percentile CI over defined resamples
    R.ciLo = ci(1); R.ciHi = ci(2);
    nB = numel(btOk);                                             % defined resamples
    R.pBoot = min(1, 2 * min((1 + sum(btOk <= 0)) / (nB + 1), (1 + sum(btOk >= 0)) / (nB + 1)));   % two-tailed bootstrap p (across-animal inference; floor ~2/k^k with k animals)
end
end

function s = contrastStat(y, isT, rows)
% Weighted mean over recordings of mean(y | target) - mean(y | other); rows = unit rows per recording.
num = 0; den = 0;                                                 % weighted sum and weights
for m = 1:numel(rows)
    r = rows{m};                                                  % units of this recording
    a = y(r(isT(r))); b = y(r(~isT(r)));                          % inside / outside the target
    w = 1 / (1 / numel(a) + 1 / numel(b));                        % precision-style weight
    num = num + w * (mean(a) - mean(b)); den = den + w;           % accumulate
end
s = num / den;                                                    % pooled within-recording difference
end

function s = bootContrast(y, isT, rec, minPerSide)
% Statistic on a bootstrap sample: groups the drawn rows by recording; NaN if no recording is usable.
[g, ~] = findgroups(rec);                                         % recordings in the sample
num = 0; den = 0;                                                 % weighted sum and weights
for k = 1:max(g)
    a = y(g == k & isT); b = y(g == k & ~isT);                    % inside / outside the target
    if numel(a) < minPerSide || numel(b) < minPerSide, continue, end   % recording unusable in this resample
    w = 1 / (1 / numel(a) + 1 / numel(b));                        % weight
    num = num + w * (mean(a) - mean(b)); den = den + w;          % accumulate
end
if den == 0, s = NaN; else, s = num / den; end                    % NaN when nothing usable
end
