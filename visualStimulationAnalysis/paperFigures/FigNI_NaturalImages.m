%% NATURAL IMAGES — intact vs scrambled (5 image pairs, PV51 and PV233)
% Stimulus: 10 images (5 intact natural scenes + a texture-scrambled version of each, same
% coarse layout and colours), 2 s on, 1 s ITI, 15 repetitions each, random order.
% Responsiveness: the figure design (figureDesign.m) — 300 ms from onset vs 300 ms pre-trial,
% max over the 10 images after the sign flip (null FPR on this layout: FigS_NullCalibration).
% Intact vs scrambled, in NI-responsive units present during the NI block:
%   per unit, per trial d = (spikes in [onset, onset+300) - spikes in [onset-300, onset)) / 0.3 s;
%   index = mean d over the 5 intact images - mean d over the 5 scrambled images (spikes/s),
%   and the same divided by the across-trial SD of the per-trial baseline rate (z units;
%   SD floored at the 25th percentile over units, so silent baselines do not explode it).
%   Per-unit p: intact/scrambled labels permuted within each pair across its 30 trials
%   (two-tailed). Selection uses all 10 images, so it favours neither label; the null
%   calibration cell below checks the selection + test pipeline gives ~5% FPR.
%   Cross-animal: hierarchical bootstrap (animal -> insertion -> unit, hierBootMatchFreq),
%   two-tailed; with 2 animals the claim rests on both animals agreeing (per-animal rows).

d = figureDesign();                                               % shared design
paperFigPrefix('figNI');                                          % saved panels are named figNI_<panel>
expsNI = [89 90 92 93 95 96 97];                                  % recordings with natural images (PV51: 89-90, PV233: 92-97)
assert(all(ismember(expsNI, d.exps)), 'FigNI: NI recordings must be in figureDesign exps.');
win    = 300;                                                     % response / baseline window (ms), as the design
nPerm  = 10000;                                                   % within-pair label permutations per unit
nBoot  = 10000;                                                   % hierarchical bootstrap replicates
cacheF = fullfile(d.paperFigs, 'FigNI_units.mat');                % per-unit table (cache)
outCsv = fullfile(d.paperFigs, 'Stats summary tables', 'NI_intact_vs_scrambled.csv');   % supplementary table
sa     = d.statsArgs;                                             % figure-design statistics arguments
arg    = @(name) sa{find(strcmp(sa(1:2:end), name)) * 2};         % value of one design argument

%% Per-recording statistics and per-unit intact - scrambled index (cached)
U = table();                                                      % one row per good unit
for ex = expsNI
    NP = loadNPclassFromTable(ex);                                % recording
    animal = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once'));   % animal ID
    o  = imageAnalysis(NP);                                       % natural-image object
    assert(~isempty(o.VST), 'FigNI: no NI stimulus in recording %d.', ex);
    R  = o.ResponseWindow('overwrite', arg('overwriteResponse'), 'durationWindow', 100);   % onsets and image index (as runStimStats)
    S  = o.StatisticsPerNeuron('overwrite', arg('overwriteStats'), 'BaseRespWindow', win, ...
        'SpatialGridMode', arg('SpatialGridMode'), 'maxCategory', arg('maxCategory'), 'applyFDR', false, ...
        'durationWindow', [], 'OffWindow', false, 'useSegments', false, 'nSegments', arg('nSegments'), ...
        'SDGWindows', arg('sdgWindows'), 'CategoryLevels', "all", 'Lock', {});   % the design's single-window test (runStimStats arguments)
    % Images: intact names end in 'Norm', scrambled in 'NormShuf(f)'; pair = name before 'Norm'
    [~, nm] = cellfun(@fileparts, o.VST.imgNames, 'UniformOutput', false);   % image names in index order
    nm = string(nm);
    isScr = contains(nm, "Shuf");                                 % scrambled images
    pairName = extractBefore(nm, "Norm");                         % pair key, e.g. 'NormalLizard'
    [pairs, ~, pairIdx] = unique(pairName, 'stable');             % pair of every image
    assert(numel(nm) == 10 && numel(pairs) == 5 && all(accumarray(pairIdx(:), isScr(:)) == 1), ...
        'FigNI: recording %d does not hold 5 intact/scrambled pairs.', ex);   % one scrambled per pair
    % Trials: onsets (ms) and image index
    on  = R.C(:, 1); img = R.C(:, 2);                             % one row per trial
    assert(numel(on) == 150 && all(accumarray(img, 1) == 15), 'FigNI: recording %d: expected 15 trials per image.', ex);
    % Spike counts of every good unit (same unit order as the statistics)
    ps  = NP.convertPhySorting2tIc(o.spikeSortingFolder);         % phy sorting
    isG = string(ps.label') == "good";                            % somatic units
    gU  = ps.ic(:, isG);                                          % their spike-sorting columns
    nU  = size(gU, 2);
    assert(nU == numel(S.pvalsResponse), 'FigNI: recording %d: %d good units, %d in the statistics.', ex, nU, numel(S.pvalsResponse));
    Mb = sum(BuildBurstMatrix(gU, round(ps.t), round(on' - win), win), 3);   % [trials x units] baseline spikes
    Mr = sum(BuildBurstMatrix(gU, round(ps.t), round(on'), win), 3);         % [trials x units] response spikes
    assert(isequal(size(Mb), [150 nU]), 'FigNI: unexpected count matrix size.');
    D  = (Mr - Mb) / (win / 1000);                                % per-trial evoked rate (spikes/s)
    sdB = std(Mb / (win / 1000), 0, 1)';                          % across-trial SD of the baseline rate (spikes/s)
    % Observed index and within-pair permutation null (same permutations for every unit)
    tI = reshape(~isScr(img), [], 1); tP = reshape(pairIdx(img), [], 1);   % per trial (columns): intact? pair?
    obsP = pairDiff(D, tI, tP);                                   % [5 x units] intact - scrambled per pair
    obs  = mean(obsP, 1)';                                        % index: mean over pairs (spikes/s)
    rng(ex);                                                      % reproducible per recording
    W = permWeights(tI, tP, nPerm);                               % [trials x nPerm] weights of permuted labelings
    nullT = (W' * D)';                                            % [units x nPerm] permuted indices (linear in the trials)
    assert(max(abs(labelWeights(tI) .' * D - obs')) < 1e-9, 'FigNI: weight form disagrees with pairDiff.');   % same statistic
    pIS = (sum(abs(nullT) >= abs(obs) - 1e-12, 2) + 1) / (nPerm + 1);   % two-tailed permutation p
    phyID = ps.phy_ID(isG); phyID = phyID(:);                     % phy IDs of the good units (column)
    T = table(repmat(ex, nU, 1), repmat(animal, nU, 1), phyID, S.pvalsResponse(:) < 0.05, ...
        statPresence(S, "NI", ""), obs, sdB, pIS, obsP', ...
        'VariableNames', {'exp','animal','phyID','responsive','present','dRate','sdBase','pIS','pairDiff'});   % per unit
    U = [U; T];                                                   %#ok<AGROW>
end
pairNames = pairs;                                                % pair order (identical in every recording, asserted below)
save(cacheF, 'U', 'pairNames');

%% Cross-animal statistics
load(cacheF, 'U', 'pairNames');
sdFloor = prctile(U.sdBase(U.present), 25);                       % SD floor (25th percentile over present units)
U.z = U.dRate ./ max(U.sdBase, sdFloor);                          % index in baseline-SD units
A = U(U.responsive & U.present, :);                               % analysable units
fprintf('NI: %d good units, %d present, %d responsive and present (%d animals, %d insertions)\n', ...
    height(U), sum(U.present), height(A), numel(unique(A.animal)), numel(unique(A.exp)));
rng(0);
Res = table();                                                    % one row per metric x animal set
for m = ["z", "dRate"]
    for an = ["all", reshape(unique(A.animal), 1, [])]            % row: a for-loop iterates over columns
        a = A; if an ~= "all", a = A(A.animal == an, :); end      % all animals or one animal
        b = hierBootMatchFreq(a.(m), nBoot, double(categorical(a.animal)), a.exp);   % (data, nrep, ANIMAL, INSERTION)
        p = min(1, 2 * min(mean(b <= 0), mean(b >= 0)));          % two-tailed bootstrap p
        Res = [Res; table(m, an, height(a), numel(unique(a.exp)), mean(b), prctile(b, 2.5), prctile(b, 97.5), p, ...
            sum(a.pIS < 0.05 & a.(m) > 0), sum(a.pIS < 0.05 & a.(m) < 0), ...
            'VariableNames', {'metric','animal','nUnits','nInsertions','mean','ciLo','ciHi','p','nPreferIntact','nPreferScrambled'})];   %#ok<AGROW>
    end
end
disp(Res);
writetable(Res, outCsv);                                          % supplementary CSV

%% Panel: intact - scrambled per unit (animal-colored), bootstrap mean and 95% CI
anAll = unique(A.animal); cmap = lines(numel(anAll));             % one color per animal
[~, ai] = ismember(A.animal, anAll);
fig = figure('Position', [100 100 360 380]); hold on;
sel = A.pIS < 0.05;                                               % units with a significant preference
swarmchart(ones(sum(~sel), 1), A.z(~sel), 14, cmap(ai(~sel), :), 'filled', 'MarkerFaceAlpha', 0.5, 'XJitterWidth', 0.5);
swarmchart(ones(sum(sel), 1), A.z(sel), 18, cmap(ai(sel), :), 'filled', 'MarkerEdgeColor', 'k', 'XJitterWidth', 0.5);
r = Res(Res.metric == "z" & Res.animal == "all", :);
errorbar(1, r.mean, r.mean - r.ciLo, r.ciHi - r.mean, 'k', 'LineWidth', 1.5, 'CapSize', 8);
yline(0, 'k:');
set(gca, 'XTick', 1, 'XTickLabel', "NI-responsive units"); ylabel('Intact - scrambled (z, baseline SD units)');
title(sprintf('p = %.3g, n = %d units, %d animals', r.p, r.nUnits, numel(anAll)));
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'IntactVsScrambled_swarm.jpg']), 'Resolution', 300);

%% Panel: per image pair (mean over units, hierarchical bootstrap 95% CI)
Pz = A.pairDiff ./ max(A.sdBase, sdFloor);                        % [units x 5] per-pair index, z units
fig = figure('Position', [100 100 520 360]); hold on;
for q = 1:numel(pairNames)
    b = hierBootMatchFreq(Pz(:, q), nBoot, double(categorical(A.animal)), A.exp);
    bar(q, mean(b), 0.6, 'FaceColor', [0.75 0.75 0.75]);
    errorbar(q, mean(b), mean(b) - prctile(b, 2.5), prctile(b, 97.5) - mean(b), 'k.', 'LineWidth', 1.2);
end
yline(0, 'k:');
set(gca, 'XTick', 1:numel(pairNames), 'XTickLabel', pairNames, 'XTickLabelRotation', 30);
ylabel('Intact - scrambled (z)'); title('Per image pair (descriptive)');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'IntactVsScrambled_perPair.jpg']), 'Resolution', 300);

%% Null calibration: selection (max over 10 images, sign flip) + intact/scrambled test
% Simulated units with the SAME evoked response to every image (H0 for intact vs scrambled),
% on recording 92's trial layout, Bernoulli spiking in 1 ms bins, gamma trial gain
% (overdispersion). Selection re-implements the design's max-statistic sign-flip test
% (max over images of the mean per-trial evoked rate). Target: ~5% of selected units with
% pIS < 0.05.
rng(7);
NP = loadNPclassFromTable(92); o = imageAnalysis(NP); R = o.ResponseWindow;   % real layout
img = R.C(:, 2); [~, nm] = cellfun(@fileparts, o.VST.imgNames, 'UniformOutput', false); nm = string(nm);
tI = reshape(~contains(nm(img), "Shuf"), [], 1); [~, ~, pIdx] = unique(extractBefore(nm, "Norm"), 'stable');
tP = reshape(pIdx(img), [], 1);                                   % per trial (columns): intact? pair?
nSim = 4000; nPermSim = 2000; nFlip = 2000;                       % simulated units, permutations, sign flips
amp  = [0 2 5 10 20];                                             % evoked rate added to every image (spikes/s)
Wsim = permWeights(tI, tP, nPermSim);                             % permuted labelings (shared by the simulated units)
w0   = labelWeights(tI);                                          % observed labeling
sel = false(nSim, 1); pSim = nan(nSim, 1);
for s = 1:nSim
    lam0 = exp(log(3) + randn);                                   % baseline rate (spikes/s), lognormal around 3
    lam1 = lam0 + amp(randi(numel(amp))) * gamrnd(5, 1/5, 150, 1);   % evoked rate per trial, same for all images
    b = binornd(win, min(lam0 / 1000, 1), 150, 1);                % baseline spikes
    rr = binornd(win, min(lam1 / 1000, 1));                       % response spikes
    Dv = (rr - b) / (win / 1000);                                 % per-trial evoked rate
    obsMax = max(accumarray(img, Dv, [], @mean));                 % observed: max over images
    fl = (2 * (rand(150, nFlip) > 0.5) - 1) .* Dv;                % sign-flipped trials
    nullMax = max(splitapply(@(x) mean(x, 1), fl, img), [], 1);   % null max over images
    sel(s) = mean(nullMax >= obsMax) < 0.05;                      % selected (responsive)
    o1 = w0' * Dv;                                                % observed intact - scrambled index
    pSim(s) = (sum(abs(Wsim' * Dv) >= abs(o1) - 1e-12) + 1) / (nPermSim + 1);   % two-tailed permutation p
end
fprintf('Null calibration NI intact-scrambled: FPR all units %.3f (n=%d); selected units %.3f (n=%d)\n', ...
    mean(pSim < 0.05), nSim, mean(pSim(sel) < 0.05), sum(sel));

%% Local functions
function P = pairDiff(D, tI, tP)
% Intact - scrambled mean of D per pair: D [trials x units], tI logical intact per trial,
% tP pair index per trial. Returns [nPairs x units].
nP = max(tP); P = zeros(nP, size(D, 2));
for q = 1:nP
    P(q, :) = mean(D(tP == q & tI, :), 1) - mean(D(tP == q & ~tI, :), 1);   % this pair's difference
end
end

function w = labelWeights(tI)
% Trial weights giving mean over pairs of (intact mean - scrambled mean) as w' * D, valid when
% every pair has the same number of intact and scrambled trials (asserted: 15 + 15 per pair).
nPairs = 5; nPer = sum(tI) / nPairs;                              % trials per image
assert(nPer == 15 && sum(~tI) == nPairs * nPer, 'FigNI: unbalanced intact/scrambled trials.');
w = (2 * double(tI(:)) - 1) / (nPer * nPairs);                    % +1/(15*5) intact, -1/(15*5) scrambled
end

function W = permWeights(tI, tP, nPerm)
% Weights of nPerm labelings with intact/scrambled labels shuffled within each pair's trials.
W = zeros(numel(tI), nPerm);                                      % [trials x permutations]
for k = 1:nPerm
    tIp = tI;                                                     % permuted labels
    for q = 1:max(tP)
        r = find(tP == q);                                        % this pair's 30 trials
        tIp(r) = tI(r(randperm(numel(r))));                       % shuffle labels within the pair
    end
    W(:, k) = labelWeights(tIp);                                  % weights of this labeling
end
end
