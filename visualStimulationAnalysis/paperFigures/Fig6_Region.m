%% FIGURE 6 (region panels) — Do responses differ between brain regions? (atlas animals: PV97, PV132, PV140, PV233)
% Uses the per-unit cache written by Fig6_Depth.m (Fig6_atlas_units.mat: measures + atlas region per unit).
% Regions come from the brain-clearing atlas (getUnitAtlasPositions). About 40% of units sit in
% "Unlabeled atlas" (ID 0): shown as a grey category, never part of a region test.
% Test (regionContrastTest): within-recording contrast of one region against the other LABELLED
% regions (recordings with >= 3 units on each side), p from region labels shuffled within each
% recording, 95% CI from a hierarchical bootstrap (animal -> recording -> unit). Tested regions
% need >= 3 usable recordings in >= 2 animals (aDVR, striatum, claustrum); the other regions are
% shown descriptively. Benjamini-Hochberg across all region tests (separate from the depth family).
% Only 4 animals: no region claim is made from one animal; PV233 has no aDVR units and PV97 few
% (see Fig6_region_counts.csv). The animal subset follows from atlas availability, not from responses.

d = figureDesign();                                               % shared design
paperFigPrefix('fig6');                                           % saved panels are named fig6_<panel>
plotLab = @(s) replace(replace(replace(s, "RG", "SB"), "SDGm", "MG"), "SDGs", "SG");   % names shown on the panels
cacheF = fullfile(d.paperFigs, 'Fig6_atlas_units.mat');           % per-unit table written by Fig6_Depth.m
L = load(cacheF, 'M', 'D');                                       % unit x stimulus measures and unit table
longNames  = ["anterior dorsal ventricular ridge", "striatum", "claustrum", "basal forebrain", "anterior dorsal cortex", "posterior dorsal cortex"];   % atlas names
shortNames = ["aDVR", "Striatum", "Claustrum", "Basal forebr.", "Ant. dors. ctx", "Post. dors. ctx"];   % panel names
regAll = [shortNames, "Unlabelled"];                              % x-axis order (grey category last)
testRegs = ["aDVR", "Striatum", "Claustrum"];                     % regions with enough animals for a test
M = labelRegions(L.M, longNames, shortNames);                     % add M.reg (short name or "Unlabelled")
D = labelRegions(L.D, longNames, shortNames);                     % same for the unit table
anAll = unique(M.animal);                                         % animals, fixed order for every panel
cmap = animalColours(anAll);                                      % fixed colour per animal
stims = ["MB", "RG", "MBR"];                                      % stimuli with enough units per region (gratings: too few responsive units)
rng(1);                                                          % reproducible permutations and bootstrap
tests = table();                                                  % one row per region test

%% Unit counts per animal x region x stimulus (what the tests can rely on)
cnt = table();                                                    % long table of counts
for s = stims
    a = M(M.stim == s & M.present, :);                            % units present during this stimulus
    G = groupsummary(a, {'animal', 'reg'}, 'sum', 'responsive');  % units and responsive units per animal x region
    G.stim = repmat(s, height(G), 1);                             % stimulus label
    cnt = [cnt; G(:, {'stim', 'animal', 'reg', 'GroupCount', 'sum_responsive'})];   %#ok<AGROW>
end
cnt.Properties.VariableNames = {'stim', 'animal', 'region', 'nUnits', 'nResponsive'};   % clear names
writetable(cnt, fullfile(d.paperFigs, 'Fig6_region_counts.csv'));  % counts table (supplement)

%% Measures: fraction responsive, z and delta rate of responsive units (MB, SB, MBR)
meas = struct('prefix', {}, 'T', {}, 'y', {}, 'binary', {});      % list of measures to test and plot
for s = stims
    a = M(M.stim == s & M.present, :);                            % present units
    meas(end+1) = struct('prefix', s + ": responsive", 'T', a, 'y', double(a.responsive), 'binary', true);   %#ok<SAGROW>
    r = a(a.responsive, :);                                       % responsive units
    meas(end+1) = struct('prefix', s + ": z", 'T', r, 'y', r.z, 'binary', false);          %#ok<SAGROW>
    meas(end+1) = struct('prefix', s + ": delta rate", 'T', r, 'y', r.dRate, 'binary', false);   %#ok<SAGROW>
end
% MB - SB z: units present in both blocks and responsive to either (as Figure 1)
mb = M(M.stim == "MB", {'exp','phyID','animal','reg','z','responsive','present'});
rg = M(M.stim == "RG", {'exp','phyID','z','responsive','present'});
rg.Properties.VariableNames(3:5) = {'zRG', 'respRG', 'presRG'};   % explicit names (no join suffixes)
J = innerjoin(mb, rg, 'Keys', {'exp','phyID'});                   % pair by key, never by row
J = J(J.present & J.presRG & (J.responsive | J.respRG), :);       % OR-set, present in both
meas(end+1) = struct('prefix', "MB - SB z", 'T', J, 'y', J.z - J.zRG, 'binary', false);
% Spatial tuning index / stripe RF (Figure 2 table) and direction selectivity (Figure 5 table)
L2 = load('W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_40-97_SpatialTuningIndex_linearlyMovingBall-rectGrid_RF_prefDir_allResp_splitHalf_presence.mat', 'tbl');
t2 = L2.tbl(L2.tbl.onOff == 1 & L2.tbl.sizeIdx == 1 & L2.tbl.lumValue == 255, :);   % white only
t2.exp = str2double(string(t2.experimentNum));                    % numeric recording ID (categorical -> value)
for st = ["linearlyMovingBall", "rectGrid"]
    lab = "MB"; if st == "rectGrid", lab = "SB"; end              % paper label
    a = innerjoin(t2(t2.stimulus == st, {'exp','phyID','L_amplitude_diff','isStripe'}), ...
        D(:, {'exp','phyID','animal','reg'}), 'Keys', {'exp','phyID'});   % keeps only the atlas animals (the Fig 2 table has all animals)
    a = a(~isnan(a.L_amplitude_diff), :);                         % tuning index defined
    meas(end+1) = struct('prefix', lab + ": spatial tuning index", 'T', a, 'y', a.L_amplitude_diff, 'binary', false);   %#ok<SAGROW>
    meas(end+1) = struct('prefix', lab + ": stripe RF", 'T', a, 'y', double(a.isStripe), 'binary', true);               %#ok<SAGROW>
end
L5 = load(fullfile(d.paperFigs, 'Fig5_directionTuning_units.mat'), 'U');
u5 = L5.U(L5.U.stim == "MB" & L5.U.responsive & L5.U.present & ~isnan(L5.U.DSIcorr), {'exp','phyID','DSIcorr','pDSI'});
u5 = innerjoin(u5, D(:, {'exp','phyID','animal','reg'}), 'Keys', {'exp','phyID'});   % keeps only the atlas animals
meas(end+1) = struct('prefix', "MB: DSI", 'T', u5, 'y', u5.DSIcorr, 'binary', false);
meas(end+1) = struct('prefix', "MB: direction-selective", 'T', u5, 'y', double(u5.pDSI < 0.05), 'binary', true);

%% Region tests (labelled units only)
for k = 1:numel(meas)
    T = meas(k).T; y = meas(k).y;                                 % this measure
    assert(numel(y) == height(T), 'Fig6_Region: %s: y and table differ in length.', meas(k).prefix);   % alignment invariant
    keep = T.reg ~= "Unlabelled";                                 % labelled units only
    for r = testRegs
        R = regionContrastTest(y(keep), T.reg(keep) == r, T.exp(keep), T.animal(keep));   % within-recording contrast
        row = struct('name', meas(k).prefix + " | " + r, 'measure', meas(k).prefix, 'region', r, 'diff', R.diff, ...
            'ciLo', R.ciLo, 'ciHi', R.ciHi, 'p', R.p, 'pBoot', R.pBoot, 'nIn', R.nIn, 'nOut', R.nOut, 'nRec', R.nRec, ...
            'nAnimals', R.nAnimals, 'note', string(R.note));      % one table row
        tests = [tests; struct2table(row)];                       %#ok<AGROW>
    end
end
ok = ~isnan(tests.p);                                             % tests that could be run
q = nan(height(tests), 1); q(ok) = bhAdjust(tests.p(ok));         % Benjamini-Hochberg over the region family
tests.qBH = q;
disp(tests(:, {'name','diff','ciLo','ciHi','p','pBoot','qBH','nIn','nOut','nRec','nAnimals','note'}));
writetable(tests, fullfile(d.paperFigs, 'Fig6_region_tests.csv'));

%% Leave-one-animal-out for tests surviving BH (q < 0.05): drop each animal, rerun the permutation test
loao = table();                                                   % one row per surviving test x dropped animal
for i = find(tests.qBH < 0.05)'
    k = find([meas.prefix] == tests.measure(i));                  % measure of this test
    T = meas(k).T; y = meas(k).y; keep = T.reg ~= "Unlabelled";   % labelled units
    for a = anAll'
        kk = keep & T.animal ~= a;                                % drop this animal
        R = regionContrastTest(y(kk), T.reg(kk) == tests.region(i), T.exp(kk), T.animal(kk), nBoot=0);   % permutation test without the animal
        loao = [loao; table(tests.name(i), a, R.diff, R.p, R.nAnimals, string(R.note), ...
            'VariableNames', {'name','dropped','diff','p','nAnimals','note'})];   %#ok<AGROW>
    end
end
disp(loao);
writetable(loao, fullfile(d.paperFigs, 'Fig6_region_LOAO.csv'));

%% Depth within aDVR only (region and depth are collinear: does depth matter inside one region?)
dtests = table();                                                 % depth-within-region tests
for s = ["MB", "RG"]
    a = M(M.stim == s & M.present & M.reg == "aDVR" & M.depthValid, :);   % aDVR units with a depth
    r = a(a.responsive, :);                                       % responsive aDVR units
    dtests = runDepthTest(dtests, s + ": responsive vs depth (aDVR only)", a, double(a.responsive));
    dtests = runDepthTest(dtests, s + ": z vs depth (aDVR only)", r, r.z);
end
disp(dtests);
writetable(dtests, fullfile(d.paperFigs, 'Fig6_region_depth_within_aDVR.csv'));

%% Calibration: shuffling the region label within each recording must reject ~5% of the time
nCal = 200; rej = 0;                                              % repetitions, rejections
a = meas(1).T; y = meas(1).y; keep = a.reg ~= "Unlabelled";      % MB responsive measure
a = a(keep, :); y = y(keep);                                      % labelled units
for k = 1:nCal
    lab = a.reg == "aDVR";                                        % true target labels
    [g, ~] = findgroups(a.exp);                                   % recordings
    for m = 1:max(g), idx = find(g == m); lab(idx) = lab(idx(randperm(numel(idx)))); end   % shuffle within recording
    Rk = regionContrastTest(y, lab, a.exp, a.animal, nPerm=300, nBoot=0);   % test on shuffled labels
    rej = rej + (Rk.p < 0.05);
end
fprintf('Fig6_Region calibration: %d/%d shuffled datasets rejected at 0.05 (expect ~%d).\n', rej, nCal, round(0.05 * nCal));

%% Panels: fraction responsive, z, delta rate by region (MB, SB, MBR)
for what = ["responsive", "z", "delta rate"]
    fig = figure('Position', [100 100 1500 380]); tl = tiledlayout(1, numel(stims), 'TileSpacing', 'compact');
    for s = stims
        nexttile(tl);
        k = find([meas.prefix] == s + ": " + what);               % measure index
        ylab = what; if what == "responsive", ylab = 'Fraction responsive'; end
        if what == "z", ylab = 'z-score (responsive units)'; end
        if what == "delta rate", ylab = '\Delta rate (spk/s, responsive units)'; end
        plotRegion(meas(k).T, meas(k).y, meas(k).prefix, tests, anAll, cmap, ylab, plotLab(s), regAll);
    end
    animalLegend(nexttile(tl, numel(stims)), anAll, cmap);        % legend in the last tile
    title(tl, 'Dots: animal means (>= 5 units); line: pooled mean; p: region vs other labelled regions, within recordings');
    exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'region_' char(replace(what, " ", "")) '.jpg']), 'Resolution', 300);
end

%% Panel: MB - SB z by region (units present in both, responsive to either)
fig = figure('Position', [100 100 520 380]);
k = find([meas.prefix] == "MB - SB z");                           % measure index
plotRegion(meas(k).T, meas(k).y, meas(k).prefix, tests, anAll, cmap, 'z_{MB} - z_{SB}', 'MB - SB', regAll);
yline(0, 'k:', 'HandleVisibility', 'off');                        % no preference
animalLegend(gca, anAll, cmap);                                   % legend
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'region_MBminusSB_z.jpg']), 'Resolution', 300);

%% Panels: localized-motion selectivity and direction selectivity by region
names = ["MB: spatial tuning index", "MB: stripe RF", "SB: spatial tuning index", "SB: stripe RF", "MB: DSI", "MB: direction-selective"];
ylabs = {'Spatial tuning index', 'Fraction stripe RFs', 'Spatial tuning index', 'Fraction stripe RFs', 'DSI (corrected)', 'Fraction direction-selective'};
fig = figure('Position', [100 100 1500 700]); tl = tiledlayout(2, 3, 'TileSpacing', 'compact');
for i = 1:numel(names)
    nexttile(tl); k = find([meas.prefix] == names(i));            % measure index
    plotRegion(meas(k).T, meas(k).y, meas(k).prefix, tests, anAll, cmap, ylabs{i}, char(names(i)), regAll);
end
animalLegend(nexttile(tl, 3), anAll, cmap);                       % legend in the top-right tile
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'region_tuning.jpg']), 'Resolution', 300);

%% Panel: depth of the units in each region (why depth and region cannot be separated freely)
fig = figure('Position', [100 100 640 420]); hold on;
Dv = D(D.depthValid, :);                                          % units with an atlas depth
for k = 1:numel(regAll)
    idx = Dv.reg == regAll(k); [~, ai] = ismember(Dv.animal(idx), anAll);   % units of this region, animal index
    jit = (rand(sum(idx), 1) - 0.5) * 0.6;                        % horizontal jitter
    scatter(k + jit, Dv.depthUm(idx), 8, cmap(ai, :), 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
    if any(idx), plot(k + [-0.35 0.35], median(Dv.depthUm(idx)) * [1 1], 'k-', 'LineWidth', 1.5); end   % median depth
end
set(gca, 'YDir', 'reverse', 'XTick', 1:numel(regAll), 'FontSize', 8, 'TickLabelInterpreter', 'tex', 'XTickLabel', ...
    arrayfun(@(r) sprintf('%s\\newline n=%d', r, sum(Dv.reg == r)), regAll, 'UniformOutput', false));   % region name over unit count (TeX \newline)
ylabel('Atlas depth below surface (\mum)'); xlim([0.4 numel(regAll) + 0.6]);
title('Unit depth by region (line: median)');
animalLegend(gca, anAll, cmap);                                   % legend
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'region_depth.jpg']), 'Resolution', 300);

%% Local functions
function T = labelRegions(T, longNames, shortNames)
% Add T.reg: short region name, or "Unlabelled" for atlas ID 0 / -1. Any other labelled region is an error.
T.reg = repmat("Unlabelled", height(T), 1);                       % default: unlabelled
[tf, loc] = ismember(string(T.region), longNames);                % atlas names known to the figure
T.reg(tf) = shortNames(loc(tf));                                  % short names
assert(all(tf == T.regionKnown), 'Fig6_Region: %d units sit in a labelled region that is not in the name list.', sum(tf ~= T.regionKnown));   % no silent mislabel
end

function plotRegion(T, y, prefix, tests, anAll, cmap, ylab, ttl, regAll)
% One panel: per region the animal means (dots, animals with >= 5 units), the pooled mean (line)
% and the region test result above the region label.
hold on;
nR = numel(regAll); nTxt = nan(nR, 1); topY = -inf;               % unit counts and highest plotted value
for k = 1:nR
    idx = T.reg == regAll(k); nTxt(k) = sum(idx);                 % units of this region
    if ~any(idx), continue, end
    for a = 1:numel(anAll)
        ia = idx & T.animal == anAll(a);                          % this animal's units in the region
        if sum(ia) < 5, continue, end                             % too few for a dot
        scatter(k + (a - (numel(anAll) + 1) / 2) * 0.14, mean(y(ia)), 20 + 2 * sqrt(sum(ia)), cmap(a, :), 'filled', ...
            'MarkerFaceAlpha', 0.85, 'HandleVisibility', 'off');  % animal mean, size ~ units
        topY = max(topY, mean(y(ia)));
    end
    plot(k + [-0.38 0.38], mean(y(idx)) * [1 1], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');   % pooled mean (dominated by the largest animal)
    topY = max(topY, mean(y(idx)));
end
if ~isfinite(topY), topY = 1; end                                 % empty panel
yl = [min(0, min(ylim)), topY * 1.35 + (topY <= 0)];             % room for the test text
ylim(yl);
for k = 1:nR                                                      % test results above the region labels
    row = tests(tests.name == prefix + " | " + regAll(k), :);     % test row (only tested regions have one)
    if isempty(row), continue, end
    if isnan(row.p), s = 'not tested'; else, s = sprintf('%+.2f\np=%.3g\npB=%.3g', row.diff, row.p, row.pBoot); end   % permutation p, bootstrap pB
    text(k, yl(2), s, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 7);
end
set(gca, 'XTick', 1:nR, 'XTickLabel', arrayfun(@(k) sprintf('%s\\newline n=%d', regAll(k), nTxt(k)), 1:nR, 'UniformOutput', false), ...
    'TickLabelInterpreter', 'tex', 'FontSize', 7);                % region name over unit count (TeX \newline = two-line labels)
xlim([0.4 nR + 0.6]); ylabel(ylab); title(ttl, 'Interpreter', 'none');
end

function dtests = runDepthTest(dtests, name, tbl, y)
% Depth test within a subset (depthCorrTest); NaN row if too few units/recordings.
row = struct('name', string(name), 'rho', NaN, 'ciLo', NaN, 'ciHi', NaN, 'p', NaN, 'n', height(tbl), 'nIns', NaN, 'nAnimals', NaN);
try
    R = depthCorrTest(tbl.depthUm, y, tbl.exp, tbl.animal);       % within-recording rank correlation
    row.rho = R.rho; row.ciLo = R.ciLo; row.ciHi = R.ciHi; row.p = R.p; row.n = R.n; row.nIns = R.nIns; row.nAnimals = R.nAnimals;
catch ME
    if ~contains(ME.message, 'fewer than 2 usable insertions'), rethrow(ME); end   % only "too few data" is tolerated
    warning('Fig6Region:tooFew', '%s: too few units per recording to test.', name);   % explicit
end
dtests = [dtests; struct2table(row)];
end

function animalLegend(ax, anAll, cmap)
% Animal colour legend on axes ax (dummy markers).
hold(ax, 'on');
h = gobjects(numel(anAll), 1);                                    % one handle per animal
for k = 1:numel(anAll)
    h(k) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', cmap(k, :), 'MarkerEdgeColor', 'none', 'DisplayName', anAll(k));
end
legend(ax, h, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'off');
end

function q = bhAdjust(p)
% Benjamini-Hochberg adjusted p-values.
[ps, o] = sort(p(:)); m = numel(ps);
qs = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));               % step-up adjustment
q = zeros(m, 1); q(o) = qs;
end
