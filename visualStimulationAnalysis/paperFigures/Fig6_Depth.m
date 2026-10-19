%% FIGURE 6 — Does a unit's depth (Z) relate to its responses?
% Depth: um below the brain surface (getUnitDepths). Source "coorZ" (hand-entered
% coor_Z minus probe position) until the collaborator's LFP surface estimate is wired
% into getUnitDepths/readLfpSurface; then set depthSource = "lfp" and rerun.
% X/Y locations come later (brain clearing).
% Test (depthCorrTest): within-insertion rank correlation between depth and a
% measure; p from depths permuted within insertions; 95% CI from a hierarchical
% bootstrap (animal -> insertion -> unit). Benjamini-Hochberg across all tests;
% leave-one-animal-out for tests that survive.

d = figureDesign();                                               % shared design
paperFigPrefix('fig6');                                           % saved panels are named fig6_<panel>
depthSource = "coorZ";                                            % "lfp" once the LFP surface estimate is available
stims = ["MB", "RG", "SDGm", "SDGs", "MBR"];                      % stimuli in the paper
cacheF = fullfile(d.paperFigs, 'Fig6_depth_units.mat');           % per-unit table (cache)

%% Depths and per-unit measures (slow: cached)
D = getUnitDepths(d.exps, source=depthSource);                    % depth per good unit
M = table();
for ex = d.exps
    NP = loadNPclassFromTable(ex);                                % recording
    for s = stims
        m = unitStimMeasures(NP, s);                              % empty if the stimulus was not recorded
        if isempty(m), continue, end
        m.exp = repmat(ex, height(m), 1);                         % recording / insertion ID
        M = [M; m];                                               %#ok<AGROW>
    end
end
M = innerjoin(M, D(:, {'exp','phyID','animal','depthUm'}), 'Keys', {'exp','phyID'});   % pair by key, never by row
assert(all(~isnan(M.depthUm)), 'Fig6: units without depth.');     % every unit must have a depth
save(cacheF, 'M', 'D', 'depthSource');

%% Tests
load(cacheF, 'M', 'D', 'depthSource');
tests = table();                                                  % one row per test
for s = stims
    a = M(M.stim == s & M.present, :);                            % units present during this stimulus
    tests = runTest(tests, s + ": responsive (yes/no)", a, double(a.responsive));
    r = a(a.responsive, :);                                       % responsive units
    tests = runTest(tests, s + ": z", r, r.z);
    tests = runTest(tests, s + ": delta rate", r, r.dRate);
end
% MB vs SB preference: units present in both blocks, responsive to either (as Figure 1)
mb = M(M.stim == "MB", {'exp','phyID','animal','depthUm','z','responsive','present'});
rg = M(M.stim == "RG", {'exp','phyID','z','responsive','present'});
rg.Properties.VariableNames(3:5) = {'zRG','respRG','presRG'};     % explicit names (no join suffixes)
J = innerjoin(mb, rg, 'Keys', {'exp','phyID'});
J = J(J.present & J.presRG & (J.responsive | J.respRG), :);       % OR-set, present in both
tests = runTest(tests, "MB - SB z (OR-set)", J, J.z - J.zRG);
% Direction selectivity (Figure 5 per-unit table)
L5 = load(fullfile(d.paperFigs, 'Fig5_directionTuning_units.mat'), 'U');
u5 = L5.U(L5.U.stim == "MB" & L5.U.responsive & L5.U.present & ~isnan(L5.U.DSIcorr), {'exp','phyID','DSIcorr','pDSI'});
u5 = innerjoin(u5, D(:, {'exp','phyID','animal','depthUm'}), 'Keys', {'exp','phyID'});
tests = runTest(tests, "MB: DSI above chance", u5, u5.DSIcorr);
tests = runTest(tests, "MB: direction-selective (yes/no)", u5, double(u5.pDSI < 0.05));
% Spatial tuning index and stripe RFs (Figure 2 table)
L2 = load('W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_40-97_SpatialTuningIndex_linearlyMovingBall-rectGrid_RF_prefDir_allResp_splitHalf_presence.mat', 'tbl');
t2 = L2.tbl(L2.tbl.onOff == 1 & L2.tbl.sizeIdx == 1 & L2.tbl.lumIdx == 1, :);
t2.exp = str2double(string(t2.experimentNum));                   % numeric recording ID (categorical -> value)
for st = ["linearlyMovingBall", "rectGrid"]
    a = innerjoin(t2(t2.stimulus == st, {'exp','phyID','L_amplitude_diff','isStripe'}), ...
        D(:, {'exp','phyID','animal','depthUm'}), 'Keys', {'exp','phyID'});
    lab = "MB"; if st == "rectGrid", lab = "SB"; end
    tests = runTest(tests, lab + ": spatial tuning index", a, a.L_amplitude_diff);
    tests = runTest(tests, lab + ": stripe RF (yes/no)", a, double(a.isStripe));
end
% Benjamini-Hochberg across the family of depth tests
ok = ~isnan(tests.p);
q = nan(height(tests), 1); q(ok) = bhAdjust(tests.p(ok));
tests.qBH = q;
disp(tests);
writetable(tests, fullfile(d.paperFigs, 'Fig6_depth_tests.csv'));

%% Leave-one-animal-out for tests surviving BH (q < 0.05)
% (rerun by hand for the surviving rows: drop each animal and call depthCorrTest)

%% Panels: z vs depth per stimulus (responsive units), animal-colored
fig = figure('Position', [100 100 1200 300]); tl = tiledlayout(1, numel(stims), 'TileSpacing', 'compact');
anAll = unique(M.animal); cmap = lines(numel(anAll));
for s = stims
    nexttile(tl); hold on;
    r = M(M.stim == s & M.present & M.responsive, :); [~, ai] = ismember(r.animal, anAll);
    scatter(r.z, r.depthUm, 12, cmap(ai, :), 'filled', 'MarkerFaceAlpha', 0.6);
    set(gca, 'YDir', 'reverse'); xlabel('z'); ylabel('Depth below surface (\mum)');
    row = tests(tests.name == s + ": z", :);
    if ~isempty(row) && ~isnan(row.rho)
        title(sprintf('%s: \\rho_w = %.2f [%.2f, %.2f], p = %.3g', s, row.rho, row.ciLo, row.ciHi, row.p));
    else
        title(s + ": too few units");
    end
end
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'z_vs_depth.jpg']), 'Resolution', 300);

%% Panels: fraction responsive by depth (200 um bins), per stimulus
fig = figure('Position', [100 100 1200 300]); tl = tiledlayout(1, numel(stims), 'TileSpacing', 'compact');
edges = floor(min(M.depthUm)/200)*200 : 200 : ceil(max(M.depthUm)/200)*200;   % common depth bins
for s = stims
    nexttile(tl);
    a = M(M.stim == s & M.present, :); b = discretize(a.depthUm, edges);
    fr = accumarray(b, a.responsive, [numel(edges)-1 1], @mean, NaN); nb = accumarray(b, 1, [numel(edges)-1 1]);
    barh(edges(1:end-1) + 100, fr, 1, 'FaceColor', [0.6 0.6 0.6]); set(gca, 'YDir', 'reverse');
    xlabel('Fraction responsive'); ylabel('Depth below surface (\mum)'); title(s);
    for k = find(nb' > 0), text(min(fr(k) + 0.02, 0.9), edges(k) + 100, sprintf('n=%d', nb(k)), 'FontSize', 7); end
    row = tests(tests.name == s + ": responsive (yes/no)", :);
    if ~isempty(row) && ~isnan(row.rho), subtitle(sprintf('\\rho_w = %.2f, p = %.3g', row.rho, row.p)); end
end
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'responsive_vs_depth.jpg']), 'Resolution', 300);

%% Local functions
function tests = runTest(tests, name, tbl, y)
% One depth test (depthCorrTest) appended as a row; NaN row if too few insertions.
row = struct('name', string(name), 'rho', NaN, 'ciLo', NaN, 'ciHi', NaN, 'p', NaN, 'n', height(tbl), 'nIns', NaN, 'nAnimals', NaN);
try
    R = depthCorrTest(tbl.depthUm, y, tbl.exp, tbl.animal);      % within-insertion rank correlation
    row.rho = R.rho; row.ciLo = R.ciLo; row.ciHi = R.ciHi; row.p = R.p; row.n = R.n; row.nIns = R.nIns; row.nAnimals = R.nAnimals;
catch ME
    if ~contains(ME.message, 'fewer than 2 usable insertions'), rethrow(ME); end   % only "too few data" is tolerated
    warning('Fig6:tooFew', '%s: too few units per insertion to test.', name);      % explicit
end
tests = [tests; struct2table(row)];
end

function q = bhAdjust(p)
% Benjamini-Hochberg adjusted p-values.
[ps, o] = sort(p(:)); m = numel(ps);
qs = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));               % step-up adjustment
q = zeros(m, 1); q(o) = qs;
end
