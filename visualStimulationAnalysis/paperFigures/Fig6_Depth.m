%% FIGURE 6 — Does a unit's depth (Z) relate to its responses?
% Depth: um below the brain surface (getUnitDepths), source "coorZ" (manipulator coor_Z
% minus probe position). The LFP surface estimate only works on test recordings that cross
% the surface; it is reported as a check on coor_Z (FigS_LfpDepthCheck), not applied (decided
% 2026-11-03). The tests below rank depth within each insertion, so they do not depend on
% the surface position; only the absolute-depth axes and bins do.
% X/Y locations come later (brain clearing).
% Test (depthCorrTest): within-insertion rank correlation between depth and a
% measure; p from depths permuted within insertions; 95% CI from a hierarchical
% bootstrap (animal -> insertion -> unit). Benjamini-Hochberg across all tests;
% leave-one-animal-out for tests that survive.

d = figureDesign();                                               % shared design
paperFigPrefix('fig6');                                           % saved panels are named fig6_<panel>
depthSource = "coorZ";                                            % manipulator depth (LFP only validates it: FigS_LfpDepthCheck)
stims = ["MB", "RG", "SDGm", "SDGs", "MBR"];                      % stimuli in the paper
plotLab = @(s) replace(replace(replace(s, "RG", "SB"), "SDGm", "MG"), "SDGs", "SG");   % names shown on the panels (as in the other figures)
cacheF = fullfile(d.paperFigs, 'Fig6_depth_units.mat');           % per-unit table (cache)

%% Depths and per-unit measures (slow: cached)
D = getUnitDepths(d.exps, source=depthSource);                    % depth per good unit
M = table();
for ex = d.exps
    NP = loadNPclassFromTable(ex);                                % recording
    for s = stims
        if ismember(s, ["SDGm", "SDGs"]) && ~ismember(ex, d.expsSDG), continue, end   % recordings 1, 8: grating protocol excluded
        m = unitStimMeasures(NP, s);                              % empty if the stimulus was not recorded
        if isempty(m), continue, end
        m.exp = repmat(ex, height(m), 1);                         % recording / insertion ID
        M = [M; m];                                               %#ok<AGROW>
    end
end
nM = height(M);                                                   % rows before pairing with depths
M = innerjoin(M, D(:, {'exp','phyID','animal','depthUm','depthValid'}), 'Keys', {'exp','phyID'});   % pair by key, never by row
assert(height(M) == nM, 'Fig6: %d units have no depth row (join dropped them).', nM - height(M));   % innerjoin drops unmatched rows silently
assert(all(~isnan(M.depthUm)), 'Fig6: units without depth.');     % every unit must have a depth
save(cacheF, 'M', 'D', 'depthSource');

%% Tests
load(cacheF, 'M', 'D', 'depthSource');
nOld = height(M);                                                 % rows in the cache
M = M(~(ismember(M.stim, ["SDGm", "SDGs"]) & ~ismember(M.exp, d.expsSDG)), :);   % gratings only from d.expsSDG (cache may predate it)
fprintf('Fig6: %d grating rows of recordings 1, 8 dropped.\n', nOld - height(M));   % report what was removed
rng(1);                                                          % reproducible permutation p-values and bootstrap CIs
tests = table();                                                  % one row per test
for s = stims
    a = M(M.stim == s & M.present, :);                            % units present during this stimulus
    tests = runTest(tests, s + ": responsive (yes/no)", a, double(a.responsive));
    r = a(a.responsive, :);                                       % responsive units
    tests = runTest(tests, s + ": z", r, r.z);
    tests = runTest(tests, s + ": delta rate", r, r.dRate);
end
% MB vs SB preference: units present in both blocks, responsive to either (as Figure 1)
mb = M(M.stim == "MB", {'exp','phyID','animal','depthUm','depthValid','z','responsive','present'});
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
t2 = L2.tbl(L2.tbl.onOff == 1 & L2.tbl.sizeIdx == 1 & L2.tbl.lumValue == 255, :);   % white only (errors on a pre-fix cache)
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

%% Shared plotting settings
% Absolute-depth axes use only insertions whose depths are all below the surface
% (getUnitDepths depthValid; recs 1, 8 have coor_Z ~ -21). The tests above use every unit:
% within-insertion ranks do not depend on the surface reference.
anAll = unique(M.animal);                                         % animals, fixed order for every panel
cmap = animalColours(anAll);                                      % fixed colour per animal, as in every other figure
exclAn = unique(M.animal(~M.depthValid));                         % animals left out of absolute-depth axes
exclNote = "";                                                    % subtitle note on excluded animals
if ~isempty(exclAn), exclNote = "depth axis excludes " + join(exclAn, ", ") + " (invalid coor\_Z)"; end
binW = 500;                                                       % depth bin (um): coarser than the coor_Z uncertainty (300-700 um)
depthLim = [0, ceil(max(M.depthUm(M.depthValid)) / binW) * binW]; % common depth axis for every panel (um)

%% Panels: z vs depth per stimulus (responsive units), animal-colored
fig = figure('Position', [100 100 1300 340]); tl = tiledlayout(1, numel(stims), 'TileSpacing', 'compact');
for s = stims
    nexttile(tl); hold on;
    r = M(M.stim == s & M.present & M.responsive & M.depthValid, :); [~, ai] = ismember(r.animal, anAll);   % plotted units
    scatter(r.z, r.depthUm, 12, cmap(ai, :), 'filled', 'MarkerFaceAlpha', 0.6);
    set(gca, 'YDir', 'reverse'); ylim(depthLim); xlabel('z-score'); ylabel('Depth below surface (\mum)');
    row = tests(tests.name == s + ": z", :);                      % test on all units (ranks within insertion)
    title(plotLab(s));                                            % stimulus (stats below, so tiles do not overlap)
    if ~isempty(row) && ~isnan(row.rho)
        subtitle(sprintf('\\rho_w = %.2f [%.2f, %.2f]\np = %.3g', row.rho, row.ciLo, row.ciHi, row.p), 'FontSize', 8);
    else
        subtitle('too few units', 'FontSize', 8);
    end
end
subtitle(tl, exclNote);                                           % note on excluded animals
animalLegend(nexttile(tl, numel(stims)), anAll, cmap);            % legend in the last tile
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'z_vs_depth.jpg']), 'Resolution', 300);

%% Panels: per-insertion rank correlation of z with depth (MB, SB) — what the test pools
% One dot per insertion (Spearman rho of depth vs z, responsive units, >= 3 units), size ~ n.
% Line and band: pooled rho_w and its 95% CI. Valid for every insertion (no surface reference).
fig = figure('Position', [100 100 800 420]); tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
for s = ["MB", "RG"]
    nexttile(tl); hold on;
    r = M(M.stim == s & M.present & M.responsive, :);             % responsive units, all insertions
    [g, ins, an] = findgroups(r.exp, r.animal);                   % insertion groups
    nI = accumarray(g, 1);                                        % units per insertion
    rhoI = splitapply(@(dd, zz) corr(dd, zz, 'Type', 'Spearman'), r.depthUm, r.z, g);   % per-insertion rank correlation
    keep = nI >= 3 & ~isnan(rhoI);                                % insertions the test uses
    [~, o] = sortrows([double(categorical(an(keep), anAll)), ins(keep)]);   % order by animal, then insertion
    rk = rhoI(keep); nk = nI(keep); ak = an(keep); ik = ins(keep);
    rk = rk(o); nk = nk(o); ak = ak(o); ik = ik(o);               % sorted
    [~, ai] = ismember(ak, anAll);                                % animal colours
    row = tests(tests.name == s + ": z", :);                      % pooled test
    if ~isempty(row) && ~isnan(row.rho)
        patch([row.ciLo row.ciHi row.ciHi row.ciLo], [0.5 0.5 numel(rk)+0.5 numel(rk)+0.5], [0.85 0.85 0.85], ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');      % 95% CI of rho_w
        xline(row.rho, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');   % pooled rho_w
    end
    xline(0, 'k:', 'HandleVisibility', 'off');                    % no association
    scatter(rk, 1:numel(rk), 12 + 6 * nk, cmap(ai, :), 'filled', 'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');
    yticks(1:numel(rk)); yticklabels(compose('exp %d (n=%d)', ik, nk)); set(gca, 'YDir', 'reverse', 'FontSize', 7);
    xlim([-1 1]); ylim([0.5 numel(rk) + 0.5]); xlabel('Spearman \rho (depth vs z)');
    lab = s; if s == "RG", lab = "SB"; end                        % paper label
    if ~isempty(row) && ~isnan(row.rho)
        title(sprintf('%s: \\rho_w = %.2f [%.2f, %.2f], p = %.3g', lab, row.rho, row.ciLo, row.ciHi, row.p));
    else
        title(lab + ": too few units");
    end
end
animalLegend(nexttile(tl, 2), anAll, cmap);                       % legend
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'z_depth_rho_per_insertion.jpg']), 'Resolution', 300);

%% Panel: MB - SB z vs depth (units present in both, responsive to either)
fig = figure('Position', [100 100 460 420]); hold on;
Jv = J(J.depthValid, :); [~, ai] = ismember(Jv.animal, anAll);    % plotted units (valid depth axis)
xline(0, 'k:', 'HandleVisibility', 'off');                        % no preference
scatter(Jv.z - Jv.zRG, Jv.depthUm, 14, cmap(ai, :), 'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse'); ylim(depthLim); xlabel('z_{MB} - z_{SB}'); ylabel('Depth below surface (\mum)');
row = tests(tests.name == "MB - SB z (OR-set)", :);               % test on all units
if ~isempty(row) && ~isnan(row.rho)
    title(sprintf('MB - SB: \\rho_w = %.2f [%.2f, %.2f], p = %.3g', row.rho, row.ciLo, row.ciHi, row.p));
end
subtitle({sprintf('%d units plotted (test: %d)', height(Jv), row.n), char(exclNote)}, 'FontSize', 8);   % counts and exclusions, two lines
animalLegend(gca, anAll, cmap);                                   % legend
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'MBminusSB_z_vs_depth.jpg']), 'Resolution', 300);

%% Panels: fraction responsive by depth (500 um bins), per stimulus
fig = figure('Position', [100 100 1200 300]); tl = tiledlayout(1, numel(stims), 'TileSpacing', 'compact');
v = M(M.depthValid, :);                                           % units with a valid depth axis
edges = floor(min(v.depthUm)/binW)*binW : binW : ceil(max(v.depthUm)/binW)*binW;   % common depth bins
for s = stims
    nexttile(tl);
    a = v(v.stim == s & v.present, :); b = discretize(a.depthUm, edges);
    fr = accumarray(b, a.responsive, [numel(edges)-1 1], @mean, NaN); nb = accumarray(b, 1, [numel(edges)-1 1]);
    barh(edges(1:end-1) + binW/2, fr, 1, 'FaceColor', [0.6 0.6 0.6]); set(gca, 'YDir', 'reverse');   % bar at bin centre
    xMax = max([fr; 0.01], [], 'omitnan') * 1.45;                 % room for the n labels right of the bars
    xlim([0 xMax]); ylim(depthLim);                               % fixed axes: bars and labels stay inside
    xlabel('Fraction responsive'); ylabel('Depth below surface (\mum)'); title(plotLab(s));   % panel name
    for k = find(nb' > 0)
        text(max(fr(k), 0) + 0.03 * xMax, edges(k) + binW/2, sprintf('n=%d', nb(k)), 'FontSize', 7);   % units in the bin
    end
    row = tests(tests.name == s + ": responsive (yes/no)", :);
    if ~isempty(row) && ~isnan(row.rho), subtitle(sprintf('\\rho_w = %.2f, p = %.3g', row.rho, row.p)); end
end
subtitle(tl, exclNote);                                           % note on excluded animals
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

function animalLegend(ax, anAll, cmap)
% Animal colour legend on axes ax (dummy markers, so data handles need no names).
hold(ax, 'on');
h = gobjects(numel(anAll), 1);                                    % one handle per animal
for k = 1:numel(anAll)
    h(k) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', cmap(k, :), 'MarkerEdgeColor', 'none', ...
        'DisplayName', anAll(k));                                 % invisible marker carrying the colour
end
legend(ax, h, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'off');   % legend outside the data
end

function q = bhAdjust(p)
% Benjamini-Hochberg adjusted p-values.
[ps, o] = sort(p(:)); m = numel(ps);
qs = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));               % step-up adjustment
q = zeros(m, 1); q(o) = qs;
end
