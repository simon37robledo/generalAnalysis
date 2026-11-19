%% FIGURE 6 (position panels) — Do responses vary with atlas x, y, z position? (atlas animals: PV97, PV132, PV140, PV233)
% Uses the per-unit cache written by Fig6_Depth.m (Fig6_atlas_units.mat). Atlas coordinates are the
% TIFF-local axes of the brain scan registration, NOT stereotaxic AP/ML/DV; the axes are labelled
% atlas x/y/z. Within one recording a unit's position varies mainly along the probe (depth) and
% across shanks, so these tests are collinear with the depth tests: read them as "is there a
% gradient beyond what depth already shows", not as independent evidence.
% Test (depthCorrTest): within-recording rank correlation between the coordinate and a measure;
% p from the coordinate permuted within recordings; 95% CI from a hierarchical bootstrap
% (animal -> recording -> unit). Benjamini-Hochberg over the 12 tests (3 axes x MB/SB x responsive/z).
% Also: supplementary comparison of atlas depth with the manipulator (coor_Z) depth per recording.

d = figureDesign();                                               % shared design
paperFigPrefix('fig6');                                           % saved panels are named fig6_<panel>
L = load(fullfile(d.paperFigs, 'Fig6_atlas_units.mat'), 'M', 'D');   % unit x stimulus measures and unit table
M = L.M; D = L.D;                                                 % shortcuts
anAll = unique(M.animal);                                         % animals, fixed order
cmap = animalColours(anAll);                                      % fixed colour per animal
axes3 = ["atlasX_um", "atlasY_um", "atlasZ_um"];                  % coordinates tested
axLab = ["atlas x", "atlas y", "atlas z"];                        % labels (not AP/ML/DV)
rng(1);                                                          % reproducible permutations and bootstrap

%% Tests: coordinate vs responsiveness / z (MB, SB)
tests = table();                                                  % one row per test
for s = ["MB", "RG"]
    a = M(M.stim == s & M.present, :);                            % present units
    r = a(a.responsive, :);                                       % responsive units
    for k = 1:3
        tests = runCoordTest(tests, s + ": responsive vs " + axLab(k), a, a.(axes3(k)), double(a.responsive));
        tests = runCoordTest(tests, s + ": z vs " + axLab(k), r, r.(axes3(k)), r.z);
    end
end
ok = ~isnan(tests.p);                                             % tests that ran
q = nan(height(tests), 1); q(ok) = bhAdjust(tests.p(ok));         % Benjamini-Hochberg over this family
tests.qBH = q;
disp(tests);
writetable(tests, fullfile(d.paperFigs, 'Fig6_position_tests.csv'));

%% Panel: forest of the within-recording rank correlations
fig = figure('Position', [100 100 620 420]); hold on;
n = height(tests); yy = (1:n)';                                   % one row per test
for i = 1:n
    if isnan(tests.rho(i)), continue, end                         % test not run
    plot([tests.ciLo(i) tests.ciHi(i)], [i i], 'k-', 'LineWidth', 1.2);   % 95% CI
    plot(tests.rho(i), i, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);  % rho_w
    text(1.02, i, sprintf('p=%.3g, q=%.2f', tests.p(i), tests.qBH(i)), 'FontSize', 7, 'VerticalAlignment', 'middle');   % exact p and BH q
end
xline(0, 'k:'); set(gca, 'YDir', 'reverse', 'YTick', yy, 'YTickLabel', tests.name, 'FontSize', 7, 'TickLabelInterpreter', 'none');
xlim([-1 1]); ylim([0.5 n + 0.5]); xlabel('\rho_w (within-recording rank correlation) [95% CI]');
title('Response vs atlas position (axes are not AP/ML/DV)');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'position_rho_forest.jpg']), 'Resolution', 300);

%% Panel: unit maps per animal (MB z-score; grey = not responsive); descriptive
fig = figure('Position', [100 100 1500 700]); tl = tiledlayout(2, numel(anAll), 'TileSpacing', 'compact');
a = M(M.stim == "MB" & M.present, :);                             % MB units present
pairs = {"atlasX_um", "atlasY_um", 'atlas x (\mum)', 'atlas y (\mum)'; "atlasZ_um", "atlasY_um", 'atlas z (\mum)', 'atlas y (\mum)'};   % two projections
for p = 1:2
    for k = 1:numel(anAll)
        nexttile(tl); hold on;
        u = a(a.animal == anAll(k), :);                           % this animal's units
        g = ~u.responsive;                                        % non-responsive
        scatter(u.(pairs{p, 1})(g), u.(pairs{p, 2})(g), 8, [0.8 0.8 0.8], 'filled');   % non-responsive units (grey)
        scatter(u.(pairs{p, 1})(~g), u.(pairs{p, 2})(~g), 12, min(u.z(~g), 10), 'filled');   % responsive units coloured by z (clipped at 10)
        axis equal; clim([0 10]); colormap(parula);
        xlabel(pairs{p, 3}); ylabel(pairs{p, 4});
        if p == 1, title(sprintf('%s (n=%d, %d responsive)', anAll(k), height(u), sum(~g))); end
    end
end
cb = colorbar(nexttile(tl, numel(anAll)), 'eastoutside'); cb.Label.String = 'MB z-score';   % colour scale
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'position_unit_maps.jpg']), 'Resolution', 300);

%% Supplement: atlas depth vs manipulator (coor_Z) depth, per recording (validation of the depth source)
Lo = load(fullfile(d.paperFigs, 'Fig6_depth_units.mat'), 'D');    % old coor_Z depths (8-animal version)
Do = Lo.D(ismember(Lo.D.exp, d.expsAtlas), {'exp','phyID','depthUm'});   % atlas recordings only
Do.Properties.VariableNames{3} = 'depthCoorZ';                    % explicit name
V = innerjoin(D(D.depthValid, {'exp','phyID','animal','depthUm'}), Do, 'Keys', {'exp','phyID'});   % pair by key
fprintf('Supplement: %d of %d atlas units paired with a coor_Z depth.\n', height(V), sum(D.depthValid));
assert(height(V) == sum(D.depthValid), 'Fig6_Position: %d units lack a coor_Z depth row.', sum(D.depthValid) - height(V));   % no silent drops
[g, ex] = findgroups(V.exp);                                      % recordings
mA = splitapply(@mean, V.depthUm, g); mC = splitapply(@mean, V.depthCoorZ, g);   % mean depth per recording, both sources
anOf = splitapply(@(s) s(1), V.animal, g);                        % animal per recording
[~, ai] = ismember(anOf, anAll);
fig = figure('Position', [100 100 480 440]); hold on;
lim = [0 max([mA; mC]) * 1.05]; plot(lim, lim, 'k:', 'HandleVisibility', 'off');   % identity line
scatter(mC, mA, 40, cmap(ai, :), 'filled'); xlim(lim); ylim(lim); axis square;
xlabel('Manipulator depth: coor\_Z - probe y (\mum)'); ylabel('Atlas depth along shank (\mum)');
title(sprintf('Mean depth per recording (n=%d)', numel(ex)));   % short title (a long one is clipped)
subtitle(sprintf('atlas - manipulator: mean %.0f \\mum, SD %.0f \\mum', mean(mA - mC), std(mA - mC)), 'FontSize', 9);   % offset between the two sources
h = gobjects(numel(anAll), 1); for k = 1:numel(anAll), h(k) = plot(NaN, NaN, 'o', 'MarkerFaceColor', cmap(k, :), 'MarkerEdgeColor', 'none', 'DisplayName', anAll(k)); end
legend(h, 'Location', 'southeast', 'FontSize', 7, 'Box', 'off');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'depth_atlas_vs_coorZ.jpg']), 'Resolution', 300);

%% Local functions
function tests = runCoordTest(tests, name, tbl, coord, y)
% One coordinate test (depthCorrTest); NaN row if too few units/recordings.
ok = ~isnan(coord) & ~isnan(y);                                   % complete rows (a unit without a position cannot be ranked)
row = struct('name', string(name), 'rho', NaN, 'ciLo', NaN, 'ciHi', NaN, 'p', NaN, 'n', sum(ok), 'nIns', NaN, 'nAnimals', NaN);
try
    R = depthCorrTest(coord(ok), y(ok), tbl.exp(ok), tbl.animal(ok));   % within-recording rank correlation
    row.rho = R.rho; row.ciLo = R.ciLo; row.ciHi = R.ciHi; row.p = R.p; row.n = R.n; row.nIns = R.nIns; row.nAnimals = R.nAnimals;
catch ME
    if ~contains(ME.message, 'fewer than 2 usable insertions'), rethrow(ME); end   % only "too few data" is tolerated
    warning('Fig6Position:tooFew', '%s: too few units per recording to test.', name);   % explicit
end
tests = [tests; struct2table(row)];
end

function q = bhAdjust(p)
% Benjamini-Hochberg adjusted p-values.
[ps, o] = sort(p(:)); m = numel(ps);
qs = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));               % step-up adjustment
q = zeros(m, 1); q(o) = qs;
end
