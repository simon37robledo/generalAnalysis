%% SUPPLEMENTARY TABLE / FIGURE — LFP brain-surface estimate vs manipulator depth (coor_Z)
% The brain surface is located on the whole-shank / depth-test recordings, the only ones
% whose sites reach above the brain (estimateLfpSurface: topmost breakpoint of the band-power
% profile vs depth; port of the collaborator's notebook). Each test is matched to an insertion
% (getLfpSurface: insertion number in the folder name, else nearest recording in time) and
% compared with that insertion's coor_Z. This is a check on coor_Z only: the depths used in
% Figure 6 stay coor_Z. The Figure 6 statistics (within-insertion rank correlation,
% depthCorrTest) are unchanged by any per-insertion surface offset.
% Output: <Paper_figs>\Stats summary tables\LFP_surface_vs_coorZ.csv and sfig6_LfpSurface_vs_coorZ.jpg

d = figureDesign();                                               % shared design (recordings, output folder)
paperFigPrefix('sfig6');                                          % Supplementary Figure 6 panels: sfig6_<panel>
outCsv = fullfile(d.paperFigs, 'Stats summary tables', 'LFP_surface_vs_coorZ.csv');   % output table
X = readtable('\\sil3\data\Large_scale_mapping_NP\Experiment_Excel.xlsx');   % master table (animal per row)

%% One getLfpSurface call per animal (its test table covers all the animal's tests)
animals = unique(string(X.Animal_ID(d.exps)), 'stable');          % animals in the paper
T = table();                                                      % all tests
for a = animals'
    ex = d.exps(find(string(X.Animal_ID(d.exps)) == a, 1));       % any recording of this animal
    NP = loadNPclassFromTable(ex);                                % recording object
    obj = linearlyMovingBallAnalysis(NP, MultipleOffsets=false);  % path-only construction (any MB session)
    S = obj.getLfpSurface(verbose=false);                         % tests of this animal (cached estimates)
    if isempty(S.tests)                                           % no test recording: reported, not tabulated
        fprintf('%s: no whole-shank / depth test.\n', a);
        continue
    end
    t = S.tests;                                                  % this animal's tests
    t.animal = repmat(a, height(t), 1);                           % animal column
    T = [T; t];                                                   %#ok<AGROW>
end
T = movevars(T, 'animal', 'Before', 1);                           % animal first
writetable(T, outCsv);                                            % supplementary table
fprintf('Wrote %s\n', outCsv);

%% Summary for the methods (valid tests only: boundary inside the column, bands agree, coor_Z known)
v = T(T.valid, :);                                                % usable tests
fprintf('%d valid tests of %d (%d animals, %d insertions).\n', height(v), height(T), ...
    numel(unique(v.animal)), numel(unique(v.matchedExp)));        % counts
fprintf('LFP surface - coor_Z: median %.0f um, range %.0f to %.0f um, mean |delta| %.0f um.\n', ...
    median(v.deltaUm), min(v.deltaUm), max(v.deltaUm), mean(abs(v.deltaUm)));   % agreement
disp(T(:, {'animal','test','matchedInsertion','shank','lfpMeanUm','lfpSdUm','coorZ','deltaUm','valid'}));

%% Panel: LFP surface vs coor_Z (animal-coloured; invalid tests as open markers)
fig = figure('Position', [100 100 520 480]); hold on;
an = unique(T.animal, 'stable'); cmap = lines(numel(an));        % one colour per animal
lims = [min([T.lfpMeanUm; T.coorZ]) max([T.lfpMeanUm; T.coorZ])] + [-100 100];   % common axis range (um)
plot(lims, lims, 'k:', 'HandleVisibility', 'off');               % identity: LFP = coor_Z
for k = 1:numel(an)
    r = T(T.animal == an(k), :);                                  % tests of one animal
    ok = r.valid;                                                 % valid tests
    errorbar(r.coorZ(ok), r.lfpMeanUm(ok), r.lfpSdUm(ok), 'o', 'Color', cmap(k, :), ...
        'MarkerFaceColor', cmap(k, :), 'DisplayName', an(k));     % mean +- SD across bands
    plot(r.coorZ(~ok), r.lfpMeanUm(~ok), 'o', 'Color', cmap(k, :), 'HandleVisibility', 'off');   % excluded tests
end
xlim(lims); ylim(lims); axis square; box off;
xlabel('Manipulator depth coor\_Z (\mum)'); ylabel('LFP surface estimate (\mum from tip)');
title('Brain surface: LFP vs manipulator'); subtitle('filled = valid, open = excluded; bars = SD across bands');
legend('Location', 'southeast');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'LfpSurface_vs_coorZ.jpg']), 'Resolution', 300);
