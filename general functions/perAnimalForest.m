function T = perAnimalForest(structPaths, groupLabels, pair, metric, params)
% perAnimalForest  Per-animal paired difference (A - B) with a within-animal bootstrap CI,
% several analysis structs (e.g. Pogona and Laudakia) on one axis. Descriptive (2026-11-18).
%
%   T = perAnimalForest(structPaths, groupLabels, pair, metric, ...)
%
%   structPaths - string array of AllExpAnalysis *_analysisStruct.mat files (paired design,
%                 presence="allBlocks"; one per group, e.g. species).
%   groupLabels - string array, same length: label of each struct ("Pogona", "Laudakia").
%   pair        - {A, B} stimulus labels as stored in TableStimComp (e.g. {'MB','RG'}).
%   metric      - 'Z-score' or 'DeltaRate'.
%   Name-value: nBoot (10000), seed (0), outName ('' = no file), xLabel ('').
%
%   Per animal: units paired by (insertion, NeurID) with innerjoin (never by row position),
%   mean of A - B over its units, and a 95% CI from hierBootMatchFreq with that animal as
%   the only level-1 group (insertions, then units, resampled). Animal colours from
%   animalColours. Output table T: group, animal, nUnits, nInsertions, meanDiff, ciLo, ciHi.

arguments
    structPaths string
    groupLabels string
    pair cell
    metric char
    params.nBoot double = 10000
    params.seed double = 0
    params.outName char = ''
    params.xLabel char = ''
end
assert(numel(structPaths) == numel(groupLabels), 'perAnimalForest: one label per struct.');   % alignment
rng(params.seed, 'twister');                                         % reproducible CIs

rows = {};                                                           % accumulator of per-animal rows
for s = 1:numel(structPaths)
    assert(isfile(structPaths(s)), 'perAnimalForest: struct not found: %s', structPaths(s));   % fail loudly
    S = load(structPaths(s), 'TableStimComp', 'params');            % pooled table + design
    assert(~(isfield(S.params, 'presence') && S.params.presence == "perStimulus"), ...
        'perAnimalForest: %s holds unpaired per-block unit sets.', structPaths(s));   % paired design only
    Tc = S.TableStimComp;
    A = Tc(Tc.stimulus == pair{1}, {'animal', 'insertion', 'NeurID', metric});   % stimulus A rows
    B = Tc(Tc.stimulus == pair{2}, {'animal', 'insertion', 'NeurID', metric});   % stimulus B rows
    assert(height(A) > 0 && height(B) > 0, 'perAnimalForest: %s lacks %s or %s rows.', structPaths(s), pair{1}, pair{2});
    J = innerjoin(A, B, 'Keys', {'animal', 'insertion', 'NeurID'});  % pair units by key
    assert(height(J) == height(A) && height(J) == height(B), ...
        'perAnimalForest: %d paired of %d / %d rows in %s; unit sets differ.', height(J), height(A), height(B), structPaths(s));   % paired design invariant
    dd = J.([metric '_A']) - J.([metric '_B']);                       % per-unit difference A - B
    animals = unique(string(J.animal));                              % animals of this struct
    for a = animals'
        k = string(J.animal) == a;                                   % units of this animal
        ins = double(categorical(string(J.insertion(k))));           % insertion index within the animal
        b = hierBootMatchFreq(dd(k), params.nBoot, ones(sum(k), 1), ins);   % one animal: insertions, then units
        rows(end+1, :) = {groupLabels(s), a, sum(k), numel(unique(ins)), mean(dd(k)), ...
            prctile(b, 2.5), prctile(b, 97.5)};                      %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'group', 'animal', 'nUnits', 'nInsertions', 'meanDiff', 'ciLo', 'ciHi'});

% ---- plot: one row per animal, groups top to bottom, gap between groups ----
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 9 1.2 + 0.55 * (height(T) + numel(structPaths))]);
ax = axes(fig); hold(ax, 'on');
y = zeros(height(T), 1); yy = 0;                                     % y position per row
for s = 1:numel(structPaths)
    k = find(T.group == groupLabels(s));                            % rows of this group
    for r = k'
        yy = yy + 1; y(r) = yy;                                      % next row
    end
    yy = yy + 1;                                                     % gap after the group
end
cols = zeros(height(T), 3);                                          % fixed animal colours, per group (species
for s = 1:numel(structPaths)                                         % share palette slots, so never one call
    k = T.group == groupLabels(s);                                   % for both: avoids the mixed-species warning)
    cols(k, :) = animalColours(T.animal(k));
end
for r = 1:height(T)
    plot(ax, [T.ciLo(r) T.ciHi(r)], [y(r) y(r)], '-', 'Color', cols(r, :), 'LineWidth', 1.5);   % 95% CI
    plot(ax, T.meanDiff(r), y(r), 'o', 'MarkerFaceColor', cols(r, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 5);   % mean
end
xline(ax, 0, 'k-');                                                  % no difference
set(ax, 'YDir', 'reverse', 'YTick', y, 'YTickLabel', compose("%s (n=%d)", T.animal, T.nUnits), 'TickLabelInterpreter', 'none');   % animal and unit count
ylim(ax, [-0.4, yy - 0.3]);                                          % room for the first group header
xl = xlim(ax);                                                       % header x position = right edge
for s = 1:numel(structPaths)
    y0 = min(y(T.group == groupLabels(s)));                          % first row of this group
    text(ax, xl(2), y0 - 0.75, groupLabels(s), 'FontWeight', 'bold', 'FontSize', 8, ...
        'HorizontalAlignment', 'right', 'Interpreter', 'none');      % group (species) header
end
if isempty(params.xLabel)
    params.xLabel = sprintf('%s - %s, %s', pair{1}, pair{2}, metric);   % default axis label
end
xlabel(ax, params.xLabel);
set(ax, 'FontSize', 8, 'FontName', 'helvetica', 'TickDir', 'out', 'Box', 'off');   % house style
if ~isempty(params.outName)
    pf = 'W:\Large_scale_mapping_NP\Paper_figs';                     % outputs go to Paper_figs
    exportgraphics(fig, fullfile(pf, [paperFigPrefix() params.outName '.jpg']), 'Resolution', 300);   % panel
    writetable(T, fullfile(pf, [paperFigPrefix() params.outName '.csv']));                           % numbers
end
end
