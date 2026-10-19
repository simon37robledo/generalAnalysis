%% SUPPLEMENTARY FIGURE — Null calibration of every responsiveness test (artificial neurons)
% Simulated units with no stimulus response (homogeneous Bernoulli spiking at the
% same rate before and during the stimulus) run through each test design on the
% trial layout of a real recording. A valid test flags ~5% of them at alpha = 0.05.
% Reads the tables written by makeStatsSummaryTables (no re-simulation):
%   Null_calibration_false_positive_rate.csv  (MB and SB designs, incl. discarded ones)
%   Null_calibration_other_stimuli.csv        (FFF, SDGs, SDGm, MBR, NI)
%   Null_calibration_delta_rate_bias.csv      (split-half delta firing rate)

d = figureDesign();                                               % shared design (output folder)
paperFigPrefix('figS');                                           % saved panels are named figS_<panel>
tblDir = fullfile(d.paperFigs, 'Stats summary tables');           % supplementary tables
F1 = readTitledCsv(fullfile(tblDir, 'Null_calibration_false_positive_rate.csv'));   % MB / SB designs
F2 = readTitledCsv(fullfile(tblDir, 'Null_calibration_other_stimuli.csv'));         % other stimuli
B  = readTitledCsv(fullfile(tblDir, 'Null_calibration_delta_rate_bias.csv'));       % delta-rate bias

fig = figure('Position', [100 100 1350 380]);
tl  = tiledlayout(1, 3, 'TileSpacing', 'compact');
% (A) MB and SB designs: used = solid, discarded = dashed
nexttile(tl); fprPanel(F1, 'MB and SB tests');
% (B) other stimuli
nexttile(tl); fprPanel(F2, 'Other stimuli');
% (C) split-half delta firing rate under H0
nexttile(tl); hold on;
names = unique(B.(1), 'stable'); cmap = lines(numel(names));
for k = 1:numel(names)
    r = B(strcmp(B.(1), names{k}), :);                            % one selection space
    errorbar(r.(2), r.(4), r.(5), '-o', 'Color', cmap(k,:), 'MarkerFaceColor', cmap(k,:), 'MarkerSize', 4, 'DisplayName', names{k});
end
yline(0, 'k:', 'HandleVisibility', 'off');                        % unbiased estimator
set(gca, 'XScale', 'log'); xlabel('Simulated firing rate (Hz)'); ylabel('Mean \Delta rate (spikes/s) \pm SEM');
title('Split-half \Delta rate under H0'); legend('Location', 'northwest', 'FontSize', 7, 'Interpreter', 'none');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'NullCalibration.jpg']), 'Resolution', 300);

% Pooled false-positive rate per design (printed for the legend / methods)
for T = {F1, F2}
    t = T{1};
    for nm = unique(t.(1), 'stable')'
        r = t(strcmp(t.(1), nm{1}), :);
        fprintf('%-75s pooled FPR %.3f (%d / %d)\n', nm{1}, sum(r.(4)) / sum(r.(3)), sum(r.(4)), sum(r.(3)));
    end
end

%% Local functions
function T = readTitledCsv(f)
% Read a table written by writeTitledCsv: title, description, blank line, header,
% rows, then free-text notes (dropped: their numeric columns are empty).
T = readtable(f, 'NumHeaderLines', 3, 'VariableNamingRule', 'preserve', 'TextType', 'char');
T = T(~isnan(T.(2)), :);                                          % data rows only (notes have no rate)
end

function fprPanel(T, ttl)
% False-positive rate vs firing rate per design, Clopper-Pearson 95% CI bars.
hold on;
names = unique(T.(1), 'stable'); cmap = lines(numel(names));
for k = 1:numel(names)
    r  = T(strcmp(T.(1), names{k}), :);                           % one design
    ci = cellfun(@(s) sscanf(s, '[%f, %f]')', r.(6), 'UniformOutput', false);   % '[lo, hi]' -> numbers
    ci = vertcat(ci{:});
    ls = '-'; if contains(names{k}, 'not used'), ls = '--'; end   % discarded designs dashed
    errorbar(r.(2), r.(5), r.(5) - ci(:,1), ci(:,2) - r.(5), ls, 'Color', cmap(k,:), 'Marker', 'o', ...
        'MarkerFaceColor', cmap(k,:), 'MarkerSize', 4, 'DisplayName', names{k});
end
yline(0.05, 'k--', '\alpha = 0.05', 'HandleVisibility', 'off');   % nominal rate
set(gca, 'XScale', 'log', 'YScale', 'log'); ylim([0.003 1.2]);    % log scale: shows 5% and ~100% together
xlabel('Simulated firing rate (Hz)'); ylabel('False-positive rate');
title(ttl); legend('Location', 'southoutside', 'FontSize', 7, 'Interpreter', 'none');
end
