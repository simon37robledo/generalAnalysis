%% SUPPLEMENTARY FIGURE — Eye movement during moving-ball (MB) and static-ball (SB) trials
% Supports Fig 1: if the eye moved during trials, MB vs SB response differences could reflect
% retinal image motion rather than the stimulus.
%
% Eye position: Elena's pupil-ellipse tables (W:\...\pupil_metrics) on the ephys clock rebuilt by
% buildVideoFrameTimes (the old videoTimeStampsSynced.mat mapping slipped after dropped frames).
% Units: pupil diameters (primary, calibration-free) and degrees = asind(shift / R), R = half the
% DLC eye-corner distance (lower bound for camera obliquity, upper bound if the eye opening is
% narrower than the globe; sensitivity R x 0.8 reported).
%
% Per trial (eyeTrialMetrics): max displacement from the window's start position in
%   stim  onset -> offset (MB ~2.4 s, SB ~0.5-1 s: not comparable between stimuli),
%   post  the first 300 ms after onset (same length for every stimulus: the comparison window),
%   base  the 300 ms before onset (spontaneous reference),
% plus saccade-like events and the trial's mean position offset from the block median.
% Positions are first smoothed with a 50 ms moving median (camera vibration; user 2026-11-19).
% Per insertion: medians over kept trials (>= 80% valid frames). Cross-animal inference:
% hierBootMatchFreq(data, nBoot, ANIMAL, INSERTION) on per-insertion values, exact p, 95% CI,
% leave-one-animal-out. MB tracking test: regression of trial eye velocity on ball velocity;
% gain = largest singular value (independent of camera orientation/mirror), compared with a
% within-insertion direction-shuffle null.
% Response vs eye movement (eyeResponseCorrelation): for each Figure 1 unit, Spearman rho between
% the trial response (rate - pre-trial baseline) and the trial's eye movement, both centred within
% stimulus condition and detrended over trials; mean rho vs 0 with the hierarchical bootstrap
% (animal -> insertion -> unit), MB - SB paired by (insertion, unit). Plus, across insertions,
% the MB - SB z difference vs the insertion's eye movement (weak: one value per insertion).
%
% Recordings: figureDesign d.exps with a complete pupil table and a valid sync (others are listed).

d = figureDesign();                                                         % shared design (recordings, output folder)
paperFigPrefix('sfig1');                                                    % panels saved as sfig1_<name>
nBoot = 100000;                                                             % bootstrap replicates (as the main tables)
nShuffle = 1000;                                                            % direction shuffles for the tracking null
rng(0);                                                                     % reproducible bootstrap and shuffles
halfCellDeg = 4.85;                                                         % half of the 9.7 deg grid cell / SB square (METHODS_BRIEF)
ballDeg = [4.9 9.7 19.4];                                                   % ball diameters used (deg)
degPerScreenPx = atand((33 / 1080) / 21.5);                                 % screen scale at the centre (0.081 deg/px, METHODS_BRIEF)

M = readtable("W:\Large_scale_mapping_NP\EyeTracking\manifest.csv", 'TextType', 'string');   % eye-video manifest
exps = intersect(d.exps, M.experimentNum)';                                 % paper recordings that have an eye video
fprintf('%d paper recordings have an eye video: %s\n', numel(exps), mat2str(exps));
F1 = load(d.structPath(d.exps, {'MB', 'RG'}), 'TableStimComp');             % Figure 1 units (MB-RG analysis struct)
TS = F1.TableStimComp;                                                      % one row per unit x stimulus (same units for MB and RG)
TS.realIns = str2double(string(TS.realInsertion));                          % within-animal insertion number (numeric categorical: never double())
TS.id = str2double(string(TS.NeurID));                                      % unit index among the recording's good units

resultsF = fullfile(d.paperFigs, 'sfig1_EyeMovement_results.mat');          % computed tables, for redrawing the panels
if ~exist('replotOnly', 'var'), replotOnly = false; end                     % set replotOnly = true before running to redraw from resultsF
if replotOnly
    load(resultsF, 'Sum', 'Trials', 'Rho', 'Still', 'St', 'InsCorr', 'G', 'P', 'U', 'notUsed', 'nBoot');   % previous run
else
%% Per recording and stimulus
Sum = table();                                                              % one row per insertion x stimulus
Rho = table();                                                              % one row per unit x stimulus: response vs eye movement
Still = table();                                                            % one row per unit x stimulus: split-half response, all vs still trials
Trials = cell(0, 1);                                                        % per-trial tables (for the CDF panel)
notUsed = strings(0, 1);                                                    % recordings left out, with the reason
for ex = exps
    try
        [E, calib] = getEyePosition(ex);                                    % eye position on the ephys clock
    catch err                                                               % no usable eye data: list it, never drop it silently
        notUsed(end+1, 1) = sprintf("%d: %s", ex, err.message); %#ok<SAGROW>
        warning('SFig1_EyeMovement:noEye', 'experimentNum %d left out: %s', ex, err.message);
        continue
    end
    F = buildVideoFrameTimes(ex);                                           % cached; for the timing uncertainty
    NP = loadNPclassFromTable(ex);                                          % recording object
    sbDurMs = NaN;                                                          % this recording's SB trial duration (for the matched MB window)
    for stim = ["SB" "MB"]                                                  % SB first: its duration sets the MB matched window
        try
            switch stim
                case "MB", vs = linearlyMovingBallAnalysis(NP);             % default: the 9x9 session
                case "SB", vs = rectGridAnalysis(NP);                       % static ball (rect grid)
            end
            if stim == "MB", [Tr, info] = eyeTrialMetrics(vs, E, calib, matchMs=sbDurMs);   % MB: also 0 -> SB-duration window
            else, [Tr, info] = eyeTrialMetrics(vs, E, calib); sbDurMs = median(Tr.durMs);   % SB: its own trial duration
            end
        catch err                                                           % stimulus absent or unsynced: listed
            notUsed(end+1, 1) = sprintf("%d %s: %s", ex, stim, err.message); %#ok<SAGROW>
            warning('SFig1_EyeMovement:noStim', 'experimentNum %d %s left out: %s', ex, stim, err.message);
            continue
        end
        n = height(Tr);                                                     % trials
        Tr.experimentNum = repmat(ex, n, 1); Tr.animal = repmat(calib.animal, n, 1); Tr.stim = repmat(stim, n, 1);   % keys
        Trials{end+1, 1} = Tr(:, {'experimentNum', 'animal', 'stim', 'post_dispDeg', 'post_kept'}); %#ok<SAGROW>   % columns common to MB and SB, for the CDF panel
        Sum = [Sum; summariseBlock(Tr, info, calib, F, stim, halfCellDeg, degPerScreenPx, nShuffle)]; %#ok<AGROW>
        u = TS(TS.animal == calib.animal & TS.realIns == calib.insertion & TS.stimulus == "MB", :);   % this insertion's Figure 1 units
        if isempty(u)                                                       % no Figure 1 unit here: nothing to correlate
            notUsed(end+1, 1) = sprintf("%d %s: no Figure 1 units (correlation only)", ex, stim); %#ok<SAGROW>
            continue
        end
        Ru = eyeResponseCorrelation(vs, Tr, u.id, stim);                    % per-unit rho, response vs eye movement
        Ru.experimentNum = repmat(ex, height(Ru), 1); Ru.animal = repmat(calib.animal, height(Ru), 1);   % keys
        zCol = TS.('Z-score');                                              % Figure 1 z-scores
        Ru.zMB = arrayfun(@(i) zCol(TS.animal == calib.animal & TS.realIns == calib.insertion & TS.id == i & TS.stimulus == "MB"), Ru.NeurID);   % MB z of the unit
        Ru.zSB = arrayfun(@(i) zCol(TS.animal == calib.animal & TS.realIns == calib.insertion & TS.id == i & TS.stimulus == "RG"), Ru.NeurID);   % SB (RG) z of the unit
        Rho = [Rho; Ru]; %#ok<AGROW>                                        % collect
        Du = stillTrialDelta(vs, Tr, u.id, stim);                           % split-half peak response: all / still / moving trials
        Du.experimentNum = repmat(ex, height(Du), 1); Du.animal = repmat(calib.animal, height(Du), 1);   % keys
        Still = [Still; Du]; %#ok<AGROW>                                    % collect
    end
end
assert(~isempty(Sum), 'SFig1_EyeMovement:noData', 'No recording has usable eye data yet.');
Trials = vertcat(Trials{:});                                                % all trials
okIns = Sum.fracKept >= 0.5;                                                % insertion kept if >= 50% of its trials are kept
for k = find(~okIns)'                                                       % list insertions with too few usable trials
    notUsed(end+1, 1) = sprintf("%d %s: only %.0f%% of trials with >= 80%% valid frames", Sum.experimentNum(k), Sum.stim(k), 100 * Sum.fracKept(k)); %#ok<SAGROW>
end
Sum.used = okIns;                                                           % flag (rows kept in the CSV)
fprintf('Left out:\n'); fprintf('  %s\n', notUsed);                         % explicit list

%% Statistics on per-insertion values
U = Sum(Sum.used, :);                                                       % usable insertion x stimulus rows
tests = {                                                                   % name, stimulus A, column A, stimulus B, column B
    'MB - SB, 300 ms after onset [pupil diam]',   "MB", 'medPostPupil',  "SB", 'medPostPupil'
    'MB - SB, 300 ms after onset [deg]',          "MB", 'medPostDeg',    "SB", 'medPostDeg'
    'MB - SB, SB-length window [pupil diam]',     "MB", 'medMatchPupil', "SB", 'medMatchPupil'
    'MB - SB, SB-length window [deg]',            "MB", 'medMatchDeg',   "SB", 'medMatchDeg'
    'MB - SB, saccade-like events [per min]',     "MB", 'eventsPerMin',  "SB", 'eventsPerMin'
    'MB - SB, trial-position SD [deg]',           "MB", 'posSDDeg',      "SB", 'posSDDeg'
    'MB 300 ms after - before onset [pupil diam]',"MB", 'medPostPupil',  "MB", 'medBasePupil'
    'SB 300 ms after - before onset [pupil diam]',"SB", 'medPostPupil',  "SB", 'medBasePupil'
    'MB tracking gain - shuffle null',            "MB", 'trackGain',     "MB", 'trackGainNull'};
St = table();                                                               % results
for i = 1:size(tests, 1)
    D = pairedValues(U, tests{i, 2}, tests{i, 3}, tests{i, 4}, tests{i, 5}); % per-insertion differences, paired by experimentNum
    St = [St; bootRow(tests{i, 1}, D, nBoot)]; %#ok<AGROW>                  % bootstrap + LOAO
end

% Response vs eye movement, per unit (Figure 1 units), hierarchical over animal -> insertion -> unit
Rho = Rho(ismember(Rho.experimentNum, U.experimentNum), :);                 % insertions with usable eye data
nNaN = nnz(isnan(Rho.rhoMovement));                                         % units without a defined rho (no spikes / too few trials)
if nNaN > 0, warning('SFig1_EyeMovement:rhoNaN', '%d unit x stimulus rows have no defined rho; excluded.', nNaN); end
for stim = ["MB" "SB"]
    for meas = ["rhoMovement" "rhoPosition"]
        r = Rho(Rho.stim == stim & isfinite(Rho.(meas)), :);                % this stimulus and measure
        D = table(r.experimentNum, r.animal, r.(meas), 'VariableNames', {'experimentNum', 'animal', 'diff'});   % rho vs 0
        St = [St; bootRow(sprintf('%s response vs eye %s: Spearman rho', stim, extractAfter(meas, "rho")), D, nBoot)]; %#ok<AGROW>
    end
end
A = Rho(Rho.stim == "MB", {'experimentNum', 'NeurID', 'animal', 'rhoMovement'}); A.Properties.VariableNames{4} = 'a';   % MB rho
B = Rho(Rho.stim == "SB", {'experimentNum', 'NeurID', 'rhoMovement'}); B.Properties.VariableNames{3} = 'b';              % SB rho
D = innerjoin(A, B, 'Keys', {'experimentNum', 'NeurID'}); D.diff = D.a - D.b; D = D(isfinite(D.diff), :);              % paired by (insertion, unit)
St = [St; bootRow('MB - SB rho(response, eye movement)', D, nBoot)];        % does MB depend more on eye movement than SB?

% Still trials: does MB > SB hold when the trials with the most eye movement are left out?
Still = Still(ismember(Still.experimentNum, U.experimentNum), :);           % insertions with usable eye data
A = Still(Still.stim == "MB", {'experimentNum', 'NeurID', 'animal', 'dAll', 'dStill', 'dMove'});   % MB split-half responses
B = Still(Still.stim == "SB", {'experimentNum', 'NeurID', 'dAll', 'dStill'});                      % SB split-half responses
B.Properties.VariableNames(3:4) = {'sbAll', 'sbStill'};                     % rename before the join
P = innerjoin(A, B, 'Keys', {'experimentNum', 'NeurID'});                   % paired by (insertion, unit)
P.diffAll = P.dAll - P.sbAll; P.diffStill = P.dStill - P.sbStill;           % MB - SB, all trials / still trials
stillTests = {'MB - SB split-half response, all trials [spk/s]', 'diffAll'; ...
              'MB - SB split-half response, still trials [spk/s]', 'diffStill'};
for i = 1:size(stillTests, 1)
    D = table(P.experimentNum, P.animal, P.(stillTests{i, 2}), 'VariableNames', {'experimentNum', 'animal', 'diff'});
    St = [St; bootRow(stillTests{i, 1}, D(isfinite(D.diff), :), nBoot)]; %#ok<AGROW>
end
D = table(P.experimentNum, P.animal, P.diffStill - P.diffAll, 'VariableNames', {'experimentNum', 'animal', 'diff'});
St = [St; bootRow('(MB - SB) still minus all trials [spk/s]', D(isfinite(D.diff), :), nBoot)];   % does excluding movement shrink the MB advantage?
D = table(P.experimentNum, P.animal, P.dStill - P.dMove, 'VariableNames', {'experimentNum', 'animal', 'diff'});
St = [St; bootRow('MB response, still minus moving trials [spk/s]', D(isfinite(D.diff), :), nBoot)];   % movement effect within MB

% Insertion level: is the MB - SB z difference larger where the eye moved more? (weak: few insertions)
G = groupsummary(Rho(Rho.stim == "MB", :), 'experimentNum', 'mean', {'zMB', 'zSB'});   % mean z per insertion
G.zDiff = G.mean_zMB - G.mean_zSB;                                          % MB - SB z
G = innerjoin(G, U(U.stim == "MB", {'experimentNum', 'animal', 'medStimPupil', 'posSDPupil', 'eventsPerMin'}), 'Keys', 'experimentNum');   % eye measures
InsCorr = table();                                                          % across-insertion Spearman correlations
for meas = ["medStimPupil" "posSDPupil" "eventsPerMin"]
    [rr, pp] = corr(G.zDiff, G.(meas), 'Type', 'Spearman', 'Rows', 'complete');   % across insertions (no hierarchy: one value per insertion)
    InsCorr = [InsCorr; table(meas, height(G), rr, pp, 'VariableNames', {'eyeMeasure', 'nInsertions', 'rho', 'p'})]; %#ok<AGROW>
end
disp(St); disp(InsCorr);

%% Tables
tblDir = fullfile(d.paperFigs, 'Stats summary tables');                     % supplementary tables
writetable(Sum, fullfile(tblDir, 'Eye_movement_per_insertion.csv'));        % per insertion x stimulus
writetable(St, fullfile(tblDir, 'Eye_movement_stats.csv'));                 % tests
writetable(table(notUsed, 'VariableNames', {'left_out'}), fullfile(tblDir, 'Eye_movement_left_out.csv'));   % exclusions
writetable(Rho, fullfile(tblDir, 'Eye_movement_vs_response_per_unit.csv'));   % per-unit correlations
writetable(P, fullfile(tblDir, 'Eye_still_trials_per_unit.csv'));           % per-unit split-half responses, all vs still trials
writetable(InsCorr, fullfile(tblDir, 'Eye_movement_vs_zdiff_insertions.csv'));   % insertion-level correlations
save(resultsF, 'Sum', 'Trials', 'Rho', 'Still', 'St', 'InsCorr', 'G', 'P', 'U', 'notUsed', 'nBoot');   % for replotOnly
end                                                                         % end of the computation (skipped with replotOnly)

%% Panels
animals = unique(U.animal);                                                 % animals present
cmap = animalColours(animals);                                              % fixed animal colours
fig = figure('Position', [100 100 1500 840]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact');

% (A) per-insertion median displacement in the 300 ms windows and the full trial
nexttile(tl); hold on;
cols = {'SB','medBaseDeg'; 'SB','medPostDeg'; 'SB','medStimDeg'; 'MB','medBaseDeg'; 'MB','medPostDeg'; 'MB','medStimDeg'};   % x positions
lab = {'SB before','SB 0-300','SB trial','MB before','MB 0-300','MB trial'};  % x labels
for e = unique(U.experimentNum)'
    yv = nan(1, size(cols, 1));                                             % this insertion's values
    for c = 1:size(cols, 1)
        r = U(U.experimentNum == e & U.stim == cols{c, 1}, :);              % its row for the stimulus
        if height(r) == 1, yv(c) = r.(cols{c, 2}); end
    end
    a = U.animal(find(U.experimentNum == e, 1));                            % its animal
    col = cmap(animals == a, :);                                            % animal colour
    plot(1:3, yv(1:3), '-', 'Color', [col 0.4]); plot(4:6, yv(4:6), '-', 'Color', [col 0.4]);   % same insertion joined
    scatter((1:6) + 0.08 * randn(1, 6), yv, 25, col, 'filled', 'MarkerFaceAlpha', 0.8);   % points
end
yline(halfCellDeg, 'k--', 'half grid cell (4.85 deg)', 'LabelHorizontalAlignment', 'left');   % reference (smallest ball, 4.9 deg, is the same height)
hA = gobjects(numel(animals), 1);                                           % legend handles, one per animal
for a = 1:numel(animals), hA(a) = scatter(nan, nan, 25, cmap(a, :), 'filled'); end   % invisible markers carrying the animal colours
legend(hA, cellstr(animals), 'Location', 'southeast', 'Box', 'off');       % animal legend
set(gca, 'YScale', 'log', 'XTick', 1:6, 'XTickLabel', lab); xtickangle(30);
ylabel('Median max. eye displacement per trial (deg)'); title('Per insertion (animal colours)');

% (B) per-trial displacement in the 300 ms after onset, cumulative per insertion
nexttile(tl); hold on;
for e = unique(U.experimentNum)'
    a = U.animal(find(U.experimentNum == e, 1)); col = cmap(animals == a, :);   % animal colour
    for stim = ["MB" "SB"]
        x = Trials.post_dispDeg(Trials.experimentNum == e & Trials.stim == stim & Trials.post_kept);   % kept trials
        x = x(isfinite(x));                                                 % no degrees without an eye-corner calibration
        if isempty(x), continue, end
        [fy, fx] = ecdf(x);                                                 % empirical CDF
        plot(fx, fy, 'Color', col, 'LineStyle', ternary(stim == "MB", '-', '--'));   % MB solid, SB dashed
    end
end
xline(halfCellDeg, 'k--'); xline(1, 'k:');                                  % references
set(gca, 'XScale', 'log'); xlim([1e-3 10]); xlabel('Max. eye displacement in 0-300 ms (deg)');   % below 0.001 deg is below the jitter floor
ylabel('Fraction of trials'); title('MB solid, SB dashed');

% (C) MB tracking: gain vs shuffle null
nexttile(tl); hold on;
r = U(U.stim == "MB", :);                                                   % MB rows
for k = 1:height(r)
    col = cmap(animals == r.animal(k), :);                                  % animal colour
    plot([1 2], [r.trackGainNull(k) r.trackGain(k)], '-', 'Color', [col 0.5]);   % same insertion
    scatter([1 2], [r.trackGainNull(k) r.trackGain(k)], 25, col, 'filled');
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'Shuffled directions', 'Observed'}); xlim([0.5 2.5]);
ylabel('Eye velocity / ball velocity (gain)'); title('Does the eye follow the ball?');
s = St(St.test == "MB tracking gain - shuffle null", :);                    % its test
if height(s) == 1, [pt, op] = formatPValue(s.p, 2 / nBoot); text(1.5, max(ylim) * 0.95, sprintf('p %s %s', op, pt), 'HorizontalAlignment', 'center'); end

% (D) per-unit correlation of the trial response with the trial's eye movement
nexttile(tl); hold on;
for si = 1:2
    stim = ["MB" "SB"]; stim = stim(si);                                    % stimulus of this column
    r = Rho(Rho.stim == stim & isfinite(Rho.rhoMovement), :);               % its units
    for a = 1:numel(animals)
        k = r.animal == animals(a);                                         % units of this animal
        scatter(si + 0.15 * randn(nnz(k), 1), r.rhoMovement(k), 12, cmap(a, :), 'filled', 'MarkerFaceAlpha', 0.5);   % one point per unit
    end
    s = St(St.test == stim + " response vs eye Movement: Spearman rho", :);   % its bootstrap
    if height(s) == 1
        errorbar(si + 0.35, s.mean, s.mean - s.ciLow, s.ciHigh - s.mean, 'k', 'LineWidth', 1.5, 'CapSize', 6);   % mean and 95% CI
        [pt, op] = formatPValue(s.p, 2 / nBoot);                            % p text
        text(si, 0.45, sprintf('%.3f [%.3f, %.3f]\np %s %s', s.mean, s.ciLow, s.ciHigh, op, pt), 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
yline(0, 'k:'); xlim([0.4 2.6]); ylim([-0.5 0.55]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'MB', 'SB'});
ylabel('Spearman \rho (response, eye movement)'); title('Within condition, detrended; Fig 1 units');

% (E) insertion level: MB - SB z difference vs how much the eye moved during MB
nexttile(tl); hold on;
for k = 1:height(G)
    scatter(G.medStimPupil(k), G.zDiff(k), 30, cmap(animals == G.animal(k), :), 'filled');   % one point per insertion
end
xlabel('Median eye displacement per MB trial (pupil diameters)'); ylabel('Mean z, MB - SB (Fig 1 units)');
c = InsCorr(InsCorr.eyeMeasure == "medStimPupil", :);                       % its correlation
title(sprintf('Per insertion: \\rho = %.2f, p = %.3g (n = %d)', c.rho, c.p, c.nInsertions));

% (F) MB - SB split-half response: all trials vs still trials (eye movement below the condition median)
nexttile(tl); hold on;
for a = 1:numel(animals)
    k = P.animal == animals(a);                                             % this animal's units
    scatter(P.diffAll(k), P.diffStill(k), 12, cmap(a, :), 'filled', 'MarkerFaceAlpha', 0.5);   % one point per unit
end
lim = prctile([P.diffAll; P.diffStill], [1 99]);                            % robust axis range
plot(lim, lim, 'k:'); xline(0, 'k:'); yline(0, 'k:'); axis([lim lim]); axis square;
xlabel('MB - SB response, all trials (spk/s)'); ylabel('MB - SB response, still trials (spk/s)');
sA = St(St.test == "MB - SB split-half response, all trials [spk/s]", :);   % its tests
sS = St(St.test == "MB - SB split-half response, still trials [spk/s]", :);
if height(sA) == 1 && height(sS) == 1
    title(sprintf('All %.2f [%.2f, %.2f]; still %.2f [%.2f, %.2f]', sA.mean, sA.ciLow, sA.ciHigh, sS.mean, sS.ciLow, sS.ciHigh), 'FontSize', 9);
end

exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'EyeMovement.jpg']), 'Resolution', 300);   % raster
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'EyeMovement.pdf']), 'ContentType', 'vector');   % vector
fprintf('Saved %s\n', fullfile(d.paperFigs, [paperFigPrefix() 'EyeMovement.jpg']));
axs = findobj(tl, 'Type', 'axes');                                          % one axes per panel
for a = axs'                                                                % single panels for the slides
    exportgraphics(a, fullfile(d.paperFigs, sprintf('%seye_%c.jpg', paperFigPrefix(), char('A' + a.Layout.Tile - 1))), 'Resolution', 300);   % sfig1_eye_A ... E
end

%% Local functions
function row = summariseBlock(Tr, info, calib, F, stim, halfCellDeg, degPerScreenPx, nShuffle)
% Per-insertion summary of one stimulus block.
k = Tr.stim_kept;                                                           % trials with >= 80% valid frames
toDeg = @(px, R) asind(min(px ./ R, 1)) + 0 * R;                            % px -> deg; NaN when R is NaN (min(NaN, 1) = 1 gave 90 deg)
row = table(Tr.experimentNum(1), calib.animal, stim, 'VariableNames', {'experimentNum', 'animal', 'stim'});   % key
row.nTrials = height(Tr); row.nKept = nnz(k); row.fracKept = mean(k);       % coverage
row.medStimDeg = median(Tr.stim_dispDeg(k));                                % full trial
row.p95StimDeg = prctile(Tr.stim_dispDeg(k), 95);
row.medPostDeg = median(Tr.post_dispDeg(Tr.post_kept));                     % 300 ms after onset
row.p95PostDeg = prctile(Tr.post_dispDeg(Tr.post_kept), 95);
row.medBaseDeg = median(Tr.base_dispDeg(Tr.base_kept));                     % 300 ms before onset
row.medPostPupil = median(Tr.post_dispPupil(Tr.post_kept));                 % same in pupil diameters
row.medBasePupil = median(Tr.base_dispPupil(Tr.base_kept));
row.medStimPupil = median(Tr.stim_dispPupil(k));
row.medMatchPupil = median(Tr.stim_dispPupil(k));                           % SB: the whole SB trial is the matched window
row.medMatchDeg = median(Tr.stim_dispDeg(k));
if ismember('match_dispPupil', Tr.Properties.VariableNames)                 % MB: 0 -> SB duration after onset
    row.medMatchPupil = median(Tr.match_dispPupil(Tr.match_kept));
    row.medMatchDeg = median(Tr.match_dispDeg(Tr.match_kept));
end
noR = 0 * calib.RPx;                                                        % NaN when there is no degree calibration (NaN > x is false, not NaN)
row.fracStimGt1Deg = mean(Tr.stim_dispDeg(k) > 1) + noR;                    % trials moving > 1 deg
row.fracStimGtHalfCell = mean(Tr.stim_dispDeg(k) > halfCellDeg) + noR;      % trials moving > half a grid cell
row.fracStimGtHalfCell_R08 = mean(toDeg(Tr.stim_dispPx(k), 0.8 * calib.RPx) > halfCellDeg) + noR;   % same with R x 0.8 (larger angles)
row.eventsPerMin = sum(Tr.stim_nEvents(k)) / (sum(Tr.durMs(k)) / 60000);    % saccade-like events per minute of stimulus (NaN without calibration)
row.baseEventsPerMin = sum(Tr.base_nEvents(Tr.base_kept)) / (nnz(Tr.base_kept) * 0.3 / 60);   % ... in the 300 ms baselines
posPx = sqrt(var(Tr.stim_meanXPx(k)) + var(Tr.stim_meanYPx(k)));            % SD of trial mean positions [px]
row.posSDDeg = toDeg(posPx, calib.RPx); row.posSDPupil = posPx / calib.pupilDiamPx;   % across-trial position spread
row.RPx = calib.RPx; row.pupilDiamPx = calib.pupilDiamPx;                   % calibration
row.noiseFloorDeg = toDeg(calib.noiseFloorPx, calib.RPx);                   % raw frame-to-frame jitter
row.smoothNoiseDeg = toDeg(info.smoothNoisePx, calib.RPx);                  % jitter after the 50 ms moving median
row.smoothMs = info.smoothMs;                                               % smoothing window used
row.fracValidFrames = calib.fracValid; row.timeUncertaintyMs = F.timeUncertaintyMs;   % data quality
row.trackGain = NaN; row.trackGainNull = NaN; row.trackP = NaN;             % MB only
if stim == "MB" && ismember('direction', Tr.Properties.VariableNames)
    B = [cos(Tr.direction(k)) sin(Tr.direction(k))] .* (Tr.speedPxS(k) * degPerScreenPx);   % ball velocity [deg/s]
    V = [Tr.stim_velXPxS(k) Tr.stim_velYPxS(k)] * (180 / pi) / calib.RPx;  % eye velocity [deg/s] (small angles)
    ok = all(isfinite([B V]), 2); B = B(ok, :); V = V(ok, :);               % complete trials
    g = @(Bm) max(svd(Bm \ V));                                             % gain: largest singular value of the 2x2 map
    row.trackGain = g(B);                                                   % observed
    null = zeros(nShuffle, 1);                                              % shuffled directions
    for s = 1:nShuffle, null(s) = g(B(randperm(size(B, 1)), :)); end        % break direction-velocity pairing
    row.trackGainNull = mean(null); row.trackP = mean(null >= row.trackGain);   % null mean, one-tailed p within insertion
end
row.noiseThrPx = info.noiseThrPx;                                           % saccade threshold used
end

function D = pairedValues(U, sA, cA, sB, cB)
% Per-insertion difference A - B, paired by experimentNum (never by row position).
A = U(U.stim == sA, {'experimentNum', 'animal', cA}); A.Properties.VariableNames{3} = 'a';
B = U(U.stim == sB, {'experimentNum', cB}); B.Properties.VariableNames{2} = 'b';
D = innerjoin(A, B, 'Keys', 'experimentNum');                               % insertions with both
D.diff = D.a - D.b;                                                         % difference
D = D(isfinite(D.diff), :);                                                 % drop undefined (e.g. no calibration), counted below
end

function row = bootRow(name, D, nBoot)
% Hierarchical bootstrap (animal -> insertion) of the mean difference, plus leave-one-animal-out.
row = table(string(name), height(D), numel(unique(D.animal)), 'VariableNames', {'test', 'nInsertions', 'nAnimals'});
row.mean = mean(D.diff); row.ciLow = NaN; row.ciHigh = NaN; row.p = NaN; row.loaoP = "";
if row.nAnimals < 2, return, end                                            % no cross-animal inference with one animal
animalId = findgroups(D.animal);                                            % level 1: animal
insId = D.experimentNum;                                                    % level 2: globally unique insertion code
bm = hierBootMatchFreq(D.diff, nBoot, animalId, insId);                     % (data, nrep, ANIMAL, INSERTION)
row.ciLow = prctile(bm, 2.5); row.ciHigh = prctile(bm, 97.5);               % 95% CI
row.p = min(1, 2 * min(mean(bm <= 0), mean(bm >= 0)));                      % two-tailed p
lo = strings(0);                                                            % leave-one-animal-out p values
for a = unique(D.animal)'
    Dk = D(D.animal ~= a, :);                                               % without animal a
    if numel(unique(Dk.animal)) < 2, lo(end+1) = a + ":NaN"; continue, end  %#ok<AGROW>
    b = hierBootMatchFreq(Dk.diff, 10000, findgroups(Dk.animal), Dk.experimentNum);   % 10k as the LOAO tables
    lo(end+1) = sprintf("%s:%.3g", a, min(1, 2 * min(mean(b <= 0), mean(b >= 0)))); %#ok<AGROW>
end
row.loaoP = strjoin(lo, "; ");                                              % p without each animal
end

function out = ternary(c, a, b)
% Inline if.
if c, out = a; else, out = b; end
end
