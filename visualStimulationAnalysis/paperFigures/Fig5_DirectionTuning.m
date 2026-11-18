%% FIGURE 5 — Direction tuning: moving ball (MB), drifting grating (SDGm), moving bar (MBR)
% Per-unit tuning from directionTuningUnits (same windows as Figures 1-4; MB at the
% RF-cell crossing). Units: responsive (p < 0.05) AND present, per stimulus.
% Direction selectivity: DSI with a direction-label permutation test per unit.
% Same tuning across stimuli: preferred-direction agreement in units selective to
% both stimuli, against a null that re-pairs units within the same insertion.
% Cross-animal statistics: hierBootMatchFreq (animal, insertion), leave-one-animal-out.

d = figureDesign();                                               % shared design
paperFigPrefix('fig5');                                           % saved panels are named fig5_<panel>
stims  = ["MB", "SDGm", "MBR"];                                   % compared stimuli
stimLab = replace(stims, "SDGm", "MG");                           % names shown on the panels (as in the other figures)
ctor   = {@linearlyMovingBallAnalysis, @StaticDriftingGratingAnalysis, @linearlyMovingBarAnalysis};  % analysis classes
cacheF = fullfile(d.paperFigs, 'Fig5_directionTuning_units.mat'); % per-unit table (cache)
nBoot  = 10000;                                                   % bootstrap / permutation size

%% Per-unit tuning for every recording (slow: cached)
U = table();                                                      % one row per unit x stimulus
for ex = d.exps
    NP = loadNPclassFromTable(ex);                                % recording
    animal = string(regexp(NP.recordingName, '(PV|SA)\d+', 'match', 'once'));   % animal ID
    for s = 1:numel(stims)
        if stims(s) == "SDGm" && ~ismember(ex, d.expsSDG)         % recordings 1 and 8: grating protocol excluded
            continue
        end
        try
            o = ctor{s}(NP);                                      % stimulus object (errors if absent)
            o.ResponseWindow;                                     % cached response window must exist
        catch
            continue                                              % stimulus not recorded here
        end
        Ts = directionTuningUnits(o, stims(s), nPerm=2000);       % per-unit tuning
        Ts.exp = repmat(ex, height(Ts), 1); Ts.animal = repmat(animal, height(Ts), 1);   % identifiers
        U = [U; Ts];                                              %#ok<AGROW>
    end
end
save(cacheF, 'U');                                                % cache for the panels below

%% Fraction of direction-selective units (per insertion), vs 5% chance
load(cacheF, 'U');
assert(~any(U.stim == "SDGm" & ~ismember(U.exp, d.expsSDG)), ...
    'Fig5: cache holds SDGm units from recordings outside d.expsSDG (1, 8): rebuild the cache.');   % stale cache
A = U(U.responsive & U.present & ~isnan(U.DSI), :);               % analysable units
rng(0);
fr = table();                                                     % one row per stimulus
for s = stims
    a = A(A.stim == s, :);                                        % this stimulus
    [g, ex, an] = findgroups(a.exp, a.animal);                    % insertions
    f = splitapply(@(p) mean(p < 0.05), a.pDSI, g);              % fraction DS per insertion
    b = hierBootMatchFreq(f, nBoot, double(categorical(an)), ex); % animal, then insertion
    fr = [fr; table(s, height(a), sum(a.pDSI < 0.05), numel(unique(an)), mean(b), prctile(b, 2.5), prctile(b, 97.5), ...
        mean(b <= 0.05), 'VariableNames', {'stim','nUnits','nDS','nAnimals','fracDS','ciLo','ciHi','pVsChance'})];   %#ok<AGROW>
end
disp(fr);
fig = figure('Position', [100 100 420 360]); hold on;
bar(1:numel(stims), fr.fracDS, 0.6, 'FaceColor', [0.7 0.7 0.7]);
errorbar(1:numel(stims), fr.fracDS, fr.fracDS - fr.ciLo, fr.ciHi - fr.fracDS, 'k.', 'LineWidth', 1.2);
yline(0.05, 'r--', 'chance');                                     % expected under H0
set(gca, 'XTick', 1:numel(stims), 'XTickLabel', stimLab); ylabel('Fraction direction-selective');   % panel names
title('Direction-selective units (p < 0.05)');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'DirectionSelective_fraction.jpg']), 'Resolution', 300);
writetable(fr, fullfile(d.paperFigs, 'Fig5_fraction_direction_selective.csv'));

%% DSI above chance per unit (animal-colored), bootstrap mean and 95% CI
% DSI minus the unit's shuffled-label mean DSI (DSIcorr): raw DSI is biased upward
% with few trials (an untuned unit scores ~0.3), so chance sits at 0 here.
% Units selective by the permutation test (pDSI < 0.05) are outlined in black.
fig = figure('Position', [100 100 570 380]); hold on;             % wider than before: room for the animal legend
anAll = unique(A.animal); cmap = animalColours(anAll);            % fixed color per animal, as in every other figure
for k = 1:numel(stims)
    a = A(A.stim == stims(k), :);
    [~, ai] = ismember(a.animal, anAll);
    sel = a.pDSI < 0.05;                                          % direction-selective units
    swarmchart(k * ones(sum(~sel), 1), a.DSIcorr(~sel), 14, cmap(ai(~sel), :), 'filled', 'MarkerFaceAlpha', 0.5, 'XJitterWidth', 0.5);   % not selective
    swarmchart(k * ones(sum(sel), 1), a.DSIcorr(sel), 18, cmap(ai(sel), :), 'filled', 'MarkerEdgeColor', 'k', 'XJitterWidth', 0.5);       % selective
    b = hierBootMatchFreq(a.DSIcorr, nBoot, double(categorical(a.animal)), a.exp);   % animal, then insertion
    errorbar(k, mean(b), mean(b) - prctile(b, 2.5), prctile(b, 97.5) - mean(b), 'k', 'LineWidth', 1.5, 'CapSize', 8);
end
yline(0, 'k:');                                                   % chance level
hAn = gobjects(numel(anAll), 1);                                  % one legend entry per animal
for k = 1:numel(anAll)
    hAn(k) = plot(NaN, NaN, 'o', 'MarkerFaceColor', cmap(k, :), 'MarkerEdgeColor', 'none', 'DisplayName', anAll(k));   % invisible marker carrying the colour
end
legend(hAn, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'off');   % animal colours, outside the data
set(gca, 'XTick', 1:numel(stims), 'XTickLabel', stimLab); ylabel('DSI - chance DSI');   % panel names
title('Direction selectivity above chance (black outline: p < 0.05)');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'DSI_swarm.jpg']), 'Resolution', 300);

%% Preferred directions of direction-selective units (screen directions)
fig = figure('Position', [100 100 900 300]);
tl = tiledlayout(1, numel(stims));
for k = 1:numel(stims)
    a = A(A.stim == stims(k) & A.pDSI < 0.05, :);                 % selective units
    nexttile(tl);
    nd = max([4; cellfun(@numel, a.dirs)]);                       % directions tested (4, or 8 in 45-degree sets)
    polarhistogram(a.prefDir, linspace(-pi/nd, 2*pi - pi/nd, nd + 1), 'FaceColor', [0.3 0.3 0.3]);   % one bin per tested direction
    title(sprintf('%s (n = %d)', stimLab(k), height(a)));         % panel name and selective units
    ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'counterclockwise';   % MB phi: 0 = up, pi/2 = left
    ax.ThetaTick = [0 90 180 270]; ax.ThetaTickLabel = {'up', 'left', 'down', 'right'};   % motion direction on the screen
end
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'PreferredDirection.jpg']), 'Resolution', 300);

%% Example direction-tuned MB unit: a priori selection and tuning curve
% Rule fixed before looking at any raster: MB units responsive AND present (A),
% direction-selective (pDSI < 0.05) and crossing their RF cell >= minCross times
% in every direction, ranked by DSIcorr (DSI above chance). The top unit is
% plotted unless exPick names another one. Tuning curve = the grid test's own
% per-trial values (ObsDiff: mean rate in the 300 ms after the ball crosses the
% RF cell minus the 300 ms pre-trial baseline), mean +/- SEM per direction.
% Only interior RF cells (interiorOnly = true, user 2026-11-03): an edge cell is
% crossed ~230 ms after onset in the entering directions and ~250 ms before offset
% in the others, so early- vs late-crossing directions also differ in time within
% the trial; an interior cell avoids that ambiguity in the example. (The top edge
% units, rec 75 phy 32 and 202, cell 81, were checked 2026-11-03: no onset response
% in late-crossing trials and left >> up at equal crossing times, i.e. directional.)
exPick   = [];                                                    % [exp phyID] overrides the rule; [] = top candidate
minCross = 5;                                                     % crossings of the RF cell required in every direction
interiorOnly = true;                                              % only RF cells off the grid edge (crossed >= 470 ms from onset and offset)
gs       = 9;                                                     % grid size of the MB test (StatisticsPerNeuron: GridSize = 9)
arrows   = char([8593 8598 8592 8601 8595 8600 8594 8599]);       % up, up-left, left, down-left, down, down-right, right, up-right
dNames   = ["up" "up-left" "left" "down-left" "down" "down-right" "right" "up-right"];   % same order (MB phi, 45-degree steps)
dIdx     = @(phi) mod(round(phi / (pi/4)), 8) + 1;                % MB phi (0 = up, pi/2 = left) -> nearest 45-degree label
M    = A(A.stim == "MB", :);                                      % analysable MB units
[rC, cC] = ind2sub([gs gs], M.cell);                              % grid row and column of each RF cell
inner = rC > 1 & rC < gs & cC > 1 & cC < gs;                      % RF cell off the grid edge
cand = M(M.pDSI < 0.05 & cellfun(@(n) all(n >= minCross), M.nPerDir) & (inner | ~interiorOnly), :);   % selective, enough crossings, interior
cand = sortrows(cand, 'DSIcorr', 'descend');                      % strongest selectivity above chance first
assert(height(cand) > 0, 'Fig5: no MB unit passes the example-unit criteria.');   % fail loudly
top  = cand(1:min(5, height(cand)), {'exp','phyID','animal','DSIcorr','pDSI','prefDir'});   % top 5 candidates
top.prefDeg  = rad2deg(top.prefDir);                              % preferred direction in degrees (MB phi)
top.prefName = dNames(dIdx(top.prefDir))';                        % nearest screen direction
top.peakRate = cellfun(@max, cand.R(1:height(top)));              % largest mean delta rate over directions (spikes/s): flags sparse units
fprintf('Example-unit candidates: %d of %d analysable MB units pass (interiorOnly = %d); top %d:\n', height(cand), height(M), interiorOnly, height(top));
disp(top);
if isempty(exPick)
    ex1 = cand(1, :);                                             % rule: top-ranked candidate
else
    ex1 = M(M.exp == exPick(1) & M.phyID == exPick(2), :);        % the chosen unit (key match)
    assert(height(ex1) == 1, 'Fig5: exPick [%d %d] is not one analysable MB unit.', exPick(1), exPick(2));   % exactly one row
    if ~any(cand.exp == exPick(1) & cand.phyID == exPick(2))
        warning('Fig5:exampleNotCandidate', 'Unit %d/%d does not pass the example-unit criteria.', exPick(1), exPick(2));   % flag a post-hoc pick
    end
end

vsEx = linearlyMovingBallAnalysis(loadNPclassFromTable(ex1.exp)); % 9x9 MB session, as in directionTuningUnits
StEx = vsEx.StatisticsPerNeuron;                                  % bare call: figure-design grid statistics (errors if stale)
spF  = sprintf('Speed%d', ex1.speed);                             % the unit's Figure 1 speed
g    = StEx.(spF);                                                % grid results at that speed
assert(StEx.gridParams.GridSize == gs, 'Fig5: grid is %dx%d, the interior rule assumes %dx%d.', ...
    StEx.gridParams.GridSize, StEx.gridParams.GridSize, gs, gs);  % edge rule matches the test's grid
win  = StEx.gridParams.GridAnalysisWindow;                        % test window after the crossing = baseline length (ms)
sDur = vsEx.ResponseWindow.(spF).stimDur;                         % stimulus duration at this speed (ms)
ps   = vsEx.dataObj.convertPhySorting2tIc(vsEx.spikeSortingFolder);   % spike sorting
isG  = string(ps.label') == 'good';                               % good units = columns of the statistics
gid  = ps.phy_ID(isG);                                            % their phy IDs
assert(gid(ex1.unitIdx) == ex1.phyID, 'Fig5: unit %d of rec %d is phy %d, the cache says phy %d.', ...
    ex1.unitIdx, ex1.exp, gid(ex1.unitIdx), ex1.phyID);           % unit ordering matches the cache
c    = ex1.cell;                                                  % RF grid cell of the tuning curve
x    = g.ObsDiff(:, c, ex1.unitIdx) * 1000;                       % per-trial response - baseline (spikes/s); NaN = cell not crossed
dr   = g.trialDirections(:);                                      % trial directions (MB phi, rad)
uD   = unique(dr(~isnan(x)));                                     % directions crossing the cell, ascending (up, left, down, right)
mR = nan(numel(uD), 1); sR = mR; nR = mR;                         % mean, SEM, crossings per direction
for k = 1:numel(uD)
    xk = x(dr == uD(k) & ~isnan(x));                              % this direction's crossings
    mR(k) = mean(xk); sR(k) = std(xk) / sqrt(numel(xk)); nR(k) = numel(xk);   % mean, SEM, n
end
[okD, ia] = ismembertol(uD, ex1.dirs{1}, 1e-4);                   % cached tuning-curve column of each direction
assert(all(okD) && max(abs(mR - ex1.R{1}(ia)'), [], 'all') < 1e-6 && isequal(nR, ex1.nPerDir{1}(ia)'), ...
    'Fig5: tuning curve of rec %d phy %d differs from the Fig5 cache (stale cache?).', ex1.exp, ex1.phyID);   % same data as DSI/pDSI
cr = g.crossingMs(~isnan(x), c);                                  % crossing times of the analysed trials (ms from onset)
if any(cr < win | cr > sDur - win)
    warning('Fig5:exampleEdgeTiming', ['Rec %d phy %d: RF-cell crossings at %.0f-%.0f ms of a %.0f ms trial; ', ...
        'its window holds onset or is cut at offset in some directions.'], ex1.exp, ex1.phyID, min(cr), max(cr), sDur);   % confound check
end
lab = string(arrows(dIdx(uD))') + " " + dNames(dIdx(uD))';        % tick labels: arrow + name
ttl = sprintf('Rec %d, phy %d (%s), cell %d: DSI %.2f, DSIcorr %.2f, p_{DSI} = %.4f, preferred %s (%.0f%s)', ...
    ex1.exp, ex1.phyID, ex1.animal, c, ex1.DSI, ex1.DSIcorr, ex1.pDSI, dNames(dIdx(ex1.prefDir)), rad2deg(ex1.prefDir), char(176));   % unit and its selectivity
fig = figure('Position', [100 100 820 380]);
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
pax = polaraxes(tl); pax.Layout.Tile = 1; hold(pax, 'on');        % polar tuning curve
rr  = max(mR, 0);                                                 % negatives clamped to 0, as in the DSI
polarplot(pax, [uD; uD(1)], [rr; rr(1)], 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');   % closed tuning curve
for k = 1:numel(uD)
    polarplot(pax, [uD(k) uD(k)], max(mR(k) + [-1 1] * sR(k), 0), 'k-', 'LineWidth', 1);     % +/- SEM
end
polarplot(pax, [ex1.prefDir ex1.prefDir], [0 max(rr)], 'r-', 'LineWidth', 2);              % vector-sum preferred direction
pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'counterclockwise'; % MB phi: 0 = up, pi/2 = left
pax.ThetaTick = rad2deg(uD'); pax.ThetaTickLabel = cellstr(lab);  % one tick per tested direction
title(pax, 'Rectified mean (spikes/s)');
ax = nexttile(tl, 2); hold(ax, 'on');                             % linear tuning curve
errorbar(ax, 1:numel(uD), mR, sR, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'CapSize', 6);   % mean +/- SEM (not rectified)
yline(ax, 0, 'k:');                                               % baseline level
set(ax, 'XTick', 1:numel(uD), 'XTickLabel', cellstr(lab + " (n=" + nR + ")")); xlim(ax, [0.5 numel(uD) + 0.5]);   % directions and crossings
xlabel(ax, 'Ball motion direction'); ylabel(ax, '\Delta rate, 300 ms after RF crossing (spikes/s)');
title(tl, ttl, 'FontSize', 9);
exportgraphics(fig, fullfile(d.paperFigs, sprintf('%sExampleUnit_rec%d_phy%d_tuning.jpg', paperFigPrefix(), ex1.exp, ex1.phyID)), 'Resolution', 300);

%% Example direction-tuned MB unit: raster per direction (desktop session: blank under matlab -batch)
% One block per direction, the trials crossing the RF cell (as the test, cf.
% gridAlignedMR): columns before 0 are the pre-trial baseline, columns from 0
% start at the ball crossing the cell. Shaded: the test's response window.
% Red ticks: stimulus offset of each trial (the window is cut there).
preB = win;                                                       % baseline length = test window (300 ms)
post = 500;                                                       % ms shown after the crossing (plotRaster_MultiExp gridPostMs)
spk  = ps.ic(:, isG); spk = spk(:, ex1.unitIdx);                  % spike-sorting column of the unit
fig = figure('Position', [100 100 460 640]); ax = axes(fig); hold(ax, 'on');
y0 = 0;                                                           % rows used by previous blocks
for k = 1:numel(uD)
    rows = dr == uD(k) & ~isnan(g.crossingMs(:, c));              % this direction's trials crossing the RF cell
    assert(sum(rows) == nR(k), 'Fig5: %d crossings in the raster, %d in the tuning curve.', sum(rows), nR(k));   % same trials
    ons  = g.trialOnsets(rows);                                   % trial onsets (ms)
    MR   = eventAlignedMR(spk, ps.t, ons, ons + g.crossingMs(rows, c), preB, post);   % [trials x (preB + post)] baseline | crossing
    [tr, tb] = find(MR > 0);                                      % trial and 1 ms bin of every spike
    plot(ax, tb - preB - 0.5, y0 + tr, 'k.', 'MarkerSize', 5);   % bin centre, 0 = crossing
    off = sDur - g.crossingMs(rows, c);                           % stimulus offset relative to the crossing (ms)
    vis = off < post;                                             % offsets inside the plotted window
    plot(ax, off(vis), y0 + find(vis), 'r|', 'MarkerSize', 4);    % offset tick per trial
    text(ax, -preB - 15, y0 + nR(k) / 2, arrows(dIdx(uD(k))), 'FontSize', 16, 'HorizontalAlignment', 'right');   % direction arrow
    y0 = y0 + nR(k);                                              % next block starts below
    yline(ax, y0 + 0.5, 'k-');                                    % block separator
end
hw = patch(ax, [0 win win 0], [0.5 0.5 y0 + 0.5 y0 + 0.5], [1 0.85 0.85], 'EdgeColor', 'none');   % response window
uistack(hw, 'bottom');                                            % shading behind the spikes
xline(ax, 0, 'k-', 'LineWidth', 1.5);                             % RF crossing (left of it: pre-trial baseline)
set(ax, 'YDir', 'reverse', 'YTick', [], 'Layer', 'top'); xlim(ax, [-preB post]); ylim(ax, [0.5 y0 + 0.5]);
xlabel(ax, 'Time from RF-cell crossing (ms); < 0: pre-trial baseline');   % rows: trials, one block per direction (arrows)
title(ax, sprintf('Rec %d, phy %d, cell %d', ex1.exp, ex1.phyID, c));
exportgraphics(fig, fullfile(d.paperFigs, sprintf('%sExampleUnit_rec%d_phy%d_raster.jpg', paperFigPrefix(), ex1.exp, ex1.phyID)), 'Resolution', 300);

%% Same tuning across stimuli
% Strict: preferred-direction agreement (mean cos of the difference) in units
% direction-selective to BOTH stimuli. Relaxed: correlation of the two tuning
% curves in units responsive and present for both. Both against a null that
% re-pairs units within the same insertion (keeps insertion-level structure).
% SDGm-MBR co-occur only in PV97: descriptive.
rng(1);
pairs = ["MB" "SDGm"; "MB" "MBR"; "SDGm" "MBR"];                  % compared stimulus pairs
ag = table();
for k = 1:size(pairs, 1)
    a1 = A(A.stim == pairs(k,1), {'exp','phyID','animal','R','dirs','prefDir','pDSI'});   % first stimulus
    a2 = A(A.stim == pairs(k,2), {'exp','phyID','R','dirs','prefDir','pDSI'});            % second stimulus
    a2.Properties.VariableNames(3:6) = {'R2','dirs2','prefDir2','pDSI2'};                 % explicit names (no join suffixes)
    J  = innerjoin(a1, a2, 'Keys', {'exp','phyID'});              % same unit, responsive + present to both (key pairing)
    S  = J(J.pDSI < 0.05 & J.pDSI2 < 0.05, :);                    % selective to both
    obsCos = NaN; pCos = NaN;                                     % strict test (needs >= 3 units)
    if height(S) >= 3
        obsCos = mean(cos(S.prefDir - S.prefDir2));               % 1 = same preferred direction
        nul = zeros(nBoot, 1);
        for b = 1:nBoot
            q = rePair(S.exp);                                    % partners re-paired within insertions
            nul(b) = mean(cos(S.prefDir - S.prefDir2(q)));
        end
        pCos = (sum(nul >= obsCos) + 1) / (nBoot + 1);            % one-tailed: more agreement than chance
    end
    obsR = NaN; pR = NaN;                                         % relaxed test
    if height(J) >= 3
        obsR = mean(arrayfun(@(i) tuneCorr(J.R{i}, J.dirs{i}, J.R2{i}, J.dirs2{i}), (1:height(J))'), 'omitnan');  % mean curve correlation
        nul = zeros(nBoot, 1);
        for b = 1:nBoot
            q = rePair(J.exp);                                    % partners re-paired within insertions
            nul(b) = mean(arrayfun(@(i) tuneCorr(J.R{i}, J.dirs{i}, J.R2{q(i)}, J.dirs2{q(i)}), (1:height(J))'), 'omitnan');
        end
        pR = (sum(nul >= obsR) + 1) / (nBoot + 1);                % one-tailed: more similar than chance
    end
    ag = [ag; table(pairs(k,1), pairs(k,2), height(J), numel(unique(J.animal)), obsR, pR, height(S), obsCos, pCos, ...
        'VariableNames', {'stimA','stimB','nBoth','nAnimals','meanCurveR','pCurveR','nSelectiveBoth','meanCosPref','pCosPref'})];   %#ok<AGROW>
end
disp(ag);
writetable(ag, fullfile(d.paperFigs, 'Fig5_same_tuning_across_stimuli.csv'));

%% Leave-one-animal-out: fraction direction-selective
lo = table();
for s = stims
    a = A(A.stim == s, :); anS = unique(a.animal);
    for k = 1:numel(anS)
        b0 = a(a.animal ~= anS(k), :);                            % drop one animal
        if numel(unique(b0.animal)) < 2, continue, end
        [g, ex, an] = findgroups(b0.exp, b0.animal);
        f = splitapply(@(p) mean(p < 0.05), b0.pDSI, g);
        b = hierBootMatchFreq(f, nBoot, double(categorical(an)), ex);
        lo = [lo; table(s, anS(k), mean(b), prctile(b, 2.5), prctile(b, 97.5), mean(b <= 0.05), ...
            'VariableNames', {'stim','dropped','fracDS','ciLo','ciHi','pVsChance'})];   %#ok<AGROW>
    end
end
disp(lo);
writetable(lo, fullfile(d.paperFigs, 'Fig5_LeaveOneAnimalOut_fractionDS.csv'));

%% Local functions
function q = rePair(ins)
% Random re-pairing of partners within each insertion (permutation null).
q = (1:numel(ins))';                                              % identity pairing
for ex = unique(ins)'
    ii = find(ins == ex);                                         % units of this insertion
    q(ii) = ii(randperm(numel(ii)));                              % shuffle partners among them
end
end

function c = tuneCorr(R1, d1, R2, d2)
% Pearson correlation of two tuning curves over their shared directions.
[ok, ia] = ismembertol(d1, d2, 1e-4);                             % directions present in both
x = R1(ok); y = R2(ia(ok)); c = NaN;                              % aligned responses
if numel(x) >= 3 && std(x) > 0 && std(y) > 0
    m = corrcoef(x, y); c = m(1, 2);                              % correlation coefficient
end
end
