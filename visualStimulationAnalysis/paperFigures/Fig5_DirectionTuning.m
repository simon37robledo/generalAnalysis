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
ctor   = {@linearlyMovingBallAnalysis, @StaticDriftingGratingAnalysis, @linearlyMovingBarAnalysis};  % analysis classes
cacheF = fullfile(d.paperFigs, 'Fig5_directionTuning_units.mat'); % per-unit table (cache)
nBoot  = 10000;                                                   % bootstrap / permutation size

%% Per-unit tuning for every recording (slow: cached)
U = table();                                                      % one row per unit x stimulus
for ex = d.exps
    NP = loadNPclassFromTable(ex);                                % recording
    animal = string(regexp(NP.recordingName, '(PV|SA)\d+', 'match', 'once'));   % animal ID
    for s = 1:numel(stims)
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
set(gca, 'XTick', 1:numel(stims), 'XTickLabel', stims); ylabel('Fraction direction-selective');
title('Direction-selective units (p < 0.05)');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'DirectionSelective_fraction.jpg']), 'Resolution', 300);
writetable(fr, fullfile(d.paperFigs, 'Fig5_fraction_direction_selective.csv'));

%% DSI above chance per unit (animal-colored), bootstrap mean and 95% CI
% DSI minus the unit's shuffled-label mean DSI (DSIcorr): raw DSI is biased upward
% with few trials (an untuned unit scores ~0.3), so chance sits at 0 here.
% Units selective by the permutation test (pDSI < 0.05) are outlined in black.
fig = figure('Position', [100 100 480 380]); hold on;
anAll = unique(A.animal); cmap = lines(numel(anAll));             % one color per animal
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
set(gca, 'XTick', 1:numel(stims), 'XTickLabel', stims); ylabel('DSI - chance DSI');
title('Direction selectivity above chance (black outline: p < 0.05)');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'DSI_swarm.jpg']), 'Resolution', 300);

%% Preferred directions of direction-selective units (screen directions)
fig = figure('Position', [100 100 900 300]);
tl = tiledlayout(1, numel(stims));
for k = 1:numel(stims)
    a = A(A.stim == stims(k) & A.pDSI < 0.05, :);                 % selective units
    nexttile(tl);
    polarhistogram(a.prefDir, linspace(-pi/8, 2*pi - pi/8, 9), 'FaceColor', [0.3 0.3 0.3]);   % 8 bins centred on 0, 45, ...
    title(sprintf('%s (n = %d)', stims(k), height(a)));
    ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'counterclockwise';   % MB phi: 0 = up, pi/2 = left
    ax.ThetaTick = [0 90 180 270]; ax.ThetaTickLabel = {'up', 'left', 'down', 'right'};   % motion direction on the screen
end
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'PreferredDirection.jpg']), 'Resolution', 300);

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
