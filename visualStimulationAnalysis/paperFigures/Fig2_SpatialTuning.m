%% FIGURE 2 — Spatial tuning of moving-ball (MB) vs static-ball (SB, rect grid) responses
% Units: responsive (permutation p < 0.05) AND present in both the MB and SB blocks, read
% from the Figure 1 statistics caches (run Fig1_MBvsRG.m first). MB uses, per unit,
% the speed with the lowest p, as in Figure 1.
% Index: L_amplitude_diff on 9x9 grid maps of the receptive field, split half
% (MB preferred direction chosen on odd repetitions and the index measured on
% even, and vice versa; SB averaged over the same halves).
% Stripes: anisotropy on the 9x9 stimulus grid vs a pixel-permutation null; MB
% takes the max over directions in the observed score and in every surrogate.
% Stripe baseline: the same pipeline on non-responsive units (p > 0.5, unitSet
% "nonResponsive"); % stripe in responsive vs baseline units is shown with CIs only.
% Group comparisons in analyzeStripeNeurons that select on the same RFs are
% descriptive (CIs only, no p-value).
% Run the cells in order; rasters must be drawn in an interactive session.

d    = figureDesign();                                            % shared design (recordings, output folder)
paperFigPrefix('fig2');                                           % saved panels are named fig2_<panel>
exps = intersect(d.exps, [40:43, 49:54, 64:66, 68:85, 87:97]);    % MB with 9 parallel trajectories (full 9x9 coverage):
                                                                  % excludes 1, 8 (single trajectory, 8 directions, SB 4x4);
                                                                  % 84, 88, 96 use their second (9x9) MB session
                                                                  % (linearlyMovingBallAnalysis picks it by default)
tuneArgs = {'indexType', "L_amplitude_diff", 'topPercent', 5, 'useRF', true, 'prefDir', true, ...
    'onOff', 1, 'sizeIdx', 1, 'lumValue', 255, 'unionResponsive', false, 'allResponsive', true, ...
    'splitHalf', true};                                           % index design (shared by both calls); white only, as Figure 1

%% Receptive fields on each half of the repetitions (slow: run once, cached per recording)
% Computes and caches ...-odd/-even RF files (and speed 2 where a unit's Figure-1
% speed is 2). SpatialTuningIndex would do this on demand; this cell only makes
% the long computation explicit and resumable.
for ex = exps
    NP = loadNPclassFromTable(ex);                                % recording
    for half = ["all", "odd", "even"]
        rectGridAnalysis(NP).CalculateReceptiveFields(AllResponsiveNeurons=false, trialHalf=half);   % SB maps
        mb = linearlyMovingBallAnalysis(NP);                      % MB analysis object
        assert(numel(unique(mb.VST.offsets)) == 9, ...            % Figure 2 needs the full 9x9 trajectory coverage
            'Fig2: recording %d has %d MB trajectories, not 9.', ex, numel(unique(mb.VST.offsets)));
        for sp = 1:numel(unique(mb.VST.speed))                    % every recorded speed
            mb.CalculateReceptiveFields(AllResponsiveNeurons=false, trialHalf=half, speed=sp);       % MB maps
        end
    end
end

%% Spatial tuning index and stripe detection (swarm MB vs SB, hierarchical bootstrap)
results = SpatialTuningIndex(exps, tuneArgs{:}, overwrite=true, PaperFig=true, ...
    plotRFs=false, plotRFunion=false, detectStripe=true, testStripeDetection=false);   % cache name does not encode every param: always overwrite

%% Paired locality score: units responsive AND present for both MB and SB, MB - SB per unit
pairedArgs = tuneArgs; pairedArgs{find(strcmp(pairedArgs(1:2:end), 'allResponsive')) * 2} = false;   % intersection mode
resPaired = SpatialTuningIndex(exps, pairedArgs{:}, overwrite=true, PaperFig=true, ...
    plotRFs=false, plotRFunion=false, detectStripe=true, testStripeDetection=false);   % paired hierarchical bootstrap of MB - SB

%% Paired MB - SB result: effect, 95% CI, p, n, animals (CSV)
% Reads the stats struct the paired call above saved (blocks{1,2}.diffs) and checks it
% against the paired cache, so a stale struct cannot be reported.
argVal = @(name) pairedArgs{find(strcmp(pairedArgs(1:2:end), name)) * 2};   % value of one name-value argument
pairName = sprintf('SpatialTuningIndex-%s-linearlyMovingBall-rectGrid-RF1-prefDir1-allResp0-tp-%d', ...
    argVal('indexType'), argVal('topPercent'));                   % statsName SpatialTuningIndex builds in paired mode
P = load(fullfile(d.paperFigs, [pairName '_statsStruct.mat']));    % saveStatsStruct output
assert(iscell(P.blocks) && numel(P.blocks) == 3 && isfield(P.blocks{1,2}, 'diffs'), ...
    'Fig2: %s has no paired-difference block.', pairName);       % expected {title, stats, unit} layout
D = P.blocks{1,2}.diffs;                                          % one row: MB - SB per unit
assert(height(D) == 1 && D.group == "MB-SB", 'Fig2: expected one MB-SB difference row.');   % paired mode only
if exist('resPaired', 'var') && isfield(resPaired, 'ps')          % paired call ran in this session
    assert(abs(resPaired.ps - D.p) < 1e-12, 'Fig2: stats struct p (%.4f) differs from resPaired.ps (%.4f).', D.p, resPaired.ps);   % same run
end
Tp = SpatialTuningIndex(exps, pairedArgs{:}, overwrite=false, plot=false, ...
    plotRFs=false, plotRFunion=false, detectStripe=true);         % paired cache (loads; a recompute matches the cell above)
Tp = Tp.tbl(Tp.tbl.onOff == 1 & Tp.tbl.sizeIdx == 1 & Tp.tbl.lumValue == 255, :);   % plotted condition (white)
J  = innerjoin(Tp(Tp.stimulus == "linearlyMovingBall", {'experimentNum','phyID','animal'}), ...
    Tp(Tp.stimulus == "rectGrid", {'experimentNum','phyID'}), 'Keys', {'experimentNum','phyID'});   % pair by key, not row
assert(height(J) == D.n, 'Fig2: paired cache has %d MB-SB units, stats struct has %d (stale struct).', height(J), D.n);   % same unit set
pairAnimals = unique(string(J.animal));                           % contributing animals
assert(numel(pairAnimals) == D.nAnimals, 'Fig2: animal count differs between cache and stats struct.');   % same animals
pairTbl = table(string(argVal('indexType')), D.n, D.nAnimals, D.nInsertions, strjoin(pairAnimals, ';'), ...
    D.mean, D.ciLo, D.ciHi, D.p, D.effect, D.effectLo, D.effectHi, D.propPos, D.propPosLo, D.propPosHi, string(P.savedOn), ...
    'VariableNames', {'index','nUnits','nAnimals','nInsertions','animals','meanDiff_MBminusSB','ciLo95','ciHi95', ...
    'pOneTailed_MBgtSB','cohenDz','cohenDzLo','cohenDzHi','propMBgtSB','propLo','propHi','statsSavedOn'});   % one row
disp(pairTbl);                                                    % p is one-tailed P(boot <= 0); the CI is two-sided 95%
writetable(pairTbl, fullfile(d.paperFigs, 'Fig2_paired_MBminusSB_locality.csv'));   % supplementary numbers

%% Stripe RFs: responsive units vs non-responsive baseline (descriptive, CIs only)
% Baseline = units present in the block with p > 0.5 for that stimulus (unitSet
% "nonResponsive"); same white-only maps, split half, per-unit speed, stripe settings
% and 1000 surrogates as the responsive set. Point = observed fraction of units with
% a stripe RF; 95% CI = hierarchical bootstrap of the 0/1 stripe flag (animal, insertion).
resResp = SpatialTuningIndex(exps, tuneArgs{:}, overwrite=false, plot=false, ...
    plotRFs=false, plotRFunion=false, detectStripe=true, testStripeDetection=false);   % responsive cache written two cells up
resBase = SpatialTuningIndex(exps, tuneArgs{:}, unitSet="nonResponsive", overwrite=true, plot=false, ...
    plotRFs=false, plotRFunion=false, detectStripe=true, testStripeDetection=false);   % baseline ('_nonResp' cache; slow)
stimLong  = ["linearlyMovingBall", "rectGrid"];                   % stimulus labels in the tables
stimShort = ["MB", "SB"];                                         % labels in the figure
setTbl    = {resResp.tbl, resBase.tbl};                           % responsive, baseline
setName   = ["responsive", "nonResponsive"];                      % unit-set labels
rng(0);                                                           % reproducible bootstrap
sfr = table();                                                    % one row per stimulus x unit set
for s = 1:2
    keysSet = cell(1, 2);                                         % (experiment, phyID) keys per set, for the disjointness check
    for k = 1:2
        t = setTbl{k};                                            % this unit set
        t = t(t.stimulus == stimLong(s) & t.onOff == 1 & t.sizeIdx == 1 & t.lumValue == 255, :);   % plotted condition (white)
        assert(height(unique(t(:, {'experimentNum','phyID'}))) == height(t), ...
            'Fig2: duplicate units in the %s %s table.', setName(k), stimShort(s));   % one row per unit
        assert(any(t.stripeBestPval < 1), 'Fig2: %s %s table has no stripe test results (cache without detectStripe).', ...
            setName(k), stimShort(s));                            % a detection run gives p < 1 for some unit
        keysSet{k} = t(:, {'experimentNum','phyID'});             % unit keys
        ins = str2double(string(t.experimentNum));                % insertion = recording number (globally unique)
        an  = findgroups(string(t.animal));                       % animal codes
        b   = hierBootMatchFreq(double(t.isStripe), 10000, an, ins);   % (data, nrep, ANIMAL, INSERTION)
        sfr = [sfr; table(stimShort(s), setName(k), height(t), sum(t.isStripe), numel(unique(an)), numel(unique(ins)), ...
            mean(t.isStripe), prctile(b, 2.5), prctile(b, 97.5), sum(t.isStripe & ins == 81), ...
            'VariableNames', {'stim','unitSet','nUnits','nStripe','nAnimals','nInsertions','fracStripe','ciLo95','ciHi95','nStripeRec81'})];   %#ok<AGROW>
    end
    nOverlap = height(innerjoin(keysSet{1}, keysSet{2}));         % units in both sets
    assert(nOverlap == 0, 'Fig2: %d %s units are in both the responsive and baseline sets.', nOverlap, stimShort(s));   % p < 0.05 vs p > 0.5
end
disp(sfr);
writetable(sfr, fullfile(d.paperFigs, 'Fig2_stripe_fraction_vs_baseline.csv'));   % supplementary numbers

fig = figure('Units', 'centimeters', 'Position', [5 5 7 6]); hold on;   % new panel
barCol = [0.25 0.25 0.25; 0.75 0.75 0.75];                        % responsive dark, baseline light
hBar   = gobjects(1, 2);                                          % bar handles for the legend
for k = 1:2
    r  = sfr(sfr.unitSet == setName(k), :);                       % this set, MB then SB
    x  = (1:2) + (k - 1.5) * 0.36;                                % grouped bar positions
    hBar(k) = bar(x, 100 * r.fracStripe, 0.34, 'FaceColor', barCol(k, :), 'EdgeColor', 'none');   % percent stripe
    errorbar(x, 100 * r.fracStripe, 100 * (r.fracStripe - r.ciLo95), 100 * (r.ciHi95 - r.fracStripe), ...
        'k.', 'LineWidth', 1, 'CapSize', 4);                      % 95% hierarchical-bootstrap CI
    for j = 1:2
        text(x(j), 100 * r.ciHi95(j), sprintf('%d/%d', r.nStripe(j), r.nUnits(j)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 6, 'FontName', 'helvetica');   % stripe units / units
    end
end
mbR = sfr(sfr.stim == "MB" & sfr.unitSet == "responsive", :);     % responsive MB row
set(gca, 'XTick', 1:2, 'XTickLabel', stimShort, 'XLim', [0.4 2.6], 'FontSize', 8, 'FontName', 'helvetica');
ylim([0, 100 * max(sfr.ciHi95) * 1.25 + eps]);                    % room for the n labels
ylabel('Units with stripe RF (% of units)');                      % percentage of units in each set
legend(hBar, {'Responsive', 'Non-responsive (p > 0.5)'}, 'Box', 'off', 'FontSize', 6, 'Location', 'northeast');   % bars only; top right stays clear of the MB labels
title(sprintf('MB stripes from rec 81: %d of %d', mbR.nStripeRec81, mbR.nStripe), 'FontSize', 7, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'StripeFraction_vs_baseline.jpg']), 'Resolution', 300);   % fig2_ panel

%% Stripe diagonal: original (BL->TR, 20-55 deg) vs mirror (BR->TL, 125-160 deg), descriptive
% Since 2026-11-12 both diagonals count as stripes (SpatialTuningIndex stripeDiagonal="both";
% the mirror was added after the BL->TR-only results had been seen). Fraction of units whose
% stripe follows each diagonal, per stimulus and unit set; 95% CI = hierarchical bootstrap of
% the 0/1 flag (animal, insertion), as the panel above. CIs only, no p-value.
diagName  = ["BLtoTR", "BRtoTL"];                                 % stripeDiagonal values
diagLabel = ["BL\rightarrowTR (20-55\circ)", "BR\rightarrowTL (125-160\circ)"];   % legend labels
rng(0);                                                           % reproducible bootstrap
sdg = table();                                                    % one row per stimulus x unit set x diagonal
for s = 1:2
    for k = 1:2
        t = setTbl{k};                                            % this unit set
        t = t(t.stimulus == stimLong(s) & t.onOff == 1 & t.sizeIdx == 1 & t.lumValue == 255, :);   % plotted condition (white)
        assert(ismember('stripeDiagonal', t.Properties.VariableNames), ...
            'Fig2: %s table has no stripeDiagonal column (cache from before 2026-11-12): rerun SpatialTuningIndex.', setName(k));
        assert(all(t.isStripe == (t.stripeDiagonal ~= "")), 'Fig2: isStripe and stripeDiagonal disagree (%s %s).', ...
            setName(k), stimShort(s));                            % every stripe has exactly one diagonal
        ins = str2double(string(t.experimentNum));                % insertion = recording number (globally unique)
        an  = findgroups(string(t.animal));                       % animal codes
        for g = 1:2
            flag = double(t.stripeDiagonal == diagName(g));       % 1 = stripe on this diagonal
            b    = hierBootMatchFreq(flag, 10000, an, ins);       % (data, nrep, ANIMAL, INSERTION)
            sdg  = [sdg; table(stimShort(s), setName(k), diagName(g), height(t), sum(flag), mean(flag), ...
                prctile(b, 2.5), prctile(b, 97.5), 'VariableNames', ...
                {'stim','unitSet','diagonal','nUnits','nStripe','fracStripe','ciLo95','ciHi95'})];   %#ok<AGROW>
        end
    end
end
disp(sdg);
writetable(sdg, fullfile(d.paperFigs, 'Fig2_stripe_diagonal_original_vs_mirror.csv'));   % supplementary numbers

fig = figure('Units', 'centimeters', 'Position', [5 5 11 6]);     % new panel: one tile per stimulus
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');                 % MB | SB
diagCol = [0.20 0.40 0.75; 0.85 0.45 0.10];                       % original blue, mirror orange
for s = 1:2
    ax = nexttile; hold(ax, 'on');                                % this stimulus
    hB = gobjects(1, 2);                                          % bar handles for the legend
    for g = 1:2
        r = sdg(sdg.stim == stimShort(s) & sdg.diagonal == diagName(g), :);   % responsive then non-responsive
        assert(height(r) == 2 && isequal(r.unitSet', setName), 'Fig2: unexpected rows for %s %s.', stimShort(s), diagName(g));
        x = (1:2) + (g - 1.5) * 0.36;                             % grouped bar positions (unit set on x)
        hB(g) = bar(ax, x, 100 * r.fracStripe, 0.34, 'FaceColor', diagCol(g, :), 'EdgeColor', 'none');   % percent of units
        errorbar(ax, x, 100 * r.fracStripe, 100 * (r.fracStripe - r.ciLo95), 100 * (r.ciHi95 - r.fracStripe), ...
            'k.', 'LineWidth', 1, 'CapSize', 4);                  % 95% hierarchical-bootstrap CI
        for j = 1:2
            text(ax, x(j), 100 * r.ciHi95(j), sprintf('%d/%d', r.nStripe(j), r.nUnits(j)), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 6, 'FontName', 'helvetica');   % stripe units / units
        end
    end
    set(ax, 'XTick', 1:2, 'XTickLabel', ["Responsive", "Non-resp."], 'XLim', [0.4 2.6], 'FontSize', 8, 'FontName', 'helvetica');
    ylim(ax, [0, 100 * max(sdg.ciHi95) * 1.3 + eps]);            % same y range in both tiles, room for labels
    title(ax, stimShort(s), 'FontSize', 8);                       % stimulus
    if s == 1, ylabel(ax, 'Units with stripe RF (% of units)'); end
    if s == 1, legend(ax, hB, diagLabel, 'Box', 'off', 'FontSize', 6, 'Location', 'northeast'); end   % MB tile: top right is free of labels
end
exportgraphics(fig, fullfile(d.paperFigs, [paperFigPrefix() 'StripeDiagonal_original_vs_mirror.jpg']), 'Resolution', 300);   % fig2_ panel

%% Stripe units vs responsive non-stripe units (depth, z, OSI/DSI, direction)
stripeArgs ={'indexType', "L_amplitude_diff", 'onOff', 1, 'sizeIdx', 1, 'lumValue', 255, 'useRF', true, ...
    'prefDir', true, 'unionResponsive', false, 'allResponsive', true, 'splitHalf', true};   % must match the call above
resStripe = analyzeStripeNeurons(exps, stripeArgs{:}, PaperFig=true);

%% Leave-one-animal-out: MB vs SB index (two-sample hierarchical bootstrap, as in SpatialTuningIndex)
rng(0);                                                           % reproducible bootstrap (unseeded until 2026-11-13: p drifted 0.330 -> 0.340 between runs)
S   = results;                                                    % SpatialTuningIndex output
tbl = S.tbl(S.tbl.onOff == 1 & S.tbl.sizeIdx == 1 & S.tbl.lumValue == 255, :);   % plotted condition (white)
animals = categories(removecats(tbl.animal));                     % animals in the figure
loao = table('Size', [numel(animals) + 1, 6], 'VariableTypes', {'string','double','double','double','double','double'}, ...
    'VariableNames', {'leftOut','nMB','nSB','meanMB','meanSB','pOneTailed'});   % one row per run
for k = 0:numel(animals)
    keep = true(height(tbl), 1);                                  % k = 0: all animals
    if k > 0, keep = tbl.animal ~= animals{k}; end                % drop one animal
    mb = tbl(keep & tbl.stimulus == "linearlyMovingBall", :);     % MB units
    sb = tbl(keep & tbl.stimulus == "rectGrid", :);               % SB units
    bMB = hierBootMatchFreq(mb.L_amplitude_diff, 10000, double(mb.animal), double(mb.experimentNum));   % animal, then insertion
    bSB = hierBootMatchFreq(sb.L_amplitude_diff, 10000, double(sb.animal), double(sb.experimentNum));
    name = "none"; if k > 0, name = string(animals{k}); end       % label of the left-out animal
    loao(k + 1, :) = {name, height(mb), height(sb), mean(mb.L_amplitude_diff), mean(sb.L_amplitude_diff), ...
        mean(bSB >= bMB)};                                        % P(SB >= MB), same convention as SpatialTuningIndex
end
disp(loao);
writetable(loao, fullfile(d.paperFigs, 'LeaveOneAnimalOut_Fig2_SpatialTuning.csv'));   % supplementary table

%% Example neurons: rasters and receptive fields (interactive session)
NP = loadNPclassFromTable(81);                                    % SB example recording
vsRe = rectGridAnalysis(NP);
vsRe.plotRaster(MergeNtrials=1, overwrite=true, AllResponsiveNeurons=false, exNeuronsPhyID=288, exNeurons=[], selectedLum=255, oneTrial=true, PaperFig=true)
colorbarLimsRG = vsRe.PlotReceptiveFields(exNeurons=15, allStimParamsCombined=false, PaperFig=true, overwrite=true);   % phy 288: position 15 with the 2026-11-05 labels (was 14)

% MB example by rule (decided 2026-11-07): among the 321 MB-responsive units of the unpaired index
% (white) with a tuning index >= the median (2.647), no stripe RF, and a grid-test RF cell off the
% 9x9 grid edge (as the Fig 5 example, Fig5_directionTuning_units.mat 'cell'), the unit with the
% highest MB z: rec 74 phy 292 (PV140, index 3.77, z 34.9, preferred direction 180 deg = down, RF
% right of centre). The legend must say "strongest unit by this rule". Replaces rec 74 phy 727
% (median rule, too sparse) and rec 75 phy 66 (RF in the screen corner).
NP = loadNPclassFromTable(74);                                    % MB example recording
vs = linearlyMovingBallAnalysis(NP);
vs.plotRaster('exNeuronsPhyID', 292, 'overwrite', true, 'MergeNtrials', 2, 'PaperFig', true, OneLuminosity='white', OneDirection='down', bin=20)
colorbarLims = vs.PlotReceptiveFields('exNeurons', 8, 'overwrite', true, 'OneDirection', 'down', 'OneLuminosity', 'white', 'PaperFig', true);   % phy 292: good-unit position 8 (checked 2026-11-07)
