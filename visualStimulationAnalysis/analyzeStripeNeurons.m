function results = analyzeStripeNeurons(exList, params)
% analyzeStripeNeurons  Compare stripe-classified neurons (MB vs SB) and
%                       look at their cortical depths.
%
%   Loads the cached SpatialTuningIndex table for the given exList /
%   parameter combination, joins it with cortical depths obtained via
%   getNeuronDepths, and produces two figures:
%
%     (1) Stripeness swarm plot:
%         MB stripeBestScore  vs  SB stripeBestScore  for the union of
%         neurons classified as stripe in EITHER stim type.  Paired
%         difference column shown; p-value via hierBoot on the per-neuron
%         differences (two-tailed).
%
%     (2) Depth swarm plot:
%         Four groups, side by side:
%             MB stripe   |  SB stripe   |  MB resp (non-stripe)  |  SB resp (non-stripe)
%         Bootstrap means and 95% CIs overlaid.  Pairwise two-tailed
%         hierBoot p-values reported for:
%             MB stripe   vs  MB resp
%             SB stripe   vs  SB resp
%             MB stripe   vs  SB stripe
%
%   PARAMETERS (must match the SpatialTuningIndex run that produced the cache)
%       indexType, onOff, sizeIdx, lumIdx        condition filter
%       useRF, prefDir, allResponsive, unionResponsive   filename-tag flags
%       nBoot      number of bootstrap resamples for hierBoot          (10000)
%       Alpha      transparency for swarm dots                          (0.4)
%       plot       generate and save figures                            (true)
%       saveDepthMap   cache the (exp,phyID) → depth lookup in tbl     (true)
%
%   OUTPUTS  (struct)
%       stripeData   long-format table fed to the swarm plot
%       depthData    long-format table for the depth plot
%       pStripe      paired-difference p-value (MB vs SB stripeness)
%       pDepth       struct of pairwise depth comparison p-values
%       figSwarm     handle to the stripeness swarm figure
%       figDepth     handle to the depth swarm figure
%       joinedTbl    full tbl with depth_um column appended (saved to disk)

arguments
    exList   double
    params.indexType        string  = "L_amplitude_diff"
    params.onOff            double  = 1
    params.sizeIdx          double  = 1
    params.lumIdx           double  = 1
    params.useRF            logical = true
    params.prefDir          logical = true
    params.allResponsive    logical = false
    params.unionResponsive  logical = false
    params.nBoot            double  = 10000
    params.Alpha            double  = 0.4
    params.plot             logical = true
    params.plotLegend       logical = false
    params.PaperFig         logical = false
    params.zStim            string  = "linearlyMovingBall"  % "linearlyMovingBall" | "rectGrid" | "both"
    params.idxStim          string  = "linearlyMovingBall"  % "linearlyMovingBall" | "rectGrid" | "both"
end

% =========================================================================
% 1.  Locate and load the cached SpatialTuningIndex results
% =========================================================================
NP_first = loadNPclassFromTable(exList(1));
vs_first = linearlyMovingBallAnalysis(NP_first);
pBase    = extractBefore(vs_first.getAnalysisFileName, 'lizards');
saveDir  = [pBase 'lizards\Combined_lizard_analysis'];

stimTypes    = ["linearlyMovingBall", "rectGrid"];
stimLabel    = strjoin(stimTypes, '-');
rfLabel      = ''; if params.useRF,           rfLabel      = '_RF';      end
prefLabel    = ''; if params.prefDir,         prefLabel    = '_prefDir'; end
allRespLabel = ''; if params.allResponsive,   allRespLabel = '_allResp'; end
unionLabel   = ''; if params.unionResponsive, unionLabel   = '_union';   end
nameOfFile   = sprintf('\\Ex_%d-%d_SpatialTuningIndex_%s%s%s%s%s.mat', ...
    exList(1), exList(end), stimLabel, rfLabel, prefLabel, allRespLabel, unionLabel);
fullPath = [saveDir nameOfFile];

if ~exist(fullPath, 'file')
    error(['SpatialTuningIndex cache not found:\n  %s\n' ...
           'Run SpatialTuningIndex with detectStripe=true and matching params first.'], ...
           fullPath);
end
S   = load(fullPath);
tbl = S.tbl;
fprintf('Loaded %s\n', fullPath);

if ~ismember('isStripe', tbl.Properties.VariableNames)
    error('Cached table lacks stripe columns. Re-run SpatialTuningIndex with detectStripe=true.');
end

% Apply the same condition filter that SpatialTuningIndex uses for plotting
idxCond = tbl.onOff == params.onOff & tbl.sizeIdx == params.sizeIdx & tbl.lumIdx == params.lumIdx;
tbl     = tbl(idxCond, :);

% =========================================================================
% 2.  Load (or compute) cortical depths and join into tbl
% =========================================================================
depthFile = fullfile(saveDir, 'NeuronDepths.mat');
if exist(depthFile, 'file')
    depthRes = load(depthFile);
    fprintf('Loaded cached depths from %s\n', depthFile);
else
    fprintf('Cached depths not found — running getNeuronDepths...\n');
    depthRes = getNeuronDepths(exList);
end

% Build (experiment, phyID) → depth lookup.
% depthRes.depthTable has Experiment + Unit (1:nGood index) + Depth_um.
% Unit index needs to be mapped to phyID via perExp.p_sort, since the
% SpatialTuningIndex table uses phyID as the neuron identifier.
depthMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ei = 1:numel(depthRes.perExp)
    pe = depthRes.perExp(ei);
    if isempty(pe.p_sort), continue; end

    lbl     = string(pe.p_sort.label');
    phy_IDg = pe.p_sort.phy_ID(lbl == 'good');

    expRows   = depthRes.depthTable.Experiment == pe.exNum;
    expDepths = depthRes.depthTable.Depth_um(expRows);
    expUnits  = depthRes.depthTable.Unit(expRows);

    for ui = 1:numel(expUnits)
        if expUnits(ui) <= numel(phy_IDg)
            depthMap(sprintf('%d_%d', pe.exNum, phy_IDg(expUnits(ui)))) = expDepths(ui);
        end
    end
end

% Append depth_um column
tbl.depth_um = nan(height(tbl), 1);
for ri = 1:height(tbl)
    ex_ri = str2double(string(tbl.experimentNum(ri)));
    key   = sprintf('%d_%d', ex_ri, tbl.phyID(ri));
    if depthMap.isKey(key)
        tbl.depth_um(ri) = depthMap(key);
    end
end
nDepthOK = sum(~isnan(tbl.depth_um));
fprintf('  Depth assigned to %d/%d rows (%d missing).\n', ...
    nDepthOK, height(tbl), height(tbl)-nDepthOK);

% =========================================================================
% 2b. Responsiveness z-score (pooled over all directions/conditions)
%     Recomputed per experiment×stimulus, then joined by phyID. Must run
%     before section 3 so the four group tables inherit the zResp column.
% =========================================================================
zp.nBoot          = params.nBoot;
zp.rngSeed        = 42;
zp.baseRespWindow = 1000;   % ms response window from onset (RG path)
zp.baselineBuffer = 200;    % ms gap before onset for baseline
zp.movingWindow   = 200;    % ms sliding window for MB peak

tbl.zResp = nan(height(tbl), 1);
zCache    = containers.Map('KeyType','char','ValueType','any');   % key 'exp_stim'

for ri = 1:height(tbl)
    ex_ri   = str2double(string(tbl.experimentNum(ri)));
    stim_ri = string(tbl.stimulus(ri));
    key     = sprintf('%d_%s', ex_ri, stim_ri);

    if ~zCache.isKey(key)
        try
            NP_ri = loadNPclassFromTable(ex_ri);
            switch stim_ri
                case "linearlyMovingBall", obj_ri = linearlyMovingBallAnalysis(NP_ri);
                case "rectGrid",           obj_ri = rectGridAnalysis(NP_ri);
                otherwise,                 obj_ri = [];
            end
            if isempty(obj_ri)
                zCache(key) = struct('z', [], 'phy', []);
            else
                [z_ri, phy_ri] = computeResponsivenessZ(obj_ri, zp);
                zCache(key) = struct('z', z_ri, 'phy', phy_ri);
            end
        catch ME
            warning('z-score failed for exp %d %s: %s', ex_ri, stim_ri, ME.message);
            zCache(key) = struct('z', [], 'phy', []);
        end
    end

    zc = zCache(key);
    if ~isempty(zc.phy)
        idx = find(zc.phy == tbl.phyID(ri), 1);
        if ~isempty(idx), tbl.zResp(ri) = zc.z(idx); end
    end
end
fprintf('  z-score assigned to %d/%d rows.\n', sum(~isnan(tbl.zResp)), height(tbl));

% =========================================================================
% 2c. Direction/orientation tuning (OSI, DSI) for MB neurons.
%     From DirectionTuning, joined by phyID. rectGrid has no direction
%     category, so SB rows stay NaN (OSI/DSI undefined there).
% =========================================================================
tbl.OSI = nan(height(tbl), 1);
tbl.DSI = nan(height(tbl), 1);

dtCache = containers.Map('KeyType','double','ValueType','any');   % key: experiment
for ri = 1:height(tbl)
    if tbl.stimulus(ri) ~= "linearlyMovingBall", continue; end    % MB only
    ex_ri = str2double(string(tbl.experimentNum(ri)));

    if ~dtCache.isKey(ex_ri)
        try
            NP_ri  = loadNPclassFromTable(ex_ri);
            obj_ri = linearlyMovingBallAnalysis(NP_ri);
            dt = DirectionTuning(obj_ri, 'save', true, 'overwrite', false);
            dtCache(ex_ri) = struct('OSI', dt.OSI, 'DSI', dt.DSI, 'phyID', dt.phyID);
        catch ME
            warning('DirectionTuning failed for exp %d: %s', ex_ri, ME.message);
            dtCache(ex_ri) = struct('OSI', [], 'DSI', [], 'phyID', []);
        end
    end

    dc = dtCache(ex_ri);
    if ~isempty(dc.phyID)
        idx = find(dc.phyID == tbl.phyID(ri), 1);
        if ~isempty(idx)
            tbl.OSI(ri) = dc.OSI(idx);
            tbl.DSI(ri) = dc.DSI(idx);
        end
    end
end
fprintf('  OSI/DSI assigned to %d MB rows.\n', sum(~isnan(tbl.OSI)));

% =========================================================================
% 3.  Split rows into the four groups of interest
% =========================================================================
mbAll    = tbl(tbl.stimulus == "linearlyMovingBall", :);
sbAll    = tbl(tbl.stimulus == "rectGrid",           :);
mbStripe = mbAll(mbAll.isStripe,  :);
sbStripe = sbAll(sbAll.isStripe,  :);
mbNonStr = mbAll(~mbAll.isStripe, :);
sbNonStr = sbAll(~sbAll.isStripe, :);

fprintf('\nStripe summary (condition onOff=%d, size=%d, lum=%d):\n', ...
    params.onOff, params.sizeIdx, params.lumIdx);
fprintf('  MB stripe: %4d / %4d responsive (%.1f%%)\n', ...
    height(mbStripe), height(mbAll), 100*height(mbStripe)/max(1,height(mbAll)));
fprintf('  SB stripe: %4d / %4d responsive (%.1f%%)\n', ...
    height(sbStripe), height(sbAll), 100*height(sbStripe)/max(1,height(sbAll)));

% =========================================================================
% 4.  STRIPENESS SWARM — MB stripe vs SB stripe (two independent groups)
%     Uses stripeBestScore of neurons that PASSED the stripe test only.
%     With allResponsive=true the two pools are independent, so we use a
%     two-sample hierBoot rather than a paired difference.
% =========================================================================

% Add value column pointing at stripeBestScore
mbStripe.value     = mbStripe.stripeBestScore;
mbStripe.insertion = mbStripe.experimentNum;
sbStripe.value     = sbStripe.stripeBestScore;
sbStripe.insertion = sbStripe.experimentNum;

% Rename stimulus labels to short forms
mbStripe.stimulus(mbStripe.stimulus == "linearlyMovingBall") = "MB";
sbStripe.stimulus(sbStripe.stimulus == "rectGrid")           = "SB";

% Long-format table — two independent groups, no NaN rows
tblStripe = [mbStripe(:, {'value','stimulus','insertion','animal'}); ...
             sbStripe(:, {'value','stimulus','insertion','animal'})];
tblStripe.NeurID = (1:height(tblStripe))';


% Two-sample hierBoot p-value
if height(mbStripe) >= 3 && height(sbStripe) >= 3
    bMB = hierBootMatchFreq(mbStripe.value, params.nBoot, ...
                   double(mbStripe.insertion), double(mbStripe.animal));
    bSB = hierBootMatchFreq(sbStripe.value, params.nBoot, ...
                   double(sbStripe.insertion), double(sbStripe.animal));
    pStripe = 2 * min(mean(bMB >= bSB), mean(bMB < bSB));   % two-tailed
else
    pStripe = NaN;
end

fprintf('\n  Stripeness MB (n=%d) vs SB (n=%d): two-sample two-tailed p = %.4f\n', ...
    height(mbStripe), height(sbStripe), pStripe);

results.stripeData = tblStripe;
results.pStripe    = pStripe;

% --- Section 4 plot: stripeness swarm ---
if params.plot
    yMaxVis = max(tblStripe.value, [], 'omitnan') * 1.15;  % ensure all dots visible
    pairs   = {'MB','SB'};
    figSwarm = plotSwarmBootstrapWithComparisons(tblStripe, pairs, pStripe, {'value'}, ...
        yLegend         = 'Stripeness (S = max/min var)', ...
        yMaxVis         = yMaxVis, ...
        diff            = false, ...
        showBothAndDiff = false, ...
        Alpha           = params.Alpha, ...
        plotMeanSem     = true, ...
        drawLines       = false);

    ax = gca;
    ax.YAxis.FontSize = 8;  ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8;  ax.XAxis.FontName = 'helvetica';
    set(figSwarm, 'Units', 'centimeters', 'Position', [10 10 5 5]);

   if params.PaperFig
        vs_first.printFig(figSwarm, sprintf('StripeNeurons-swarm-%s', params.indexType), ...
            PaperFig = params.PaperFig);
    end

    results.figSwarm = figSwarm;
end

% =========================================================================
% 5.  DEPTH PLOT (4 groups: stripe vs responsive non-stripe for MB and SB)
% =========================================================================
groupOrder = {'MB stripe', 'SB stripe', 'MB resp', 'SB resp'};
groupRows  = {mbStripe,    sbStripe,    mbNonStr,  sbNonStr};

% Build long-format depth table — one row per neuron per group
depthData = table([], categorical([]), categorical([]), categorical([]), [], ...
    'VariableNames', {'depth','group','insertion','animal','NeurID'});
for gi = 1:numel(groupOrder)
    rows = groupRows{gi};
    if isempty(rows), continue; end
    keep = ~isnan(rows.depth_um);
    if ~any(keep), continue; end
    addTbl = table(rows.depth_um(keep), ...
        categorical(repmat(groupOrder(gi), sum(keep), 1)), ...
        rows.experimentNum(keep), ...
        rows.animal(keep), ...
        rows.phyID(keep), ...
        'VariableNames', {'depth','group','insertion','animal','NeurID'});
    depthData = [depthData; addTbl]; %#ok<AGROW>
end

fprintf('\nDepth groups (rows with depth):\n');
for gi = 1:numel(groupOrder)
    fprintf('  %-10s n=%d\n', groupOrder{gi}, sum(depthData.group == groupOrder{gi}));
end

results.depthData = depthData;

% Pairwise hierBoot (two-tailed) — three comparisons of interest
results.pDepth.MBstripe_vs_MBresp   = compareGroupsHB(depthData, 'depth', 'MB stripe', 'MB resp',   params.nBoot);
results.pDepth.SBstripe_vs_SBresp   = compareGroupsHB(depthData, 'depth', 'SB stripe', 'SB resp',   params.nBoot);
results.pDepth.MBstripe_vs_SBstripe = compareGroupsHB(depthData, 'depth', 'MB stripe', 'SB stripe', params.nBoot);

fprintf('\nDepth comparisons (two-tailed hierBoot p):\n');
fprintf('  MB stripe  vs  MB resp   : p = %.4f\n', results.pDepth.MBstripe_vs_MBresp);
fprintf('  SB stripe  vs  SB resp   : p = %.4f\n', results.pDepth.SBstripe_vs_SBresp);
fprintf('  MB stripe  vs  SB stripe : p = %.4f\n', results.pDepth.MBstripe_vs_SBstripe);

% --- Section 5 plot: depth swarm with significance brackets ---
if params.plot
    figDepth = figure('Units','centimeters','Position',[5 5 14 8]);
    hold on;
    nG = numel(groupOrder);

    % Global animal → color mapping (consistent across all four groups)
    allAnimalNames = categories(removecats(depthData.animal));
    nAnimals       = numel(allAnimalNames);
    animalCmap     = lines(max(nAnimals, 1));

    for gi = 1:nG
        rows = depthData(depthData.group == groupOrder{gi}, :);
        if isempty(rows), continue; end

        [~, animalIdx] = ismember(string(rows.animal), allAnimalNames);
        dotColors = animalCmap(animalIdx, :);

        swarmchart(gi * ones(height(rows), 1), rows.depth, ...
            16, dotColors, 'filled', ...
            'XJitter',         'density', ...
            'XJitterWidth',    0.35, ...
            'MarkerFaceAlpha', 0.65, ...
            'MarkerEdgeColor', 'none');

        try
            bM = hierBootMatchFreq(rows.depth, params.nBoot, ...
                          double(rows.insertion), double(rows.animal));
            m  = mean(bM);
            ci = prctile(bM, [2.5 97.5]);
            errorbar(gi, m, m-ci(1), ci(2)-m, 'k', ...
                'LineWidth', 1.5, 'CapSize', 8);
            plot(gi, m, 'ko', 'MarkerFaceColor','w', ...
                'MarkerSize', 7, 'LineWidth', 1.5);
        catch ME
            warning('hierBoot failed for %s: %s', groupOrder{gi}, ME.message);
        end
    end

    % Significance brackets (reversed Y-axis: smaller depth = top of plot)
    allDepths   = depthData.depth(~isnan(depthData.depth));
    yMin        = min(allDepths);
    yMax        = max(allDepths);
    yRange      = yMax - yMin;
    bracketStep = yRange * 0.06;
    tickH       = yRange * 0.015;

    comps = {
        1, 2, results.pDepth.MBstripe_vs_SBstripe, 1;
        1, 3, results.pDepth.MBstripe_vs_MBresp,   2;
        2, 4, results.pDepth.SBstripe_vs_SBresp,   3;
        };
   
    % Significance brackets — only drawn when p < 0.05
    for ci = 1:size(comps, 1)
        g1 = comps{ci,1};  g2 = comps{ci,2};
        pV = comps{ci,3};  lv = comps{ci,4};

        if isnan(pV) || pV >= 0.05
            continue;   % skip non-significant comparisons — no clutter
        end

        yB = yMin - lv * bracketStep;
        plot([g1 g2], [yB yB],        'k-', 'LineWidth', 1);
        plot([g1 g1], [yB, yB+tickH], 'k-', 'LineWidth', 1);
        plot([g2 g2], [yB, yB+tickH], 'k-', 'LineWidth', 1);
        text((g1+g2)/2, yB+tickH*2.5, sigStars(pV), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize', 11, 'FontWeight', 'bold');
    end

    % Reserve vertical space only for the highest significant bracket level,
    % plus one extra step of padding so the label isn't clipped at the edge.
    maxSigLevel = 0;
    for ci = 1:size(comps, 1)
        if ~isnan(comps{ci,3}) && comps{ci,3} < 0.05
            maxSigLevel = max(maxSigLevel, comps{ci,4});
        end
    end
    % Extra 2 bracket-steps above the highest significant bracket so the
    % bracket line itself and any tick are never at the very edge.
    ylim([yMin - (maxSigLevel + 2) * bracketStep,  yMax + bracketStep]);


    % Animal legend
    lgdH = gobjects(nAnimals, 1);
    for ai = 1:nAnimals
        lgdH(ai) = plot(nan, nan, 'o', ...
            'Color', animalCmap(ai,:), 'MarkerFaceColor', animalCmap(ai,:), ...
            'MarkerSize', 6, 'DisplayName', allAnimalNames{ai});
    end
    
    if params.plotLegend
        legend(lgdH, 'Location','best', 'FontSize',7, 'Box','off');
    end

    set(gca, 'XTick',1:nG, 'XTickLabel',groupOrder, ...
        'YDir','reverse', 'XLim',[0.5 nG+0.5], ...
        'FontSize',8, 'FontName','helvetica');
    ylabel('Depth (\mum)', 'FontSize', 9);
    xtickangle(20);
    title('Depth: stripe vs responsive non-stripe', 'FontSize', 9);
    box on;  hold off;


    if params.PaperFig
        vs_first.printFig(figDepth, sprintf('StripeNeurons_depth_%s', params.indexType), ...
            PaperFig = params.PaperFig);
    end

    
    results.figDepth = figDepth;
end

% =========================================================================
% 7.  STRIPE DIRECTION DISTRIBUTION  (MB only) — flat histogram
%     stripeBestDir = index of the highest-SCORE passing direction.
%     Converted to radians, then matched to the 4 cardinal ball directions.
%     Convention: 0 = up, pi/2 = right, pi = down, 3*pi/2 = left.
% =========================================================================
if ~isempty(mbStripe) && ismember('stripeBestDir', mbStripe.Properties.VariableNames)

    dirCache     = containers.Map('KeyType','double','ValueType','any');
    stripeDirRad = nan(height(mbStripe), 1);

    for ri = 1:height(mbStripe)
        ex_ri  = str2double(string(mbStripe.experimentNum(ri)));
        dirIdx = mbStripe.stripeBestDir(ri);
        if isnan(dirIdx), continue; end

        if ~dirCache.isKey(ex_ri)
            try
                NP_ri    = loadNPclassFromTable(ex_ri);
                obj_ri   = linearlyMovingBallAnalysis(NP_ri);
                S_rf_ri  = obj_ri.CalculateReceptiveFields;
                rfSpeed  = S_rf_ri.params.speed;
                NR_ri    = obj_ri.ResponseWindow;
                NV_ri    = NR_ri.(sprintf('Speed%d', rfSpeed)).NeuronVals;
                uDirs_ri = unique(NV_ri(1,:,6));   % sorted unique directions (rad)
                dirCache(ex_ri) = uDirs_ri;
            catch ME
                warning('Could not load directions for exp %d: %s', ex_ri, ME.message);
                dirCache(ex_ri) = [];
            end
        end

        uDirsRad = dirCache(ex_ri);
        if ~isempty(uDirsRad) && dirIdx >= 1 && dirIdx <= numel(uDirsRad)
            stripeDirRad(ri) = uDirsRad(dirIdx);
        end
    end

    results.stripeDirRad = stripeDirRad;

    if params.plot
        validDir = stripeDirRad(~isnan(stripeDirRad));
        if ~isempty(validDir)
            % Reference directions and their cardinal labels
            dirRadRef = [0,    pi/2,    pi,    3*pi/2];
            dirNames  = {'up', 'right', 'down', 'left'};

            % Count neurons per direction (tolerance match, wrapped to [0,2pi))
            counts = zeros(1, numel(dirRadRef));
            for k = 1:numel(dirRadRef)
                counts(k) = sum(abs(mod(validDir, 2*pi) - dirRadRef(k)) < 1e-3);
            end

            figDir = figure('Units','centimeters','Position',[5 5 8 6]);
            bar(counts, 0.6, 'FaceColor', [0.85 0.2 0.2], 'EdgeColor', 'none');
            set(gca, 'XTick', 1:numel(dirNames), 'XTickLabel', dirNames, ...
                'FontSize', 8, 'FontName', 'helvetica');
            ylabel('Number of stripe neurons', 'FontSize', 9);
            xlabel('Ball direction', 'FontSize', 9);
            title(sprintf('MB stripe best direction (n=%d)', numel(validDir)), ...
                'FontSize', 9);
            box off;


            if params.PaperFig
                vs_first.printFig(figDir, sprintf('StripeNeurons_dirDist_%s', params.indexType), ...
                    PaperFig = params.PaperFig);
            end
 
            results.figDir = figDir;
        end
    end
end


% =========================================================================
% 8.  RESPONSIVENESS Z-SCORE: stripe vs responsive non-stripe
% =========================================================================
% Choose which stimulus/stimuli to show, since RG has few stripe neurons
switch params.zStim
    case "linearlyMovingBall"
        groupOrder = {'MB stripe', 'MB resp'};
        groupRows  = {mbStripe,    mbNonStr};
        zComps     = { 1, 2, [], 1 };          % MB stripe vs MB resp
    case "rectGrid"
        groupOrder = {'SB stripe', 'SB resp'};
        groupRows  = {sbStripe,    sbNonStr};
        zComps     = { 1, 2, [], 1 };          % SB stripe vs SB resp
    case "both"
        groupOrder = {'MB stripe', 'SB stripe', 'MB resp', 'SB resp'};
        groupRows  = {mbStripe,    sbStripe,    mbNonStr,  sbNonStr};
        zComps     = { 1, 3, [], 1;  2, 4, [], 2 };   % MB & SB stripe-vs-resp
end

zData = table([], categorical([]), categorical([]), categorical([]), [], ...
    'VariableNames', {'z','group','insertion','animal','NeurID'});
for gi = 1:numel(groupOrder)
    rows = groupRows{gi};
    if isempty(rows) || ~ismember('zResp', rows.Properties.VariableNames), continue; end
    keep = ~isnan(rows.zResp);
    if ~any(keep), continue; end
    zData = [zData; table(rows.zResp(keep), ...
        categorical(repmat(groupOrder(gi), sum(keep), 1)), ...
        rows.experimentNum(keep), rows.animal(keep), rows.phyID(keep), ...
        'VariableNames', {'z','group','insertion','animal','NeurID'})]; %#ok<AGROW>
end
results.zData = zData;

results.pZ = struct();
for ci = 1:size(zComps,1)
    gA = groupOrder{zComps{ci,1}};
    gB = groupOrder{zComps{ci,2}};
    pV = compareGroupsHB(zData, 'z', gA, gB, params.nBoot);
    zComps{ci,3} = pV;
    results.pZ.(matlab.lang.makeValidName([gA '_vs_' gB])) = pV;
    fprintf('  z: %s vs %s : p = %.4f\n', gA, gB, pV);
end

if params.plot && ~isempty(zData)
    figZ = figure('Units','centimeters','Position',[5 5 14 8]);
    hold on;
    nG = numel(groupOrder);
    allAnimalNames = categories(removecats(zData.animal));
    nAnimals       = numel(allAnimalNames);
    animalCmap     = lines(max(nAnimals,1));

    for gi = 1:nG
        rows = zData(zData.group == groupOrder{gi}, :);
        if isempty(rows), continue; end
        [~, aIdx] = ismember(string(rows.animal), allAnimalNames);
        swarmchart(gi*ones(height(rows),1), rows.z, 16, animalCmap(aIdx,:), 'filled', ...
            'XJitter','density', 'XJitterWidth',0.35, ...
            'MarkerFaceAlpha',0.65, 'MarkerEdgeColor','none');
        try
            bM = hierBootMatchFreq(rows.z, params.nBoot, double(rows.insertion), double(rows.animal));
            m = mean(bM);  ci = prctile(bM,[2.5 97.5]);
            errorbar(gi, m, m-ci(1), ci(2)-m, 'k', 'LineWidth',1.5, 'CapSize',8);
            plot(gi, m, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',1.5);
        catch ME
            warning('hierBoot failed for %s: %s', groupOrder{gi}, ME.message);
        end
    end

    % Significance brackets — normal Y-axis: brackets ABOVE the data
    yMax = max(zData.z);  yMin = min(zData.z);  zR = yMax - yMin;
    bStep = zR*0.10;  tickH = zR*0.025;
    % comps = {
    %     1, 2, results.pZ.MBstripe_vs_SBstripe, 1;
    %     1, 3, results.pZ.MBstripe_vs_MBresp,   2;
    %     2, 4, results.pZ.SBstripe_vs_SBresp,   3;
    % };
    maxSig = 0;
    for ci = 1:size(zComps,1)
        g1=zComps{ci,1}; g2=zComps{ci,2}; pV=zComps{ci,3}; lv=zComps{ci,4};
        if isnan(pV) || pV >= 0.05, continue; end
        yB = yMax + lv*bStep;                       % above data
        plot([g1 g2],[yB yB],        'k-','LineWidth',1);
        plot([g1 g1],[yB, yB-tickH], 'k-','LineWidth',1);   % ticks point DOWN to data
        plot([g2 g2],[yB, yB-tickH], 'k-','LineWidth',1);
        text((g1+g2)/2, yB, sigStars(pV), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
        maxSig = max(maxSig, lv);
    end
    ylim([yMin - bStep, yMax + (maxSig+1.5)*bStep]);

    if params.plotLegend
        lgdH = gobjects(nAnimals,1);
        for ai = 1:nAnimals
            lgdH(ai) = plot(nan,nan,'o','Color',animalCmap(ai,:), ...
                'MarkerFaceColor',animalCmap(ai,:),'MarkerSize',6, ...
                'DisplayName',allAnimalNames{ai});
        end
        legend(lgdH,'Location','best','FontSize',7,'Box','off');
    end

    set(gca,'XTick',1:nG,'XTickLabel',groupOrder, ...
        'XLim',[0.5 nG+0.5],'FontSize',8,'FontName','helvetica');
    ylabel('Responsiveness (z-score)','FontSize',9);
    xtickangle(20);
    title('Responsiveness: stripe vs non-stripe','FontSize',9);
    box on; hold off;

    if params.PaperFig
        vs_first.printFig(figZ, sprintf('StripeNeurons-zscore-%s', params.indexType), ...
            PaperFig = params.PaperFig);
    end
    results.figZ = figZ;
end

% =========================================================================
% 9.  TUNING INDEX: MB stripe vs MB responsive non-stripe (OSI and DSI)
%     MB-only — rectGrid has no direction tuning.
%     Columns: [OSI stripe | OSI resp | DSI stripe | DSI resp]
% =========================================================================
tuneOrder  = {'OSI stripe', 'OSI resp', 'DSI stripe', 'DSI resp'};
tuneSrc    = {mbStripe,     mbNonStr,   mbStripe,     mbNonStr};
tuneMetric = {'OSI',        'OSI',      'DSI',        'DSI'};

tuneData = table([], categorical([]), categorical([]), categorical([]), [], ...
    'VariableNames', {'val','group','insertion','animal','NeurID'});
for gi = 1:numel(tuneOrder)
    rows = tuneSrc{gi};
    mcol = tuneMetric{gi};
    if isempty(rows) || ~ismember(mcol, rows.Properties.VariableNames), continue; end
    vals = rows.(mcol);
    keep = ~isnan(vals);
    if ~any(keep), continue; end
    tuneData = [tuneData; table(vals(keep), ...
        categorical(repmat(tuneOrder(gi), sum(keep), 1)), ...
        rows.experimentNum(keep), rows.animal(keep), rows.phyID(keep), ...
        'VariableNames', {'val','group','insertion','animal','NeurID'})]; %#ok<AGROW>
end
results.tuneData = tuneData;

results.pTune.OSI = compareGroupsHB(tuneData, 'val', 'OSI stripe', 'OSI resp', params.nBoot);
results.pTune.DSI = compareGroupsHB(tuneData, 'val', 'DSI stripe', 'DSI resp', params.nBoot);

fprintf('\nTuning index (two-tailed hierBoot p):\n');
fprintf('  OSI stripe vs resp : p = %.4f\n', results.pTune.OSI);
fprintf('  DSI stripe vs resp : p = %.4f\n', results.pTune.DSI);

if params.plot && ~isempty(tuneData)
    figTune = figure('Units','centimeters','Position',[5 5 12 8]);
    hold on;
    nG = numel(tuneOrder);
    allAnimalNames = categories(removecats(tuneData.animal));
    nAnimals       = numel(allAnimalNames);
    animalCmap     = lines(max(nAnimals,1));

    for gi = 1:nG
        rows = tuneData(tuneData.group == tuneOrder{gi}, :);
        if isempty(rows), continue; end
        [~, aIdx] = ismember(string(rows.animal), allAnimalNames);
        swarmchart(gi*ones(height(rows),1), rows.val, 16, animalCmap(aIdx,:), 'filled', ...
            'XJitter','density', 'XJitterWidth',0.35, ...
            'MarkerFaceAlpha',0.65, 'MarkerEdgeColor','none');
        try
            bM = hierBootMatchFreq(rows.val, params.nBoot, double(rows.insertion), double(rows.animal));
            m = mean(bM);  ci = prctile(bM,[2.5 97.5]);
            errorbar(gi, m, m-ci(1), ci(2)-m, 'k', 'LineWidth',1.5, 'CapSize',8);
            plot(gi, m, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',1.5);
        catch ME
            warning('hierBoot failed for %s: %s', tuneOrder{gi}, ME.message);
        end
    end

    % Brackets: 1 vs 2 (OSI), 3 vs 4 (DSI) — disjoint spans, both at level 1
    yMax = max(tuneData.val);  yMin = min([0; tuneData.val]);  vR = yMax - yMin;
    bStep = vR*0.10;  tickH = vR*0.025;
    comps = { 1,2,results.pTune.OSI,1;  3,4,results.pTune.DSI,1 };
    maxSig = 0;
    for ci = 1:size(comps,1)
        g1=comps{ci,1}; g2=comps{ci,2}; pV=comps{ci,3}; lv=comps{ci,4};
        if isnan(pV) || pV >= 0.05, continue; end
        yB = yMax + lv*bStep;
        plot([g1 g2],[yB yB],        'k-','LineWidth',1);
        plot([g1 g1],[yB, yB-tickH], 'k-','LineWidth',1);
        plot([g2 g2],[yB, yB-tickH], 'k-','LineWidth',1);
        text((g1+g2)/2, yB, sigStars(pV), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
        maxSig = max(maxSig, lv);
    end
    ylim([yMin, yMax + (maxSig+1.5)*bStep]);

    if params.plotLegend
        lgdH = gobjects(nAnimals,1);
        for ai = 1:nAnimals
            lgdH(ai) = plot(nan,nan,'o','Color',animalCmap(ai,:), ...
                'MarkerFaceColor',animalCmap(ai,:),'MarkerSize',6, ...
                'DisplayName',allAnimalNames{ai});
        end
        legend(lgdH,'Location','best','FontSize',7,'Box','off');
    end

    set(gca,'XTick',1:nG,'XTickLabel',tuneOrder, ...
        'XLim',[0.5 nG+0.5],'FontSize',8,'FontName','helvetica');
    ylabel('Selectivity index','FontSize',9);
    xtickangle(20);
    title('Direction/orientation tuning: MB stripe vs non-stripe','FontSize',9);
    box on; hold off;

    if params.PaperFig
        vs_first.printFig(figTune, sprintf('StripeNeurons-tuning-%s', params.indexType), ...
            PaperFig = params.PaperFig);
    end
    results.figTune = figTune;
end

% =========================================================================
% 10.  SPATIAL TUNING INDEX: stripe vs responsive non-stripe
%      Uses params.indexType at each neuron's preferred direction (same for
%      both groups — symmetric). NOT stripeBestDir, to keep the comparison
%      apples-to-apples (non-stripe neurons have no stripe direction).
% =========================================================================
switch params.idxStim
    case "linearlyMovingBall"
        idxOrder = {'MB stripe', 'MB resp'};
        idxRows  = {mbStripe,    mbNonStr};
        idxComps = { 1, 2, [], 1 };
    case "rectGrid"
        idxOrder = {'SB stripe', 'SB resp'};
        idxRows  = {sbStripe,    sbNonStr};
        idxComps = { 1, 2, [], 1 };
    case "both"
        idxOrder = {'MB stripe', 'SB stripe', 'MB resp', 'SB resp'};
        idxRows  = {mbStripe,    sbStripe,    mbNonStr,  sbNonStr};
        idxComps = { 1, 3, [], 1;  2, 4, [], 2 };
end

idxData = table([], categorical([]), categorical([]), categorical([]), [], ...
    'VariableNames', {'val','group','insertion','animal','NeurID'});
for gi = 1:numel(idxOrder)
    rows = idxRows{gi};
    if isempty(rows), continue; end
    vals = rows.(params.indexType);          % the tuning index, already at prefDir
    keep = ~isnan(vals);
    if ~any(keep), continue; end
    idxData = [idxData; table(vals(keep), ...
        categorical(repmat(idxOrder(gi), sum(keep), 1)), ...
        rows.experimentNum(keep), rows.animal(keep), rows.phyID(keep), ...
        'VariableNames', {'val','group','insertion','animal','NeurID'})]; %#ok<AGROW>
end
results.idxData = idxData;

results.pIdx = struct();
for ci = 1:size(idxComps,1)
    gA = idxOrder{idxComps{ci,1}};
    gB = idxOrder{idxComps{ci,2}};
    pV = compareGroupsHB(idxData, 'val', gA, gB, params.nBoot);
    idxComps{ci,3} = pV;
    results.pIdx.(matlab.lang.makeValidName([gA '_vs_' gB])) = pV;
    fprintf('  idx: %s vs %s : p = %.4f\n', gA, gB, pV);
end

if params.plot && ~isempty(idxData) && any(~isnan(idxData.val))
    figIdx = figure('Units','centimeters','Position',[5 5 12 8]);
    hold on;
    nG = numel(idxOrder);
    allAnimalNames = categories(removecats(idxData.animal));
    nAnimals       = numel(allAnimalNames);
    animalCmap     = lines(max(nAnimals,1));

    for gi = 1:nG
        rows = idxData(idxData.group == idxOrder{gi}, :);
        if isempty(rows), continue; end
        [~, aIdx] = ismember(string(rows.animal), allAnimalNames);
        swarmchart(gi*ones(height(rows),1), rows.val, 16, animalCmap(aIdx,:), 'filled', ...
            'XJitter','density', 'XJitterWidth',0.35, ...
            'MarkerFaceAlpha',0.65, 'MarkerEdgeColor','none');
        vv = rows.val(~isnan(rows.val));
        if numel(vv) < 3, continue; end
        try
            bM = hierBootMatchFreq(vv, params.nBoot, ...
                double(rows.insertion(~isnan(rows.val))), ...
                double(rows.animal(~isnan(rows.val))));
            m = mean(bM);  ci2 = prctile(bM,[2.5 97.5]);
            errorbar(gi, m, m-ci2(1), ci2(2)-m, 'k', 'LineWidth',1.5, 'CapSize',8);
            plot(gi, m, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',1.5);
        catch ME
            warning('hierBoot failed for %s: %s', idxOrder{gi}, ME.message);
        end
    end

    allVals = idxData.val(~isnan(idxData.val));
    yMax = max(allVals);  yMin = min([0; allVals]);  vR = yMax - yMin;
    bStep = vR*0.10;  tickH = vR*0.025;  maxSig = 0;
    for ci = 1:size(idxComps,1)
        g1 = idxComps{ci,1};  g2 = idxComps{ci,2};
        pV = idxComps{ci,3};  lv = idxComps{ci,4};
        if isnan(pV) || pV >= 0.05, continue; end
        yB = yMax + lv*bStep;
        plot([g1 g2],[yB yB],        'k-','LineWidth',1);
        plot([g1 g1],[yB, yB-tickH], 'k-','LineWidth',1);
        plot([g2 g2],[yB, yB-tickH], 'k-','LineWidth',1);
        text((g1+g2)/2, yB, sigStars(pV), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',11,'FontWeight','bold');
        maxSig = max(maxSig, lv);
    end
    ylim([yMin, yMax + (maxSig+1.5)*bStep]);

    if params.plotLegend
        lgdH = gobjects(nAnimals,1);
        for ai = 1:nAnimals
            lgdH(ai) = plot(nan,nan,'o','Color',animalCmap(ai,:), ...
                'MarkerFaceColor',animalCmap(ai,:),'MarkerSize',6, ...
                'DisplayName',allAnimalNames{ai});
        end
        legend(lgdH,'Location','best','FontSize',7,'Box','off');
    end

    set(gca,'XTick',1:nG,'XTickLabel',idxOrder, ...
        'XLim',[0.5 nG+0.5],'FontSize',8,'FontName','helvetica');
    ylabel(params.indexType, 'FontSize',9, 'Interpreter','none');
    xtickangle(20);
    title('Spatial tuning index: stripe vs non-stripe','FontSize',9);
    box on; hold off;

    if params.PaperFig
        vs_first.printFig(figIdx, sprintf('StripeNeurons-tuningIndex-%s', params.indexType), ...
            PaperFig = params.PaperFig);
    end
    results.figIdx = figIdx;
end

% =========================================================================
% 6.  Save the joined table to disk for downstream use
% =========================================================================
results.joinedTbl = tbl;
joinedPath = fullfile(saveDir, sprintf('StripeNeurons_joinedTbl_%s.mat', params.indexType));
save(joinedPath, '-struct', 'results');
fprintf('\nJoined table + results saved to:\n  %s\n', joinedPath);

end


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================
function pVal = compareGroupsHB(dataTbl, valCol, grpA, grpB, nBoot)
rowsA = dataTbl(dataTbl.group == grpA, :);
rowsB = dataTbl(dataTbl.group == grpB, :);
rowsA = rowsA(~isnan(rowsA.(valCol)), :);     % drop NaN
rowsB = rowsB(~isnan(rowsB.(valCol)), :);
if height(rowsA) < 3 || height(rowsB) < 3
    pVal = NaN;
    return
end
try
    bA = hierBootMatchFreq(rowsA.(valCol), nBoot, double(rowsA.insertion), double(rowsA.animal));
    bB = hierBootMatchFreq(rowsB.(valCol), nBoot, double(rowsB.insertion), double(rowsB.animal));
    pOne = mean(bA >= bB);
    pVal = 2 * min(pOne, 1 - pOne);
catch
    pVal = NaN;
end
end

function s = sigStars(p)
% Convert a p-value to a significance star string for bracket annotation.
if     isnan(p) || p >= 0.05,  s = 'ns';
elseif p < 0.001,               s = '***';
elseif p < 0.01,                s = '**';
else,                           s = '*';
end
end

function [zScore, phyIDg] = computeResponsivenessZ(vsObj, zp)
% Pooled responsiveness z-score per good unit (all trials, all directions).
% Adapted from the sign-flip permutation z in StatisticsPerNeuronPerCategory:
%   z = (mean(response-baseline) - nullMean) / sd(baseline)
% computed identically for every neuron so stripe/non-stripe are comparable.

    p      = vsObj.dataObj.convertPhySorting2tIc(vsObj.spikeSortingFolder);
    label  = string(p.label');
    goodU  = p.ic(:, label == 'good');
    phyIDg = p.phy_ID(label == 'good');

    rw = vsObj.ResponseWindow;

    % Baseline window: ends baselineBuffer ms before onset, <= baseRespWindow long
    preStim   = vsObj.VST.interTrialDelay*1000 - zp.baselineBuffer;
    baseDur   = min(preStim, zp.baseRespWindow);
    baseStart = baseDur + zp.baselineBuffer;
    if baseDur <= 0
        error('Non-positive baseline duration (interTrialDelay too short).');
    end

    isMB = isequal(vsObj.stimName, 'linearlyMovingBall');

    if isMB
        nSpeeds  = numel(unique(vsObj.VST.speed));
        respList = cell(nSpeeds,1);
        baseList = cell(nSpeeds,1);
        for s = 1:nSpeeds
            fName = sprintf('Speed%d', s);
            tt = rw.(fName).C(:,1)';
            sd = rw.(fName).stimDur;
            Mr = BuildBurstMatrix(goodU, round(p.t), round(tt), round(sd));
            mrMov       = movmean(Mr, zp.movingWindow, 3, 'Endpoints','discard');
            respList{s} = max(mrMov, [], 3);                       % peak window per trial
            Mb = BuildBurstMatrix(goodU, round(p.t), round(tt - baseStart), baseDur);
            baseList{s} = mean(Mb, 3);
        end
        responsesFull = vertcat(respList{:});
        baselinesFull = vertcat(baseList{:});
    else
        tt = rw.C(:,1)';
        Mr = BuildBurstMatrix(goodU, round(p.t), round(tt), zp.baseRespWindow);
        responsesFull = mean(Mr, 3);
        Mb = BuildBurstMatrix(goodU, round(p.t), round(tt - baseStart), baseDur);
        baselinesFull = mean(Mb, 3);
    end

    DiffFull = responsesFull - baselinesFull;     % [nTrials x nNeurons]
    sdBase   = std(baselinesFull, 0, 1);          % [1 x nNeurons]
    nTrials  = size(DiffFull, 1);

    ObsStat = mean(DiffFull, 1);                  % [1 x nNeurons]

    rng(zp.rngSeed, 'twister');                   % reproducible sign-flip null
    signs    = 2*randi(2, nTrials, zp.nBoot) - 3; % {-1,+1}
    nullDist = (signs' * DiffFull) / nTrials;     % [nBoot x nNeurons]
    nullMean = mean(nullDist, 1);

    zScore = (ObsStat - nullMean) ./ sdBase;
    zScore(sdBase == 0) = 0;
    zScore = zScore(:);                            % column, aligned with phyIDg
end