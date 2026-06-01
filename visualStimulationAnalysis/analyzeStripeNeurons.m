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
% 4.  STRIPENESS SWARM (MB vs SB, paired diff)
% =========================================================================
% Union of stripe-classified neurons across stim types.  For each unique
% (experiment, phyID) in the union, look up the stripeBestScore in BOTH
% stim types' rows — a neuron classified MB-stripe still has an SB row
% (just with isStripe=false and a low score), so the paired comparison is
% well defined.
keysMB  = [double(string(mbStripe.experimentNum)), mbStripe.phyID];
keysSB  = [double(string(sbStripe.experimentNum)), sbStripe.phyID];
allKeys = unique([keysMB; keysSB], 'rows');
nUnion  = size(allKeys, 1);

mbScores   = nan(nUnion, 1);
sbScores   = nan(nUnion, 1);
animalsU   = strings(nUnion, 1);
insertions = nan(nUnion, 1);

mbExpVec = str2double(string(mbAll.experimentNum));
sbExpVec = str2double(string(sbAll.experimentNum));

for ri = 1:nUnion
    ex_ri = allKeys(ri, 1);
    ph_ri = allKeys(ri, 2);

    rMB = mbAll(mbExpVec == ex_ri & mbAll.phyID == ph_ri, :);
    rSB = sbAll(sbExpVec == ex_ri & sbAll.phyID == ph_ri, :);

    if ~isempty(rMB)
        mbScores(ri)   = rMB.stripeBestScore(1);
        animalsU(ri)   = string(rMB.animal(1));
        insertions(ri) = ex_ri;
    end
    if ~isempty(rSB)
        sbScores(ri) = rSB.stripeBestScore(1);
        if animalsU(ri) == ""
            animalsU(ri)   = string(rSB.animal(1));
            insertions(ri) = ex_ri;
        end
    end
end

% Long-format table for plotSwarmBootstrapWithComparisons
% (one row per (neuron, stim) combination)
nLong = sum(~isnan(mbScores)) + sum(~isnan(sbScores));
valueL   = nan(nLong, 1);
stimL    = strings(nLong, 1);
insL     = nan(nLong, 1);
animL    = strings(nLong, 1);
neurIDL  = nan(nLong, 1);
k = 1;
for ri = 1:nUnion
    if ~isnan(mbScores(ri))
        valueL(k)  = mbScores(ri);
        stimL(k)   = "MB";
        insL(k)    = insertions(ri);
        animL(k)   = animalsU(ri);
        neurIDL(k) = ri;
        k = k + 1;
    end
    if ~isnan(sbScores(ri))
        valueL(k)  = sbScores(ri);
        stimL(k)   = "SB";
        insL(k)    = insertions(ri);
        animL(k)   = animalsU(ri);
        neurIDL(k) = ri;
        k = k + 1;
    end
end

tblStripe = table(valueL, categorical(stimL), categorical(insL), ...
    categorical(animL), neurIDL, ...
    'VariableNames', {'value','stimulus','insertion','animal','NeurID'});

% Paired hierBoot on the per-neuron differences (only neurons with both)
bothMask = ~isnan(mbScores) & ~isnan(sbScores);
diffs    = mbScores(bothMask) - sbScores(bothMask);
insDiff  = insertions(bothMask);
[~, ~, animIdxDiff] = unique(animalsU(bothMask));

if numel(diffs) >= 3
    bootDiff = hierBoot(diffs, params.nBoot, insDiff, animIdxDiff);
    % two-tailed p
    pStripe  = 2 * min(mean(bootDiff <= 0), mean(bootDiff >= 0));
else
    pStripe = NaN;
end

fprintf('\n  Stripeness MB vs SB: paired n=%d, two-tailed p = %.4f\n', ...
    sum(bothMask), pStripe);

results.stripeData = tblStripe;
results.pStripe    = pStripe;
results.unionKeys  = allKeys;
results.mbScores   = mbScores;
results.sbScores   = sbScores;

if params.plot
    pairs = {'MB','SB'};
    figSwarm = plotSwarmBootstrapWithComparisons(tblStripe, pairs, pStripe, {'value'}, ...
        yLegend         = 'Stripeness (S = max/min var)', ...
        diff            = false, ...
        Alpha           = params.Alpha, ...
        plotMeanSem     = true, ...
        drawLines       = false, ...
        showBothAndDiff = true);

    ax = gca;
    ax.YAxis.FontSize = 8;  ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8;  ax.XAxis.FontName = 'helvetica';
    set(figSwarm, 'Units', 'centimeters', 'Position', [10 10 7 5]);

    swarmPdf = fullfile(saveDir, sprintf('StripeNeurons_swarm_%s.pdf', params.indexType));
    exportgraphics(figSwarm, swarmPdf, 'ContentType', 'vector');
    fprintf('  Saved: %s\n', swarmPdf);
    results.figSwarm = figSwarm;
end

% =========================================================================
% 5.  DEPTH PLOT (4 groups, custom swarm + bootstrap mean ± 95% CI)
% =========================================================================
groupOrder = {'MB stripe', 'SB stripe', 'MB resp', 'SB resp'};
groupRows  = {mbStripe,    sbStripe,    mbNonStr,  sbNonStr};
groupCols  = [0.85 0.20 0.20;   % MB stripe  — dark red
              0.20 0.40 0.85;   % SB stripe  — dark blue
              0.95 0.65 0.60;   % MB resp    — light red
              0.60 0.75 0.95];  % SB resp    — light blue

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

% Pairwise hierBoot (two-tailed) for the three comparisons of interest
results.pDepth.MBstripe_vs_MBresp   = compareGroupsHB(depthData, 'MB stripe', 'MB resp',   params.nBoot);
results.pDepth.SBstripe_vs_SBresp   = compareGroupsHB(depthData, 'SB stripe', 'SB resp',   params.nBoot);
results.pDepth.MBstripe_vs_SBstripe = compareGroupsHB(depthData, 'MB stripe', 'SB stripe', params.nBoot);

fprintf('\nDepth comparisons (two-tailed hierBoot p):\n');
fprintf('  MB stripe  vs  MB resp   : p = %.4f\n', results.pDepth.MBstripe_vs_MBresp);
fprintf('  SB stripe  vs  SB resp   : p = %.4f\n', results.pDepth.SBstripe_vs_SBresp);
fprintf('  MB stripe  vs  SB stripe : p = %.4f\n', results.pDepth.MBstripe_vs_SBstripe);

if params.plot
    figDepth = figure('Units','centimeters','Position',[5 5 11 8]);
    hold on;
    rng(7, 'twister');   % reproducible jitter

    nG = numel(groupOrder);
    for gi = 1:nG
        rows = depthData(depthData.group == groupOrder{gi}, :);
        if isempty(rows), continue; end

        % Jittered scatter
        x = gi + 0.25 * (rand(height(rows),1) - 0.5);
        scatter(x, rows.depth, 14, groupCols(gi,:), 'filled', ...
            'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');

        % Bootstrap mean and 95% CI
        try
            insNum  = double(rows.insertion);
            animNum = double(rows.animal);
            bM = hierBoot(rows.depth, params.nBoot, insNum, animNum);
            m  = mean(bM);
            ci = prctile(bM, [2.5 97.5]);
            errorbar(gi, m, m - ci(1), ci(2) - m, 'k', ...
                'LineWidth', 1.2, 'CapSize', 7);
            plot(gi, m, 'ko', 'MarkerFaceColor', 'w', ...
                'MarkerSize', 6, 'LineWidth', 1.2);
        catch ME
            warning('hierBoot failed for %s: %s', groupOrder{gi}, ME.message);
        end
    end

    set(gca, ...
        'XTick',       1:nG, ...
        'XTickLabel',  groupOrder, ...
        'YDir',        'reverse', ...     % deeper = lower on the page
        'XLim',        [0.5, nG+0.5], ...
        'FontSize',    8, ...
        'FontName',    'helvetica');
    ylabel('Cortical depth (\mum)', 'FontSize', 9);
    xtickangle(20);
    title('Depth: stripe vs responsive non-stripe', 'FontSize', 9);
    box on;
    hold off;

    depthPdf = fullfile(saveDir, sprintf('StripeNeurons_depth_%s.pdf', params.indexType));
    exportgraphics(figDepth, depthPdf, 'ContentType', 'vector');
    fprintf('  Saved: %s\n', depthPdf);
    results.figDepth = figDepth;
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
%  LOCAL FUNCTION
% =========================================================================
function pVal = compareGroupsHB(depthData, grpA, grpB, nBoot)
% Two-sample hierarchical bootstrap p-value (two-tailed) for the
% difference in means between two groups in depthData.

rowsA = depthData(depthData.group == grpA, :);
rowsB = depthData(depthData.group == grpB, :);

if height(rowsA) < 3 || height(rowsB) < 3
    pVal = NaN;
    return
end

try
    bA = hierBoot(rowsA.depth, nBoot, double(rowsA.insertion), double(rowsA.animal));
    bB = hierBoot(rowsB.depth, nBoot, double(rowsB.insertion), double(rowsB.animal));
    pOne = mean(bA >= bB);
    pVal = 2 * min(pOne, 1 - pOne);    % two-tailed
catch
    pVal = NaN;
end

end