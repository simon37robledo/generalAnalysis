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
        dotColors = animalCmap(animalIdx, :);   % [nDots × 3]

        swarmchart(gi * ones(height(rows), 1), rows.depth, ...
            16, dotColors, 'filled', ...
            'XJitter',         'density', ...
            'XJitterWidth',    0.35, ...
            'MarkerFaceAlpha', 0.65, ...
            'MarkerEdgeColor', 'none');

        try
            bM = hierBoot(rows.depth, params.nBoot, ...
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

    % ---- Significance brackets ----------------------------------------
    % YDir='reverse': small depth values are at the TOP of the plot.
    % Brackets are placed above the shallowest neuron → y < yMin.
    % Ticks point downward (toward data) → y increases from bracket toward data.
    allDepths   = depthData.depth(~isnan(depthData.depth));
    yMin        = min(allDepths);
    yMax        = max(allDepths);
    yRange      = yMax - yMin;
    bracketStep = yRange * 0.06;    % vertical gap between bracket levels
    tickH       = yRange * 0.015;  % tick length

    % Each row: [g1, g2, pVal, bracketLevel]
    % Level 1 is closest to the data, level 3 is furthest above
    comps = {
        1, 2, results.pDepth.MBstripe_vs_SBstripe, 1;  % MB stripe vs SB stripe
        1, 3, results.pDepth.MBstripe_vs_MBresp,   2;  % MB stripe vs MB resp
        2, 4, results.pDepth.SBstripe_vs_SBresp,   3;  % SB stripe vs SB resp
    };

    for ci = 1:size(comps, 1)
        g1   = comps{ci,1};
        g2   = comps{ci,2};
        pVal = comps{ci,3};
        lv   = comps{ci,4};

        yB = yMin - lv * bracketStep;   % bracket line (above data in reversed axis)

        plot([g1 g2], [yB yB], 'k-', 'LineWidth', 1);              % horizontal bar
        plot([g1 g1], [yB, yB+tickH], 'k-', 'LineWidth', 1);       % left tick  ↓ toward data
        plot([g2 g2], [yB, yB+tickH], 'k-', 'LineWidth', 1);       % right tick ↓ toward data
        text((g1+g2)/2, yB-tickH, sigStars(pVal), ...              % label above bar
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize', 9, 'FontWeight', 'bold');
    end

    % Expand ylim to reveal brackets above the shallowest data point
    ylim([yMin - (size(comps,1)+1)*bracketStep,  yMax + bracketStep]);

    % Animal legend
    lgdH = gobjects(nAnimals, 1);
    for ai = 1:nAnimals
        lgdH(ai) = plot(nan, nan, 'o', ...
            'Color', animalCmap(ai,:), 'MarkerFaceColor', animalCmap(ai,:), ...
            'MarkerSize', 6, 'DisplayName', allAnimalNames{ai});
    end
    legend(lgdH, 'Location','best', 'FontSize',7, 'Box','off');

    set(gca, ...
        'XTick',      1:nG, ...
        'XTickLabel', groupOrder, ...
        'YDir',       'reverse', ...
        'XLim',       [0.5, nG+0.5], ...
        'FontSize',   8, ...
        'FontName',   'helvetica');
    ylabel('Cortical depth (\mum)', 'FontSize', 9);
    xtickangle(20);
    title('Depth: stripe vs responsive non-stripe', 'FontSize', 9);
    box on;  hold off;

    depthPdf = fullfile(saveDir, ...
        sprintf('StripeNeurons_depth_%s.pdf', params.indexType));
    exportgraphics(figDepth, depthPdf, 'ContentType','vector');
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
%  LOCAL FUNCTIONS
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
    bA = hierBootMatchFreq(rowsA.depth, nBoot, double(rowsA.insertion), double(rowsA.animal));
    bB = hierBootMatchFreq(rowsB.depth, nBoot, double(rowsB.insertion), double(rowsB.animal));
    pOne = mean(bA >= bB);
    pVal = 2 * min(pOne, 1 - pOne);    % two-tailed
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