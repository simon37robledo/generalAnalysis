function T = plotRG_OnOff(exList, params)
% plotRG_OnOff  Static ball (SB, rectGrid): onset- vs offset-aligned responses.
%
%   T = plotRG_OnOff(exList, Name=Value)
%
%   Units: SB-responsive (permutation p < alpha) and present throughout the SB
%   block. Per unit, trials are the split-half held-out trials of DeltaRateCV
%   (splitHalfMR): the position x {ON, OFF} window chosen on one half of the
%   repetitions, shown with the other half's trials, and vice versa. The same
%   trials are aligned to stimulus onset and to stimulus offset (onset + that
%   recording's stimulus duration), each after the same pre-trial baseline.
%
%   Figures (saved with printFig, i.e. with the current paperFigPrefix):
%     RG-OnOff-PSTH    onset vs offset PSTH: per-recording mean of units,
%                      then mean +/- SEM across recordings (as plotPSTH_MultiExp)
%     RG-OnOff-Raster  two panels, same rows (units), sorted by
%                      mean z(0-sortWinMs after onset) - mean z(0-sortWinMs after offset)
%
%   Z-score: (rate - mean baseline) / max(SD baseline, floor), with per-trial
%   baseline rates of all trials and a floor at the sdFloorPrctile percentile
%   of the non-zero baseline SDs (as plotRaster_MultiExp).
%
%   OUTPUT
%     T  table, one row per unit: exp, phyID, prefSeg (1 = ON, 2 = OFF window
%        of the responsiveness test), onMeanZ, offMeanZ, onMinusOff, and the
%        z-scored PSTHs (onZ, offZ)

arguments
    exList double                                                      % recording IDs
    params.binWidth       double  = 20                                  % bin width (ms)
    params.preBase        double  = 200                                 % pre-trial baseline shown (ms)
    params.postMs         double  = 500                                 % window after onset / offset (ms)
    params.sortWinMs      double  = 300                                 % window for the ON - OFF sorting index (ms, = test window)
    params.alpha          double  = 0.05                                % responsiveness threshold
    params.sdFloorPrctile double  = 25                                  % z-score SD floor percentile
    params.climPrctile    double  = 90                                  % raster colour limit percentile
    params.PaperFig       logical = false                               % save the figures
end
assert(mod(params.preBase, params.binWidth) == 0 && mod(params.postMs, params.binWidth) == 0, ...
    'plotRG_OnOff: preBase and postMs must be multiples of binWidth.');  % bins aligned to t = 0
edges = 0 : params.binWidth : params.preBase + params.postMs;          % bin edges (ms) of baseline | after event
nBins = numel(edges) - 1;                                              % bins per PSTH

onR = []; offR = []; bM = []; bS = [];                                 % per-unit rate PSTHs and baselines
expU = []; phyU = []; segU = [];                                       % unit identity and preferred window
vsLast = [];                                                           % an analysis object for printFig
for ex = exList
    NP = loadNPclassFromTable(ex);                                     % recording
    try
        obj = rectGridAnalysis(NP);                                    % SB analysis object
        noRG = isempty(obj.VST);                                       % stimulus not presented
    catch
        noRG = true;                                                   % no SB in this recording
    end
    if noRG
        fprintf('  exp %d: no SB - skipped.\n', ex);
        continue
    end
    Stats = obj.StatisticsPerNeuron;                                   % cached statistics (bare call; stale caches error)
    verifySplitHalfStore(Stats, sprintf('exp %d RG', ex));             % stored split reproduces DeltaRateCV
    assert(~isempty(Stats.ObsDiffSeg) && size(Stats.ObsDiffSeg, 3) == 2, ...
        'plotRG_OnOff: exp %d SB stats are not ON/OFF (offWindowStims must include RG).', ex);  % ON/OFF design required
    pres = statPresence(Stats, "RG", "");                              % present throughout the SB block
    units = find(Stats.pvalsResponse(:) < params.alpha & pres(:))';    % responsive and present
    p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);  % spike sorting
    lab    = string(p_sort.label');                                    % unit labels
    goodU  = p_sort.ic(:, lab == "good");                              % good units (columns match Stats)
    phyIDs = p_sort.phy_ID(lab == "good");                             % phy cluster IDs
    assert(size(goodU, 2) == numel(Stats.pvalsResponse), ...
        'plotRG_OnOff: exp %d has %d good units but %d p-values.', ex, size(goodU, 2), numel(Stats.pvalsResponse));  % alignment
    for u = units
        [MRon, br] = splitHalfMR(Stats, u, goodU(:, u), p_sort.t, "onset",  params.preBase, params.postMs);  % baseline | from onset
        MRoff      = splitHalfMR(Stats, u, goodU(:, u), p_sort.t, "offset", params.preBase, params.postMs);  % baseline | from offset
        onR  = [onR;  heldOutRate(MRon,  edges, params.binWidth)];     %#ok<AGROW> spikes/s
        offR = [offR; heldOutRate(MRoff, edges, params.binWidth)];     %#ok<AGROW> spikes/s
        bM = [bM; mean(br)]; bS = [bS; std(br)];                       %#ok<AGROW> across-trial baseline (all trials)
        expU = [expU; ex]; phyU = [phyU; phyIDs(u)]; segU = [segU; Stats.prefSeg(u)];  %#ok<AGROW> identity
    end
    fprintf('  exp %d: %d SB-responsive unit(s).\n', ex, numel(units));
    vsLast = obj;                                                      % keep for printFig
end
assert(~isempty(onR), 'plotRG_OnOff: no SB-responsive units in these recordings.');  % nothing to plot

% ---- Z-score with the data-driven SD floor (as plotRaster_MultiExp) ----
sdFloor = prctile(bS(bS > 0), params.sdFloorPrctile);                  % floor for silent baselines
den  = max(bS, sdFloor);                                               % denominator per unit
onZ  = (onR  - bM) ./ den;                                             % [units x bins]
offZ = (offR - bM) ./ den;                                             % [units x bins]

% ---- Sorting index: ON minus OFF, in the test window after each event ----
binLeft = edges(1:end-1) - params.preBase;                             % left edge of each bin re. event (ms)
inWin   = binLeft >= 0 & binLeft < params.sortWinMs;                   % bins of the test window
onMeanZ  = mean(onZ(:, inWin), 2);                                     % mean z after onset
offMeanZ = mean(offZ(:, inWin), 2);                                    % mean z after offset
onMinusOff = onMeanZ - offMeanZ;                                       % sorting index
[~, ord] = sort(onMinusOff, 'ascend');                                 % row 1 (bottom) = most OFF-dominated, top = most ON-dominated

T = table(expU, phyU, segU, onMeanZ, offMeanZ, onMinusOff, onZ, offZ, ...
    'VariableNames', {'exp','phyID','prefSeg','onMeanZ','offMeanZ','onMinusOff','onZ','offZ'});  % per-unit output
fprintf('SB ON vs OFF: %d units, %d recordings. Preferred test window: ON %d, OFF %d (descriptive).\n', ...
    height(T), numel(unique(expU)), sum(segU == 1), sum(segU == 2));
fprintf('  ON > OFF (index > 0): %d units; OFF > ON: %d units.\n', sum(onMinusOff > 0), sum(onMinusOff < 0));

% ---- PSTH: per-recording mean, then mean +/- SEM across recordings ----
recs = unique(expU);                                                   % contributing recordings
recOn = nan(numel(recs), nBins); recOff = nan(numel(recs), nBins);     % per-recording mean traces
for k = 1:numel(recs)
    recOn(k, :)  = mean(onZ(expU == recs(k), :), 1);                   % mean of this recording's units
    recOff(k, :) = mean(offZ(expU == recs(k), :), 1);
end
fig1 = figure; set(fig1, 'Units', 'centimeters', 'Position', [5 5 9 7]);  % PSTH figure
ax = axes(fig1); hold(ax, 'on');
xE = edges - params.preBase;                                           % bin edges re. event (ms)
cols = [0 0.447 0.741; 0.85 0.325 0.098];                              % onset blue, offset orange
names = {sprintf('Onset (n=%d)', numel(recs)), sprintf('Offset (n=%d)', numel(recs))};  % legend labels
D = {recOn, recOff}; h = gobjects(1, 2);                               % traces and handles
for c = 1:2
    m = mean(D{c}, 1); se = std(D{c}, 0, 1) / sqrt(size(D{c}, 1));     % mean and SEM across recordings
    [xs, us] = stairs(xE, [m + se, m(end) + se(end)]);                 % upper SEM outline (steps)
    [~,  ls] = stairs(xE, [m - se, m(end) - se(end)]);                 % lower SEM outline
    fill(ax, [xs' fliplr(xs')], [us' fliplr(ls')], cols(c, :), 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % SEM band
    h(c) = stairs(ax, xE, [m m(end)], 'Color', cols(c, :), 'LineWidth', 1.5);  % mean PSTH (bins as steps)
end
xline(ax, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');      % event (onset or offset)
xlabel(ax, 'Time re. SB onset / offset [ms]'); ylabel(ax, 'Z-score');
legend(ax, h, names, 'Location', 'best', 'FontSize', 7);
title(ax, sprintf('SB ON vs OFF: %d units, N = %d recordings', height(T), numel(recs)), 'FontSize', 9);
set(ax, 'FontName', 'Helvetica', 'FontSize', 8, 'TickDir', 'out', 'Box', 'off'); hold(ax, 'off');

% ---- Raster: same rows in both panels, sorted by ON - OFF ----
fig2 = figure; set(fig2, 'Units', 'centimeters', 'Position', [5 5 12 9]);  % raster figure
tl = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');  % two panels
allZ = [onZ(:); offZ(:)];                                              % shared colour scale
cHi = prctile(allZ, params.climPrctile); if cHi <= 0, cHi = 1; end     % upper colour limit (z)
tc = (edges(1:end-1) + params.binWidth / 2 - params.preBase) / 1000;   % bin centres (s)
lbl = {'Onset', 'Offset'}; Z = {onZ(ord, :), offZ(ord, :)};            % panels
for c = 1:2
    axc = nexttile(tl);                                                % panel
    imagesc(axc, tc, 1:height(T), Z{c}, [0 cHi]);                      % units x time
    colormap(axc, flipud(gray)); axc.YDir = 'normal';                  % white = baseline, row 1 at the bottom
    xline(axc, 0, 'k', 'LineWidth', 1.5);                              % event
    title(axc, sprintf('SB %s (n=%d)', lbl{c}, height(T)), 'FontSize', 8);
    xlabel(axc, sprintf('Time re. %s (s)', lower(lbl{c})));
    if c == 1, ylabel(axc, 'Unit (sorted by ON - OFF, ON-dominated at top)'); end
    set(axc, 'FontName', 'Helvetica', 'FontSize', 8, 'TickDir', 'out', 'Box', 'off');
end
cb = colorbar(axc); cb.Label.String = 'Z-score';                       % shared colour scale
nOnTop = sum(onMinusOff > 0);                                          % rows with ON > OFF
for c = 1:2
    yline(nexttile(tl, c), height(T) - nOnTop + 0.5, 'r-', 'LineWidth', 0.8);  % boundary between OFF- and ON-dominated rows
end

if params.PaperFig && ~isempty(vsLast)
    figure(fig1);                                                      % printFig prints the CURRENT figure: make the PSTH current
    vsLast.printFig(fig1, 'RG-OnOff-PSTH', PaperFig = true);           % name gets the current paperFigPrefix
    figure(fig2);                                                      % then the raster
    vsLast.printFig(fig2, 'RG-OnOff-Raster', PaperFig = true);
end
end

function r = heldOutRate(MR, edges, binWidth)
% Mean of the two held-out PSTHs (spikes/s), as DeltaRateCV averages both splits.
r = zeros(2, numel(edges) - 1);                                        % one PSTH per selection half
for h = 1:2
    st = repmat(1:size(MR{h}, 2), size(MR{h}, 1), 1);                  % column index of every 1 ms position
    r(h, :) = histcounts(st(logical(MR{h})), edges) / (binWidth * size(MR{h}, 1)) * 1000;  % spikes/s
end
r = mean(r, 1);                                                        % average of the two
end
