%% SUPPLEMENTARY — Figure 2 spatial tuning index at top-cell thresholds of 5-30%
% The index compares the mean response of the top X% of the 9x9 receptive-field cells
% with the remaining cells (L_amplitude_diff, SpatialTuningIndex/tuningIndices). Figure 2
% uses X = 5; earlier runs used 30. To show the result does not depend on that choice,
% the index is recomputed for X = 5, 10, ..., 30 with every other setting identical to
% Figure 2 (same recordings and arguments, read from Fig2_SpatialTuning.m's header cell).
%   unpaired: all MB-responsive vs all SB-responsive units (one panel per threshold;
%             one-tailed P(SB >= MB), two-sample hierarchical bootstrap, as Figure 2's LOAO)
%   paired  : units responsive to both, MB - SB per unit (the same stats struct Figure 2
%             reports from: mean, 95% CI, one-tailed p)
% Figure 2's per-unit caches (read by analyzeStripeNeurons and Figure 6) do not encode the
% threshold, so they are backed up before the sweep and restored, byte-identical, after it.

pf  = fileparts(mfilename('fullpath'));                           % paperFigures folder
txt = splitlines(string(fileread(fullfile(pf, 'Fig2_SpatialTuning.m'))));   % Figure 2 script
c2  = find(startsWith(strtrim(txt), "%%"), 2);                     % first two cell headers
eval(strjoin(txt(1:c2(2)-1), newline));                           % Figure 2 header: d, exps, tuneArgs
paperFigPrefix('sfig2');                                          % Supplementary Figure 2 panels: sfig2_<panel>
tps   = 5:5:30;                                                   % thresholds (% of grid cells)
nBoot = 10000;                                                    % bootstrap replicates (as Figure 2)
setArg = @(a, name, v) [a(1:find(strcmp(a(1:2:end), name))*2-1), {v}, a(find(strcmp(a(1:2:end), name))*2+1:end)];   % replace one name-value

%% Back up Figure 2's per-unit caches
cacheDir = fullfile(extractBefore(d.paperFigs, 'Paper_figs'), 'lizards', 'Combined_lizard_analysis');   % W:\...\Combined_lizard_analysis
cacheF = string(fullfile(cacheDir, compose('Ex_%d-%d_SpatialTuningIndex_linearlyMovingBall-rectGrid_RF_prefDir%s_splitHalf_presence.mat', ...
    exps(1), exps(end), ["_allResp"; ""])));                      % unpaired and paired caches (2x1 string)
bakDir = fullfile(tempdir, 'FigS_TuningIndexThreshold_cacheBackup');   % backup folder
if ~isfolder(bakDir), mkdir(bakDir); end
for k = 1:numel(cacheF)
    assert(isfile(cacheF(k)), 'Cache not found: %s (run Fig2_SpatialTuning.m first).', cacheF(k));   % Figure 2 must exist
    copyfile(cacheF(k), bakDir);                                  % keep the Figure 2 version
end

%% Sweep over thresholds (caches are restored below even if the sweep fails)
R = table();                                                      % one row per threshold
try
for tp = tps
    a = setArg(tuneArgs, 'topPercent', tp);                       % Figure 2 arguments, this threshold
    % unpaired: all responsive units per stimulus (panel saved)
    U = SpatialTuningIndex(exps, a{:}, overwrite=true, PaperFig=true, plotRFs=false, ...
        plotRFunion=false, detectStripe=false);                   % stripes do not depend on the threshold
    t  = U.tbl(U.tbl.onOff == 1 & U.tbl.sizeIdx == 1 & U.tbl.lumValue == 255, :);   % plotted condition (white)
    mb = t(t.stimulus == "linearlyMovingBall", :);  sb = t(t.stimulus == "rectGrid", :);
    rng(0);                                                       % reproducible bootstrap
    bMB = hierBootMatchFreq(mb.L_amplitude_diff, nBoot, double(mb.animal), double(mb.experimentNum));   % (data, nrep, ANIMAL, INSERTION)
    bSB = hierBootMatchFreq(sb.L_amplitude_diff, nBoot, double(sb.animal), double(sb.experimentNum));
    pU  = mean(bSB >= bMB);                                       % one-tailed P(SB >= MB), as Figure 2
    % paired: units responsive to both (stats struct as Figure 2 reports it)
    ap = setArg(a, 'allResponsive', false);                       % intersection mode
    SpatialTuningIndex(exps, ap{:}, overwrite=true, PaperFig=true, plotRFs=false, ...
        plotRFunion=false, detectStripe=false);                   % writes the paired stats struct
    pairName = sprintf('SpatialTuningIndex-L_amplitude_diff-linearlyMovingBall-rectGrid-RF1-prefDir1-allResp0-tp-%d', tp);
    P = load(fullfile(d.paperFigs, [pairName '_statsStruct.mat'])); % saveStatsStruct output
    D = P.blocks{1,2}.diffs;                                      % MB - SB row
    assert(height(D) == 1 && D.group == "MB-SB", 'Threshold %d: expected one MB-SB difference row.', tp);
    R = [R; table(tp, height(mb), height(sb), mean(mb.L_amplitude_diff), mean(sb.L_amplitude_diff), pU, ...
        D.n, D.nAnimals, D.mean, D.ciLo, D.ciHi, D.p, ...
        'VariableNames', {'topPercent','nMB','nSB','meanMB','meanSB','pUnpaired_SBgeMB', ...
        'nPaired','nAnimalsPaired','pairedMeanMBminusSB','pairedCiLo','pairedCiHi','pPaired_MBgtSB'})];   %#ok<AGROW>
    close all;                                                    % free memory between thresholds
end
catch ME
    restoreCaches(cacheF, bakDir);                                % never leave Figure 2's caches modified
    rethrow(ME);                                                  % then fail loudly
end
restoreCaches(cacheF, bakDir);                                    % Figure 2 version back, byte-identical
disp(R);
writetable(R, fullfile(d.paperFigs, 'Stats summary tables', 'Fig2_tuning_index_by_threshold.csv'));   % supplementary table

%% Local functions
function restoreCaches(cacheF, bakDir)
% Copy the backed-up Figure 2 caches back and check they are byte-identical.
for k = 1:numel(cacheF)
    [~, nm, ex] = fileparts(cacheF(k));                           % cache file name
    bak = fullfile(bakDir, nm + ex);                              % its backup
    copyfile(bak, cacheF(k), 'f');                                % put the Figure 2 version back
    assert(isequal(readBytes(cacheF(k)), readBytes(bak)), 'Cache restore failed: %s', cacheF(k));   % fail loudly if not identical
end
fprintf('Figure 2 caches restored (byte-identical).\n');
end

function b = readBytes(f)
% Raw bytes of a file (binary-safe comparison).
fid = fopen(f, 'r'); assert(fid > 0, 'Cannot open %s.', f);       % open read-only
b = fread(fid, Inf, '*uint8'); fclose(fid);                       % all bytes
end
