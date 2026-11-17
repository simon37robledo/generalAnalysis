%% FIGURE 7 — Laudakia vulgaris (SA): replication of the movement and locality results
% Same statistics design as the Pogona figures (figureDesign: windows, white only,
% presence >= 0.7, max-statistic permutation test, split-half delta rate, hierarchical
% bootstrap). Units: BombCell 'good' with PV's saved parameter set (relabelled 2026-10-31;
% SA5 re-sorted with Kilosort4 and labelled 2026-11-04).
% Recordings (inventory 2026-11-03, onset-TTL lags normal: MB 69-75 ms, SB 46-52, SDG 88):
%   MB vs SB : 28-32 (SA5), 44, 45, 47, 48 (SA8), 56 (SA9), 98 (SA11) — 4 animals, 11 insertions
%              (46 has no MB; 46 and 55 use a non-grid static protocol). SA11 (98) was dropped on
%              2026-11-12 (post hoc, 1 unit entered the comparison, 7 of 32 present in both blocks) and
%              put back on 2026-11-18 (user): removing a post-hoc exclusion.
% SA5 labels use minSNR 8 instead of 18 (user, 2026-11-12, post hoc): its NP 1.0 probe gives lower SNR
% (median 14 vs 41 on NP 2.0); 8 gives SA5 the pass rate NP 2.0 recordings have at 18 (91%).
%   gratings : 28-32 (SA5) and 56 (SA9 ins 2; 55's grating sync fails) — 2 animals, 6 insertions:
%              small n, the leave-one-animal-out rows are each animal alone (consistency, not inference).
% SA5 (2023 protocol, Neuropixels 1.0 probe): single-trajectory MB (8 directions, 100 px; grid test
% as in PV recordings 1 and 8), 4x4 SB, gratings at 1 Hz (square wave, angles 0:45:360, so 0 and
% 360 are the same direction and stay separate categories). Its trials are timed by TTL + the
% recording's measured diode lag (getSyncedDiodeTriggers 'trialOnsetFromTTL', cached on Y:; lags
% FFF 39-44, SB 41-42, MB 73-74, SDG static onset 78-79 ms).
% Comparisons: MB vs SB (Fig 1); SDGm vs SDGs (movement, Fig 3); MB vs SDGm and
% SB vs SDGs (locality, Fig 4). Metrics: z-score and fraction responsive.

d = figureDesign();                                               % shared design (same as Pogona)
paperFigPrefix('fig7');                                           % saved panels are named fig7_<panel>
expsMBSB = [28:32 44 45 47 48 56 98];                            % SA5 (single-trajectory MB, 4x4 SB) + 9x9 MB and SB, sync checked; SA11 (98) back in (2026-11-18)
expsSDG  = [28:32 56];                                            % SA recordings with usable gratings (SA5, SA9 ins 2)

%% MB vs SB (4 animals): statistics and leave-one-animal-out
stims = {'MB', 'RG'};
AllExpAnalysis(expsMBSB, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims));
leaveOneAnimalOut(d.structPath(expsMBSB, stims), figLabel='Figure 7 MB-SB', outName='LeaveOneAnimalOut_SA_MB-RG');
AllExpAnalysis(expsMBSB, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims), allPresentUnits=true);   % selection-free check (2026-11-17): every unit present in both blocks
leaveOneAnimalOut(strrep(d.structPath(expsMBSB, stims), '_analysisStruct.mat', '_allPresent_analysisStruct.mat'), ...
    figLabel='Figure 7 MB-SB (all present units)', outName='LeaveOneAnimalOut_SA_MB-RG_allPresent');   % its own LOAO table

%% Gratings (2 animals, small n): movement and locality
for pr = {{'SDGm', 'SDGs'}, {'MB', 'SDGm'}, {'RG', 'SDGs'}}
    stims = pr{1};                                                % compared pair
    AllExpAnalysis(expsSDG, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims));   % hierarchical bootstrap over 2 animals
    leaveOneAnimalOut(d.structPath(expsSDG, stims), figLabel=sprintf('Figure 7 %s-%s', stims{:}), ...
        outName=sprintf('LeaveOneAnimalOut_SA_%s-%s', stims{:}));  % with 2 animals: each animal alone
    AllExpAnalysis(expsSDG, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims), allPresentUnits=true);   % selection-free check (2026-11-17)
    leaveOneAnimalOut(strrep(d.structPath(expsSDG, stims), '_analysisStruct.mat', '_allPresent_analysisStruct.mat'), ...
        figLabel=sprintf('Figure 7 %s-%s (all present units)', stims{:}), ...
        outName=sprintf('LeaveOneAnimalOut_SA_%s-%s_allPresent', stims{:}));   % its own LOAO table
end

%% Sensitivity (2026-11-18): units present in each block separately, compared unpaired
% Probe drift loses units between blocks; here each stimulus keeps the units present during its
% own block (two-sample hierarchical bootstrap). The paired design above stays primary.
cmpAll = {{'MB','RG'}, expsMBSB; {'SDGm','SDGs'}, expsSDG; {'MB','SDGm'}, expsSDG; {'RG','SDGs'}, expsSDG};
for k = 1:size(cmpAll, 1)
    stims = cmpAll{k,1};                                          % compared pair
    AllExpAnalysis(cmpAll{k,2}, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims), presence="perStimulus");   % responsive units
    AllExpAnalysis(cmpAll{k,2}, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims), presence="perStimulus", ...
        allPresentUnits=true);                                    % every unit present in its block
end

%% Per-animal MB - SB (and grating) differences, Pogona and Laudakia on one axis (descriptive)
pairsF = {{'MB','RG'}, d.exps, expsMBSB, 'Z-score', 'MB - SB z-score', 'MB-SB_Zscore';
          {'MB','RG'}, d.exps, expsMBSB, 'DeltaRate', 'MB - SB \Delta rate (spikes/s)', 'MB-SB_DeltaRate';
          {'SDGm','SDGs'}, d.expsSDG, expsSDG, 'Z-score', 'MG - SG z-score', 'SDGm-SDGs_Zscore';
          {'MB','SDGm'}, d.expsSDG, expsSDG, 'Z-score', 'MB - MG z-score', 'MB-SDGm_Zscore';
          {'RG','SDGs'}, d.expsSDG, expsSDG, 'DeltaRate', 'SB - SG \Delta rate (spikes/s)', 'SB-SDGs_DeltaRate'};
for k = 1:size(pairsF, 1)
    perAnimalForest([string(d.structPath(pairsF{k,2}, pairsF{k,1})), string(d.structPath(pairsF{k,3}, pairsF{k,1}))], ...
        ["Pogona" "Laudakia"], pairsF{k,1}, pairsF{k,4}, xLabel=pairsF{k,5}, outName=['PerAnimal_' pairsF{k,6}]);   % fig7_PerAnimal_*.jpg/.csv
end

%% Population PSTHs (batch): split-half trials, as Figures 1, 3 and 4
sdgW = d.statsArgs{find(strcmp(d.statsArgs, 'sdgWindows')) + 1};   % [static drifting] test windows (ms)
psthMB = d.setArg(d.psthArgs, 'postStim', 2300);                  % MB onset segment (ms), as Figure 1
stemMB = fullfile(d.paperFigs, [paperFigPrefix() 'PSTH-MB-RG_allStims_cv_mbOnset_onOff']);   % panel name (no recording range in it)
anNames = ["SA5", "SA8", "SA9", "SA11"];                          % per animal (why SA8 differs); drawn first, then renamed,
anExps  = {28:32, 44:48, 56, 98};                                 % because the pooled panel below reuses the same file name
for k = 1:numel(anNames)
    if isfile([stemMB '.jpg']), delete([stemMB '.jpg']); end      % remove the stale copy so a failed save cannot pass
    exA = intersect(expsMBSB, anExps{k});                         % this animal's MB vs SB recordings
    plotPSTH_MultiExp(exA, psthMB{:}, overwrite=true, stimTypes={"MB","RG"}, ...
        unionResponsive=true, requireAllStims=true, splitPanels=true, capMs=struct('RG', 500), offsetSegMs=500);
    assert(isfile([stemMB '.jpg']), 'Fig7: per-animal PSTH %s was not saved.', anNames(k));   % fail loudly
    movefile([stemMB '.jpg'], [stemMB '_' char(anNames(k)) '.jpg'], 'f');   % keep it as ..._SA5.jpg etc.
end
plotPSTH_MultiExp(expsMBSB, psthMB{:}, overwrite=true, stimTypes={"MB","RG"}, ...
    unionResponsive=true, requireAllStims=true, splitPanels=true, capMs=struct('RG', 500), offsetSegMs=500);   % MB vs SB, all animals
plotPSTH_MultiExp(expsSDG, d.psthArgsSDG{:}, overwrite=true, stimTypes={"SDGm","SDGs"}, unionResponsive=true);   % MG vs SG
psthMG = d.setArg(d.setArg(d.psthArgsSDG, 'postStim', 2400), 'capMs', struct('SDGs', sdgW(1), 'SDGm', sdgW(2)));
plotPSTH_MultiExp(expsSDG, psthMG{:}, overwrite=true, stimTypes={"MB","SDGm"}, unionResponsive=true, requireAllStims=true);   % MB vs MG
plotPSTH_MultiExp(expsSDG, d.psthArgsStatic{:}, overwrite=true, stimTypes={"RG","SDGs"}, unionResponsive=true, requireAllStims=true);   % SB vs SG

%% Population rasters (interactive session): one row per unit, shared z scale
plotRaster_MultiExp(expsMBSB, d.rasterArgs{:}, overwrite=true, stimTypes=["MB","RG"], unionUnits=true, ...
    postStimPer=struct('MB', 2300, 'RG', 500), offsetSegMs=500);  % as Figure 1
plotRaster_MultiExp(expsSDG, d.rasterArgs{:}, overwrite=true, stimTypes=["SDGm","SDGs"], unionUnits=true, ...
    postStimPer=struct('SDGs', sdgW(1), 'SDGm', sdgW(2)));         % as Figure 3
plotRaster_MultiExp(expsSDG, d.rasterArgsSDG{:}, overwrite=true, stimTypes=["MB","SDGm"], unionUnits=true, ...
    postStimPer=struct('MB', 2400, 'SDGm', sdgW(2)));              % as Figure 4
plotRaster_MultiExp(expsSDG, d.rasterArgsStatic{:}, overwrite=true, stimTypes=["RG","SDGs"], unionUnits=true);   % as Figure 4
