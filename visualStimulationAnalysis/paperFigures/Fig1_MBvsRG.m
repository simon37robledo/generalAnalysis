%% FIGURE 1 — Moving ball (MB) vs static ball (SB, rect grid RG)
% Statistics design: figureDesign.m (shared by every figure). Run the cells in
% order; rasters must be drawn in an interactive MATLAB session (they render
% blank under matlab -batch).

d     = figureDesign();                                           % shared design
paperFigPrefix('fig1');                                           % saved panels are named fig1_<panel>
stims = {'MB', 'RG'};                                             % compared stimuli
fig1Struct = d.structPath(d.exps, stims);                         % analysis struct written by the stats cell

%% Statistics (permutation test, swarm panels, analysis struct)
% Positions in the good-unit list (labels of 2026-11-05: BombCell without whole-recording presence):
markUnits = { 1, 'MB', 'PV233', 6 ;                               % -> X in the MB column (phy 357, rec 97)
             15, 'RG', 'PV140', 9 };                              % -> + in the RG column (phy 288, rec 81; was 14 with the old labels)
args = [d.statsArgs, {'ComparePairs', stims, 'Lock', d.lock(stims), 'markUnits', markUnits}];  % full call
[tempTableMW, summaryStats] = AllExpAnalysis(d.exps, args{:});    % per-recording stats + pooled comparison

%% Example neurons (from the struct written above)
T = extractExampleNeurons(fig1Struct, plot=true, nPer=15, PaperFig=true);

%% Leave-one-animal-out sensitivity
Tloao = leaveOneAnimalOut(fig1Struct);                            % writes LeaveOneAnimalOut.xlsx/.csv next to the struct

%% Supplementary CSV tables (Paper_figs\Stats summary tables\)
statsTablePaths = makeStatsSummaryTables();                       % needs the struct and the LOAO .xlsx above

%% Population PSTH: onset | offset, split-half trials (figureDesign)
% Two panels with a shared z-score axis. Each panel: aligned to stimulus onset (MB 0-2.3 s,
% the shortest sweep is 2.34 s; SB 0-500 ms, the shortest SB is ~490 ms), a 100 ms gap,
% then 500 ms aligned to each recording's own measured offset (SB lasts ~500 ms in 24
% recordings and ~1000 ms in 20; MB sweeps 2.34-2.39 s, 9.6 s at the slow speed of recs
% 49/51/53). Same held-out trials in both segments.
paperFigPrefix('fig1');                                           % main-figure panel names
psthMB = d.setArg(d.psthArgs, 'postStim', 2300);                  % onset segment (ms), within every MB sweep
plotPSTH_MultiExp(d.exps, psthMB{:}, overwrite=true, stimTypes={"MB","RG"}, ...
    unionResponsive=true, requireAllStims=true, splitPanels=true, capMs=struct('RG', 500), offsetSegMs=500);

%% Population raster (interactive session): onset | offset, split-half trials
paperFigPrefix('fig1');                                           % main-figure panel names
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes=["MB","RG"], unionUnits=true, ...
    postStimPer=struct('MB', 2300, 'RG', 500), offsetSegMs=500);  % onset segments as the PSTH, then 500 ms after each offset

%% FP/FN screen on the raster cache (diagnostic, not a panel)
% Uses the stored-trial raster cache of 6 Nov (crossing-aligned MB, stored preferred
% category), the one the screen was designed for.
D = findFalseNegAndPos('W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_1-97_n44_Raster_MB-RG_union.mat');
Dfn = sortrows(D(D.FN_for ~= "", :), 'Severity', 'descend');      % strongest missed responses first
plotFalseNegPosNeurons(Dfn(1:min(15, height(Dfn)), :));           % inspect the top candidates

%% Supplementary Figure 1: MB aligned to the crossing of its selected grid cell (split-half)
% The cell x direction is the one with the largest response on the other half of the
% repetitions; t = 0 is the ball's closest approach to that cell's centre.
paperFigPrefix('sfig1');                                          % supplementary panel names
psthCross   = d.setArg(d.psthArgs, 'mbAlign', "crossing");        % MB crossing-aligned, SB onset-aligned (0-500 ms)
rasterCross = d.setArg(d.rasterArgs, 'mbAlign', "crossing");
plotPSTH_MultiExp(d.exps, psthCross{:}, overwrite=true, stimTypes={"MB","RG"}, ...
    unionResponsive=true, requireAllStims=true);
plotRaster_MultiExp(d.exps, rasterCross{:}, overwrite=true, stimTypes=["MB","RG"], unionUnits=true);  % interactive session

%% Supplementary Figure 1: SB onset vs offset (SB-responsive units, split-half trials)
paperFigPrefix('sfig1');                                          % supplementary panel names
Tonoff = plotRG_OnOff(d.exps, PaperFig=true);                     % PSTH + raster sorted by ON - OFF (interactive session)

%% Supplementary Figure 1: SB OFF responses, aligned to each recording's stimulus offset
paperFigPrefix('sfig1');                                          % supplementary panel names
plotPSTH_MultiExp(d.exps, d.psthArgs{:}, overwrite=true, stimTypes={"RG"}, alignTo="offset");
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes="RG", unionUnits=false, alignTo="offset");

%% Example single-neuron rasters and receptive fields
paperFigPrefix('fig1');                                           % main-figure panel names
for ex = [74]                                                     % RG example recording
    NP = loadNPclassFromTable(ex);
    vsRe = rectGridAnalysis(NP);
    vsRe.plotRaster(MergeNtrials=1, overwrite=true, AllResponsiveNeurons=false, exNeuronsPhyID=768, exNeurons=[], selectedLum=255, oneTrial=true, PaperFig=true)   % phy 768 (was position 111)
    colorbarLimsRG = vsRe.PlotReceptiveFields(exNeurons=246, allStimParamsCombined=false, PaperFig=true, overwrite=true);   % phy 768: position 246 with the 2026-11-05 labels (was 111)
end
for ex = [97]                                                     % MB example recording
    NP = loadNPclassFromTable(ex);
    vs = linearlyMovingBallAnalysis(NP);
    vs.plotRaster('exNeuronsPhyID', 357, 'overwrite', true, 'MergeNtrials', 1, 'PaperFig', true, OneLuminosity='white', OneDirection='up', bin=20)
    colorbarLims = vs.PlotReceptiveFields('exNeurons', 1, 'overwrite', true, 'OneDirection', 'up', 'OneLuminosity', 'white', 'PaperFig', true);   % phy 357: position 1 (unchanged)
end
