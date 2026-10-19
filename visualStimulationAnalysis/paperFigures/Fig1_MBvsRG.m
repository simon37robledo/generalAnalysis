%% FIGURE 1 — Moving ball (MB) vs static ball (SB, rect grid RG)
% Statistics design: figureDesign.m (shared by every figure). Run the cells in
% order; rasters must be drawn in an interactive MATLAB session (they render
% blank under matlab -batch).

d     = figureDesign();                                           % shared design
paperFigPrefix('fig1');                                           % saved panels are named fig1_<panel>
stims = {'MB', 'RG'};                                             % compared stimuli
fig1Struct = d.structPath(d.exps, stims);                         % analysis struct written by the stats cell

%% Statistics (permutation test, swarm panels, analysis struct)
markUnits = { 1, 'MB', 'PV233', 6 ;                               % -> X in the MB column
             14, 'RG', 'PV140', 9 };                              % -> + in the RG column
args = [d.statsArgs, {'ComparePairs', stims, 'Lock', d.lock(stims), 'markUnits', markUnits}];  % full call
[tempTableMW, summaryStats] = AllExpAnalysis(d.exps, args{:});    % per-recording stats + pooled comparison

%% Example neurons (from the struct written above)
T = extractExampleNeurons(fig1Struct, plot=true, nPer=15, PaperFig=true);

%% Leave-one-animal-out sensitivity
Tloao = leaveOneAnimalOut(fig1Struct);                            % writes LeaveOneAnimalOut.xlsx/.csv next to the struct

%% Supplementary CSV tables (Paper_figs\Stats summary tables\)
statsTablePaths = makeStatsSummaryTables();                       % needs the struct and the LOAO .xlsx above

%% Population PSTH (MB aligned to the RF crossing, RG to onset)
plotPSTH_MultiExp(d.exps, d.psthArgs{:}, overwrite=true, stimTypes={"MB","RG"}, ...
    unionResponsive=true, requireAllStims=true);

%% Population raster (interactive session) — RG panel shows the ON period only
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes=["MB","RG"], unionUnits=true);

%% FP/FN screen on the raster cache
D = findFalseNegAndPos('W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_1-97_n44_Raster_MB-RG_union.mat');
Dfn = sortrows(D(D.FN_for ~= "", :), 'Severity', 'descend');      % strongest missed responses first
plotFalseNegPosNeurons(Dfn(1:min(15, height(Dfn)), :));           % inspect the top candidates

%% Supplementary: RG OFF responses, aligned to each recording's stimulus offset
plotPSTH_MultiExp(d.exps, d.psthArgs{:}, overwrite=true, stimTypes={"RG"}, alignTo="offset");
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes="RG", unionUnits=false, alignTo="offset");

%% Example single-neuron rasters and receptive fields
for ex = [74]                                                     % RG example recording
    NP = loadNPclassFromTable(ex);
    vsRe = rectGridAnalysis(NP);
    vsRe.plotRaster(MergeNtrials=1, overwrite=true, AllResponsiveNeurons=false, exNeuronsPhyID=[], exNeurons=111, selectedLum=255, oneTrial=true, PaperFig=true)
    colorbarLimsRG = vsRe.PlotReceptiveFields(exNeurons=111, allStimParamsCombined=false, PaperFig=true, overwrite=true);
end
for ex = [97]                                                     % MB example recording
    NP = loadNPclassFromTable(ex);
    vs = linearlyMovingBallAnalysis(NP);
    vs.plotRaster('exNeuronsPhyID', 357, 'overwrite', true, 'MergeNtrials', 1, 'PaperFig', true, OneLuminosity='white', OneDirection='up', bin=20)
    colorbarLims = vs.PlotReceptiveFields('exNeurons', 1, 'overwrite', true, 'OneDirection', 'up', 'OneLuminosity', 'white', 'PaperFig', true);
end
