%% FIGURE 3 — Moving vs static global stimuli (drifting vs static grating)
% Statistics design: figureDesign.m. SDGm is tested against the blank pre-trial
% baseline (not the static-grating phase). Category analyses (angles, temporal
% and spatial frequency) and the natural-image postures come in the next round
% (StatisticsPerNeuronPerCategory still needs the same fixes). NV is not in the paper.

d     = figureDesign();                                           % shared design
paperFigPrefix('fig3');                                           % saved panels are named fig3_<panel>
stims = {'SDGm', 'SDGs'};                                         % compared stimuli
fig3Struct = d.structPath(d.exps, stims);                         % analysis struct written by the stats cell

%% Statistics
args = [d.statsArgs, {'ComparePairs', stims, 'Lock', d.lock(stims)}];
[tempTableMW, summaryStats] = AllExpAnalysis(d.exps, args{:});    % SDGm vs SDGs

%% Leave-one-animal-out sensitivity
Tloao = leaveOneAnimalOut(fig3Struct, figLabel='Figure 3', outName='LeaveOneAnimalOut_SDGm-SDGs');

%% Population PSTH
plotPSTH_MultiExp(d.exps, d.psthArgsSDG{:}, overwrite=true, stimTypes={"SDGm","SDGs"}, unionResponsive=true);

%% Population raster (interactive session): static then drifting phase of each trial
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes="SDG", unionUnits=true, ...
    useCompleteWindow=true, preBase=500);
