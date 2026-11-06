%% FIGURE 3 — Moving vs static global stimuli (drifting vs static grating)
% Statistics design: figureDesign.m. SDGm is tested against the blank pre-trial
% baseline (not the static-grating phase). Category analyses (angles, temporal
% and spatial frequency) come in the next round
% (StatisticsPerNeuronPerCategory still needs the same fixes). NV is not in the paper.
% Recordings: d.expsSDG (recordings 1 and 8 excluded: different grating protocol, figureDesign).

d     = figureDesign();                                           % shared design
paperFigPrefix('fig3');                                           % saved panels are named fig3_<panel>
stims = {'SDGm', 'SDGs'};                                         % compared stimuli
fig3Struct = d.structPath(d.expsSDG, stims);                         % analysis struct written by the stats cell

%% Statistics
args = [d.statsArgs, {'ComparePairs', stims, 'Lock', d.lock(stims)}];
[tempTableMW, summaryStats] = AllExpAnalysis(d.expsSDG, args{:});    % SDGm vs SDGs

%% Leave-one-animal-out sensitivity
Tloao = leaveOneAnimalOut(fig3Struct, figLabel='Figure 3', outName='LeaveOneAnimalOut_SDGm-SDGs');

%% Population PSTH
plotPSTH_MultiExp(d.expsSDG, d.psthArgsSDG{:}, overwrite=true, stimTypes={"SDGm","SDGs"}, unionResponsive=true);

%% Population raster (interactive session): static then drifting phase of each trial
plotRaster_MultiExp(d.expsSDG, d.rasterArgs{:}, overwrite=true, stimTypes="SDG", unionUnits=true, ...
    useCompleteWindow=true, preBase=500);
