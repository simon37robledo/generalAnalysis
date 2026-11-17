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

%% Selection-free check (2026-11-13): every unit present in both blocks, no responsiveness filter
% The main panels keep units responsive to either grating (92 SDGm vs 8 SDGs), which builds part
% of SDGm > SDGs into the unit set. Same test on all present units; outputs carry '_allPresent'.
[tempTableAll, summaryStatsAll] = AllExpAnalysis(d.expsSDG, args{:}, allPresentUnits=true);   % all present units
fig3StructAll = strrep(fig3Struct, '_analysisStruct.mat', '_allPresent_analysisStruct.mat');    % its analysis struct
TloaoAll = leaveOneAnimalOut(fig3StructAll, figLabel='Figure 3 (all present units)', outName='LeaveOneAnimalOut_SDGm-SDGs_allPresent');

%% Population PSTH
plotPSTH_MultiExp(d.expsSDG, d.psthArgsSDG{:}, overwrite=true, stimTypes={"SDGm","SDGs"}, unionResponsive=true);

%% Population raster (interactive session): split-half trials, one panel per test
% Static grating 0-1 s from static onset, drifting grating 0-2 s from motion onset
% (each phase's own test window and selection). The 200 ms baseline fits the
% shortest inter-trial interval (400 ms, recordings 40-43); the old combined
% static-then-drifting panel (preBase 500) cannot show split-half trials, because
% the two tests select their best condition separately.
sdgW = d.statsArgs{find(strcmp(d.statsArgs, 'sdgWindows')) + 1};   % [static drifting] test windows (ms)
plotRaster_MultiExp(d.expsSDG, d.rasterArgs{:}, overwrite=true, stimTypes=["SDGm","SDGs"], unionUnits=true, ...
    postStimPer=struct('SDGs', sdgW(1), 'SDGm', sdgW(2)));         % drifting | static, as every other panel (was SG | MG);
                                                                  % each drawn over its test window (SDGm 1923 ms, was 2000)
