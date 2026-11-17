%% FIGURE 4 — Local vs global stimuli
% Statistics design: figureDesign.m.
% Static: SB (RG) vs static grating (SDGs), paired within neurons in every
%   recording of d.expsSDG holding both (recordings 1 and 8 excluded from every grating
%   comparison: different grating protocol, see figureDesign). FFF was dropped from this contrast:
%   it exists in only 19 recordings and 3 units respond to it (chance level), so
%   the 3-way design reduced the analysable set to 4 units.
% Moving: MB vs moving bar (MBR) and MB vs drifting grating (SDGm) as the
%   inferential tests (MBR and SDGm co-occur only in PV97). MBR vs SDGm is a
%   descriptive within-PV97 panel only (one animal: no animal-level claim).

d = figureDesign();                                               % shared design
paperFigPrefix('fig4');                                           % saved panels are named fig4_<panel>

%% Static locality: RG vs SDGs (statistics)
stimsS = {'RG', 'SDGs'};                                          % pairwise: more units than the 3-way with FFF
args = [d.statsArgs, {'ComparePairs', stimsS, 'Lock', d.lock(stimsS)}];
AllExpAnalysis(d.expsSDG, args{:});
leaveOneAnimalOut(d.structPath(d.expsSDG, stimsS), figLabel='Figure 4 static', outName='LeaveOneAnimalOut_RG-SDGs');
AllExpAnalysis(d.expsSDG, args{:}, allPresentUnits=true);       % selection-free check (2026-11-17): every unit present in both blocks
leaveOneAnimalOut(strrep(d.structPath(d.expsSDG, stimsS), '_analysisStruct.mat', '_allPresent_analysisStruct.mat'), ...
    figLabel='Figure 4 static (all present units)', outName='LeaveOneAnimalOut_RG-SDGs_allPresent');   % its own LOAO file

%% Static locality: PSTH and raster (RG ON period; interactive session for the raster)
plotPSTH_MultiExp(d.expsSDG, d.psthArgsStatic{:}, overwrite=true, stimTypes={"RG","SDGs"}, unionResponsive=true, requireAllStims=true);
plotRaster_MultiExp(d.expsSDG, d.rasterArgsStatic{:}, overwrite=true, stimTypes=["RG","SDGs"], unionUnits=true);   % 0-1 s

%% Moving locality: MB vs MBR (statistics)
stimsB = {'MB', 'MBR'};
args = [d.statsArgs, {'ComparePairs', stimsB, 'Lock', d.lock(stimsB)}];
AllExpAnalysis(d.exps, args{:});
leaveOneAnimalOut(d.structPath(d.exps, stimsB), figLabel='Figure 4 MB-MBR', outName='LeaveOneAnimalOut_MB-MBR');
AllExpAnalysis(d.exps, args{:}, allPresentUnits=true);          % selection-free check (2026-11-17)
leaveOneAnimalOut(strrep(d.structPath(d.exps, stimsB), '_analysisStruct.mat', '_allPresent_analysisStruct.mat'), ...
    figLabel='Figure 4 MB-MBR (all present units)', outName='LeaveOneAnimalOut_MB-MBR_allPresent');   % its own LOAO file

%% Moving locality: MB vs SDGm (statistics)
stimsG = {'MB', 'SDGm'};
args = [d.statsArgs, {'ComparePairs', stimsG, 'Lock', d.lock(stimsG)}];
AllExpAnalysis(d.expsSDG, args{:});
leaveOneAnimalOut(d.structPath(d.expsSDG, stimsG), figLabel='Figure 4 MB-SDGm', outName='LeaveOneAnimalOut_MB-SDGm');
AllExpAnalysis(d.expsSDG, args{:}, allPresentUnits=true);       % selection-free check (2026-11-17)
leaveOneAnimalOut(strrep(d.structPath(d.expsSDG, stimsG), '_analysisStruct.mat', '_allPresent_analysisStruct.mat'), ...
    figLabel='Figure 4 MB-SDGm (all present units)', outName='LeaveOneAnimalOut_MB-SDGm_allPresent');   % its own LOAO file

%% Moving locality: PSTHs and rasters (all onset-aligned, split-half trials)
% MB, MBR and SDGm aligned to stimulus onset (SDGm: motion onset). The ball and the
% bar cross each RF at a different time, so population PSTHs are spread over the sweep;
% the per-unit rasters show the timing.
paperFigPrefix('fig4');                                           % main-figure panel names
psthSweep = d.setArg(d.psthArgs, 'postStim', 2400);               % MB / MBR sweeps (~2.4 s)
plotPSTH_MultiExp(d.exps, psthSweep{:}, overwrite=true, stimTypes={"MB","MBR"}, unionResponsive=true, requireAllStims=true);
sdgW = d.statsArgs{find(strcmp(d.statsArgs, 'sdgWindows')) + 1};   % [static drifting] test windows (ms)
psthMG = d.setArg(d.setArg(d.psthArgsSDG, 'postStim', 2400), 'capMs', struct('SDGs', sdgW(1), 'SDGm', sdgW(2)));  % MB 0-2.4 s, MG over its 1923 ms test window (was 2000)
plotPSTH_MultiExp(d.expsSDG, psthMG{:}, overwrite=true, stimTypes={"MB","SDGm"}, unionResponsive=true, requireAllStims=true);
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes=["MB","MBR"], unionUnits=true, useCompleteWindow=true);  % whole sweeps + 500 ms
plotRaster_MultiExp(d.expsSDG, d.rasterArgsSDG{:}, overwrite=true, stimTypes=["MB","SDGm"], unionUnits=true, ...
    postStimPer=struct('MB', 2400, 'SDGm', sdgW(2)));              % MB 0-2.4 s, MG over its 1923 ms test window (was 2000)

%% Supplementary Figure 4: MB aligned to the crossing of its selected grid cell (split-half)
paperFigPrefix('sfig4');                                          % supplementary panel names
psthCrossB = d.setArg(d.psthArgs, 'mbAlign', "crossing");         % as the 6 Nov panels, split-half trials
psthCrossG = d.setArg(d.psthArgsSDG, 'mbAlign', "crossing");
plotPSTH_MultiExp(d.exps, psthCrossB{:}, overwrite=true, stimTypes={"MB","MBR"}, unionResponsive=true, requireAllStims=true);
plotPSTH_MultiExp(d.expsSDG, psthCrossG{:}, overwrite=true, stimTypes={"MB","SDGm"}, unionResponsive=true, requireAllStims=true);
rasterCrossB = d.setArg(d.rasterArgs, 'mbAlign', "crossing");
rasterCrossG = d.setArg(d.rasterArgsSDG, 'mbAlign', "crossing");
plotRaster_MultiExp(d.exps, rasterCrossB{:}, overwrite=true, stimTypes=["MB","MBR"], unionUnits=true, useCompleteWindow=true);
plotRaster_MultiExp(d.expsSDG, rasterCrossG{:}, overwrite=true, stimTypes=["MB","SDGm"], unionUnits=true);  % MB 500 ms after the crossing
paperFigPrefix('fig4');                                           % back to main-figure names

%% Supplementary (descriptive, PV97 only): MBR vs SDGm
expsPV97 = 49:54;                                                 % the only recordings with both MBR and SDG
stimsP = {'MBR', 'SDGm'};
args = [d.statsArgs, {'ComparePairs', stimsP, 'Lock', d.lock(stimsP)}];
AllExpAnalysis(expsPV97, args{:});                                % one animal: insertion-level resampling only
