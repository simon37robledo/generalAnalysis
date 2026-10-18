%% FIGURE 4 — Local vs global stimuli
% Statistics design: figureDesign.m.
% Static: 3-way SB (RG) vs static grating (SDGs) vs full-field flash (FFF); every
%   FFF recording also has RG and SDG, so all three contrasts are paired within
%   the same neurons (19 recordings / 4 animals).
% Moving: MB vs moving bar (MBR) and MB vs drifting grating (SDGm) as the
%   inferential tests (MBR and SDGm co-occur only in PV97). MBR vs SDGm is a
%   descriptive within-PV97 panel only (one animal: no animal-level claim).

d = figureDesign();                                               % shared design

%% Static locality: RG vs SDGs vs FFF (statistics)
stimsS = {'RG', 'SDGs', 'FFF'};
args = [d.statsArgs, {'ComparePairs', stimsS, 'Lock', d.lock(stimsS)}];
AllExpAnalysis(d.exps, args{:});
leaveOneAnimalOut(d.structPath(d.exps, stimsS), figLabel='Figure 4 static', outName='LeaveOneAnimalOut_RG-SDGs-FFF');

%% Static locality: PSTH and raster (RG ON period; interactive session for the raster)
plotPSTH_MultiExp(d.exps, d.psthArgs{:}, overwrite=true, stimTypes={"RG","SDGs","FFF"}, unionResponsive=true, requireAllStims=true);
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes=["RG","SDGs","FFF"], unionUnits=true);

%% Moving locality: MB vs MBR (statistics)
stimsB = {'MB', 'MBR'};
args = [d.statsArgs, {'ComparePairs', stimsB, 'Lock', d.lock(stimsB)}];
AllExpAnalysis(d.exps, args{:});
leaveOneAnimalOut(d.structPath(d.exps, stimsB), figLabel='Figure 4 MB-MBR', outName='LeaveOneAnimalOut_MB-MBR');

%% Moving locality: MB vs SDGm (statistics)
stimsG = {'MB', 'SDGm'};
args = [d.statsArgs, {'ComparePairs', stimsG, 'Lock', d.lock(stimsG)}];
AllExpAnalysis(d.exps, args{:});
leaveOneAnimalOut(d.structPath(d.exps, stimsG), figLabel='Figure 4 MB-SDGm', outName='LeaveOneAnimalOut_MB-SDGm');

%% Moving locality: PSTHs and rasters
% MB is aligned to the RF crossing; MBR and SDGm to stimulus onset. The bar
% crosses each RF at a different time, so the MBR population PSTH is smeared;
% the per-unit raster over the whole sweep (useCompleteWindow) shows it better.
plotPSTH_MultiExp(d.exps, d.psthArgs{:}, overwrite=true, stimTypes={"MB","MBR"}, unionResponsive=true, requireAllStims=true);
plotPSTH_MultiExp(d.exps, d.psthArgs{:}, overwrite=true, stimTypes={"MB","SDGm"}, unionResponsive=true, requireAllStims=true);
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes=["MB","MBR"], unionUnits=true, useCompleteWindow=true);
plotRaster_MultiExp(d.exps, d.rasterArgs{:}, overwrite=true, stimTypes=["MB","SDGm"], unionUnits=true);

%% Supplementary (descriptive, PV97 only): MBR vs SDGm
expsPV97 = 49:54;                                                 % the only recordings with both MBR and SDG
stimsP = {'MBR', 'SDGm'};
args = [d.statsArgs, {'ComparePairs', stimsP, 'Lock', d.lock(stimsP)}];
AllExpAnalysis(expsPV97, args{:});                                % one animal: insertion-level resampling only
