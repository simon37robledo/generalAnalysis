%% FIGURE 7 — Laudakia vulgaris (SA): replication of the movement and locality results
% Same statistics design as the Pogona figures (figureDesign: windows, white only,
% presence >= 0.7, max-statistic permutation test, split-half delta rate, hierarchical
% bootstrap). Units: BombCell 'good' with PV's saved parameter set (relabelled 2026-10-31;
% SA5 re-sorted with Kilosort4 and labelled 2026-11-04).
% Recordings (inventory 2026-11-03, onset-TTL lags normal: MB 69-75 ms, SB 46-52, SDG 88):
%   MB vs SB : 28-32 (SA5), 44, 45, 47, 48 (SA8), 56 (SA9), 98 (SA11) — 4 animals, 11 insertions
%              (46 has no MB; 46 and 55 use a non-grid static protocol)
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
expsMBSB = [28:32 44 45 47 48 56 98];                            % SA5 (single-trajectory MB, 4x4 SB) + 9x9 MB and SB, sync checked
expsSDG  = [28:32 56];                                            % SA recordings with usable gratings (SA5, SA9 ins 2)

%% MB vs SB (4 animals): statistics and leave-one-animal-out
stims = {'MB', 'RG'};
AllExpAnalysis(expsMBSB, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims));
leaveOneAnimalOut(d.structPath(expsMBSB, stims), figLabel='Figure 7 MB-SB', outName='LeaveOneAnimalOut_SA_MB-RG');

%% Gratings (2 animals, small n): movement and locality
for pr = {{'SDGm', 'SDGs'}, {'MB', 'SDGm'}, {'RG', 'SDGs'}}
    stims = pr{1};                                                % compared pair
    AllExpAnalysis(expsSDG, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims));   % hierarchical bootstrap over 2 animals
    leaveOneAnimalOut(d.structPath(expsSDG, stims), figLabel=sprintf('Figure 7 %s-%s', stims{:}), ...
        outName=sprintf('LeaveOneAnimalOut_SA_%s-%s', stims{:}));  % with 2 animals: each animal alone
end
