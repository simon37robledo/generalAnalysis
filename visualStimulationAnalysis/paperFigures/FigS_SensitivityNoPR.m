%% SUPPLEMENTARY — sensitivity to BombCell's whole-recording presence criterion (Figures 1 and 7)
% BombCell labels a unit 'good' only if it fires in most of the WHOLE recording (minPresenceRatio).
% Some recordings are hours long while the stimuli occupy their first part, so units present during
% the stimuli but lost later are labelled MUA. Here every recording (all animals, both species) is
% relabelled from its own saved BombCell metrics and parameters with only minPresenceRatio = 0
% (scratch relabelNoPR.m, 2026-11-05: the unchanged parameters reproduced every current label;
% good units 1063 -> 2602). Units still count only if present throughout every compared block
% (per-block presence >= 0.7), as in the main analysis. Everything else is the main design.
% Labels: <sorting parent>\noPR_<sorting folder>; caches: <stim>_Analysis_noPR; outputs:
% Paper_figs\Sensitivity_noPR (analysisVariant). Main labels, caches and figures are not touched.

analysisVariant('noPR');                                          % all objects built below use the noPR labels/caches
try
    d = figureDesign();                                           % shared design; d.paperFigs -> Paper_figs\Sensitivity_noPR
    paperFigPrefix('figS_noPR');                                  % saved panels are named figS_noPR_<panel>

    %% Figure 1 (Pogona): MB vs SB with leave-one-animal-out
    stims = {'MB', 'RG'};
    AllExpAnalysis(d.exps, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims));
    leaveOneAnimalOut(d.structPath(d.exps, stims), figLabel='Sensitivity noPR Figure 1 MB-SB', ...
        outName='LeaveOneAnimalOut_noPR_Fig1');

    %% Figure 7 (Laudakia): same recordings as Fig7_SA.m
    expsMBSB = [28:32 44 45 47 48 56 98];                         % SA5, SA8, SA9 ins 2, SA11
    expsSDG  = [28:32 56];                                        % SA recordings with usable gratings
    AllExpAnalysis(expsMBSB, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims));
    leaveOneAnimalOut(d.structPath(expsMBSB, stims), figLabel='Sensitivity noPR Figure 7 MB-SB', ...
        outName='LeaveOneAnimalOut_noPR_SA_MB-RG');
    for pr = {{'SDGm', 'SDGs'}, {'MB', 'SDGm'}, {'RG', 'SDGs'}}
        stims = pr{1};                                            % compared pair
        AllExpAnalysis(expsSDG, d.statsArgs{:}, 'ComparePairs', stims, 'Lock', d.lock(stims));
        leaveOneAnimalOut(d.structPath(expsSDG, stims), figLabel=sprintf('Sensitivity noPR Figure 7 %s-%s', stims{:}), ...
            outName=sprintf('LeaveOneAnimalOut_noPR_SA_%s-%s', stims{:}));   % with 2 animals: each animal alone
    end
catch ME
    analysisVariant('');                                          % never leave the session in the variant
    rethrow(ME);
end
analysisVariant('');                                              % back to the main analysis
