%% DATASET SUMMARY — units responsive to at least one stimulus
% Cross-stimulus, so it belongs to no single figure. Reads each family's cached
% statistics with a bare call (never recomputes): run Fig1, Fig3 and Fig4 first
% so every family's cache carries the figureDesign statistics. A unit counts as
% tested for a family only if present throughout that family's block.

d = figureDesign();                                               % shared recording list
paperFigPrefix('summary');                                           % saved panels are named summary_<panel>
[Tresp, Sresp] = AllExpResponsiveUnion(d.exps, ...
    Stimuli        = {'MB','RG','MBR','SDGm','SDGs','FFF'}, ...   % families in the paper
    excludeStim    = struct('SDGm', setdiff(d.exps, d.expsSDG), 'SDGs', setdiff(d.exps, d.expsSDG)), ...   % gratings of recordings 1, 8 not analysed
    threshold      = 0.05, ...                                    % same cutoff as the figures
    useTtest       = false, ...                                   % permutation p
    computeMissing = false, ...                                   % cache-only: never writes stats
    overwrite      = true, ...                                    % rebuild the pooled union table
    plotFig        = true, ...
    PaperFig       = true);

fprintf('%d / %d testable units (%d recorded) respond to >=1 stimulus\n', ...
    Sresp.totalUnion, Sresp.totalTestable, Sresp.totalSomatic);

%% Laudakia (Figure 7 recordings, SA11 included): the stimuli Figure 7 analyses (2026-11-19)
expsSA    = [28:32 44 45 47 48 56 98];                            % as Fig7_SA.m: SA5, SA8, SA9, SA11
expsSAsdg = [28:32 56];                                           % recordings with usable gratings (SA5, SA9 insertion 2)
[TrespSA, SrespSA] = AllExpResponsiveUnion(expsSA, ...
    Stimuli        = {'MB','RG','SDGm','SDGs'}, ...               % families tested in Figure 7
    excludeStim    = struct('SDGm', setdiff(expsSA, expsSAsdg), 'SDGs', setdiff(expsSA, expsSAsdg)), ...   % gratings only where Figure 7 uses them
    threshold      = 0.05, ...                                    % same cutoff as the figures
    useTtest       = false, ...                                   % permutation p
    computeMissing = false, ...                                   % cache-only: never writes stats
    overwrite      = true, ...                                    % rebuild the pooled union table
    plotFig        = true, ...
    PaperFig       = true);
fprintf('Laudakia: %d / %d testable units (%d recorded) respond to >=1 stimulus\n', ...
    SrespSA.totalUnion, SrespSA.totalTestable, SrespSA.totalSomatic);
