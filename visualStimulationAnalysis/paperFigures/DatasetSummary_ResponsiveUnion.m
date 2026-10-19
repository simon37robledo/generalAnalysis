%% DATASET SUMMARY — units responsive to at least one stimulus
% Cross-stimulus, so it belongs to no single figure. Reads each family's cached
% statistics with a bare call (never recomputes): run Fig1, Fig3 and Fig4 first
% so every family's cache carries the figureDesign statistics. A unit counts as
% tested for a family only if present throughout that family's block.

d = figureDesign();                                               % shared recording list
paperFigPrefix('summary');                                           % saved panels are named summary_<panel>
[Tresp, Sresp] = AllExpResponsiveUnion(d.exps, ...
    Stimuli        = {'MB','RG','MBR','SDGm','SDGs','FFF'}, ...   % families in the paper
    threshold      = 0.05, ...                                    % same cutoff as the figures
    useTtest       = false, ...                                   % permutation p
    computeMissing = false, ...                                   % cache-only: never writes stats
    overwrite      = true, ...                                    % rebuild the pooled union table
    plotFig        = true, ...
    PaperFig       = true);

fprintf('%d / %d testable units (%d recorded) respond to >=1 stimulus\n', ...
    Sresp.totalUnion, Sresp.totalTestable, Sresp.totalSomatic);
