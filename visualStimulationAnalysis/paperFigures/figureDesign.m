function d = figureDesign()
% figureDesign  Shared statistics design for every paper figure.
%
%   d = figureDesign()
%
%   Each stimulus has ONE StatisticsPerNeuron cache per recording, so every
%   figure must request identical statistics for a given stimulus; otherwise a
%   later figure silently recomputes (and overwrites) the cache an earlier
%   figure was built from. All figure scripts take their AllExpAnalysis
%   arguments from here.
%
%   Design (validated against a null simulation, see checkMovingWindowNullFPR
%   and checkLayoutNullFPR):
%     - 300 ms response windows and 300 ms pre-trial baseline (SpatialGridMode
%       forces BaseRespWindow = 300 for every stimulus in the call).
%     - MB: spatial-grid test (9 x 9 cells, window from each ball crossing).
%     - RG: ON window plus OFF window from each recording's offset.
%     - MBR: 9 equal segments across the sweep, matched-length baseline.
%     - Others (SDGs, SDGm, FFF, NI): single fixed window from onset;
%       SDGm against the blank pre-trial baseline.
%     - Responsiveness = best condition (maxCategory), permutation p < 0.05,
%       no FDR; units present throughout every compared block only.
%     - Luminosity locked to 255 for MB and RG (the only stimuli with levels).
%
%   Fields
%     exps        recording list used by every figure
%     paperFigs   figure / struct output folder
%     statsArgs   name-value cell for AllExpAnalysis (append ComparePairs, Lock)
%     lock        @(stims) per-stimulus Lock cell for AllExpAnalysis
%     structPath  @(exps, stims) path of the *_analysisStruct.mat AllExpAnalysis writes
%     rasterArgs  name-value cell shared by the plotRaster_MultiExp calls
%     psthArgs    name-value cell shared by the plotPSTH_MultiExp calls

d.exps      = [1, 8, 40:43, 49:54, 64:66, 68:85, 87:97];          % all recordings (stimuli absent from one are skipped)
d.paperFigs = 'W:\Large_scale_mapping_NP\Paper_figs';             % where AllExpAnalysis saves figures and structs

d.statsArgs = { ...
    'overwrite',          true, ...                               % rebuild the pooled tables
    'PaperFig',           true, ...                               % save figures and the analysis struct
    'overwriteResponse',  false, ...                              % response windows are design-independent
    'overwriteStats',     true, ...                               % recompute per-recording statistics
    'useFDR',             false, ...                              % p < 0.05 per unit, no FDR
    'SpatialGridMode',    true, ...                               % MB grid test; forces 300 ms windows
    'maxCategory',        true, ...                               % best condition
    'useTtest',           false, ...                              % permutation p-values (primary criterion)
    'offWindowStims',     {'RG'}, ...                             % ON + OFF windows
    'segmentStims',       {'MBR'}, ...                            % segmented test
    'nSegments',          9};                                     % 9 segments (as the 9 x 9 grid)

d.lock       = @lockFor;                                          % per-stimulus Lock cell
d.structPath = @(exps, stims) fullfile(d.paperFigs, ...
    sprintf('Ex_%d-%d_Combined_%s_analysisStruct.mat', exps(1), exps(end), strjoin(stims, '-')));  % AllExpAnalysis naming

d.rasterArgs = { ...
    'TakeTopPercentTrials', [], 'PaperFig', true, 'zScore', true, 'postStim', 500, ...
    'prefCatPSTH', false, 'statPrefCatPSTH', true, 'sdFloorPrctile', 25, 'speed', "best", ...
    'useCompleteWindow', false, 'MarkFPnFN', false, 'useTtest', false, 'sortBy', "meanZScore"};  % shared raster settings
d.psthArgs = { ...
    'zScore', true, 'TakeTopPercentTrials', [], 'PaperFig', true, 'byDepth', false, ...
    'binWidth', 50, 'postStim', 500};                             % shared PSTH settings
end

function L = lockFor(stims)
% Lock luminosity 255 for MB and RG (white stimuli); no lock for the others.
L = cell(1, numel(stims));                                        % one entry per stimulus
for k = 1:numel(stims)
    if ismember(stims{k}, {'MB', 'RG'})
        L{k} = {"luminosities", 255};                             % white ball / white squares only
    else
        L{k} = {};                                                % no luminosity levels to lock
    end
end
end
