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
%     - Others (FFF, NI): single 300 ms window from onset.
%     - Gratings: SDGs 1000 ms from static onset, SDGm 1923 ms from motion
%       onset (one period of the slowest temporal frequency: 300 ms missed
%       the slow, sustained drifting responses), both vs a blank pre-trial
%       baseline of the same length capped by the ITI (900 ms; 400 ms in
%       recordings 40-43). The per-category round must use the same windows.
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
    'nSegments',          9, ...                                  % 9 segments (as the 9 x 9 grid)
    'sdgWindows',         [1000 1923]};                           % gratings [static drifting] (ms): drifting = one period
                                                                  % of the slowest TF (1/0.52 Hz); static = the longest
                                                                  % window inside every static phase (1000-1500 ms)

d.lock       = @lockFor;                                          % per-stimulus Lock cell
d.structPath = @(exps, stims) fullfile(d.paperFigs, ...
    sprintf('Ex_%d-%d_Combined_%s_analysisStruct.mat', exps(1), exps(end), strjoin(stims, '-')));  % AllExpAnalysis naming

d.rasterArgs = { ...
    'TakeTopPercentTrials', [], 'PaperFig', true, 'zScore', true, 'postStim', 500, ...
    'prefCatPSTH', false, 'statPrefCatPSTH', true, 'sdFloorPrctile', 25, 'speed', "best", ...
    'useCompleteWindow', false, 'MarkFPnFN', false, 'useTtest', false, 'sortBy', "meanZScore", ...
    'binWidth', 20};                                              % shared raster settings; 20 ms columns (10 ms rows were noisy)
d.psthArgs = { ...
    'zScore', true, 'TakeTopPercentTrials', [], 'PaperFig', true, 'byDepth', false, ...
    'binWidth', 50, 'postStim', 500};                             % shared PSTH settings
% PSTH bin rule: 50 ms bins for windows <= 1 s, 100 ms for 2 s windows (fewer
% recordings and slower responses: smaller bins are noise).
sdgW = d.statsArgs{find(strcmp(d.statsArgs, 'sdgWindows')) + 1}; % [static drifting] test windows (ms)
d.psthArgsSDG = setArg(setArg(d.psthArgs, 'postStim', 2000), 'binWidth', 100);   % drifting-grating comparisons (0-2 s)
d.psthArgsSDG = [d.psthArgsSDG, {'capMs', struct('SDGs', sdgW(1))}];   % SDGs drawn only over its static test window
d.psthArgsStatic = [setArg(d.psthArgs, 'postStim', 1000), {'capMs', struct('SDGs', sdgW(1))}];  % RG vs SDGs (0-1 s)
end

function a = setArg(a, name, value)
% Replace the value of one name-value pair in a cell of arguments.
k = find(strcmp(a(1:2:end), name), 1);                            % pair index
assert(~isempty(k), 'figureDesign: argument %s not found.', name);
a{2*k} = value;                                                   % new value
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
