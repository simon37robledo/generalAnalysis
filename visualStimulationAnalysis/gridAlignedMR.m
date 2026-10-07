function [MR, baseRates, info] = gridAlignedMR(Stats, u, spikeCol, tMs, preBase, postLen)
% gridAlignedMR  Trial spike matrix of one unit, aligned the way the MB spatial-grid test sees it.
%
%   [MR, baseRates, info] = gridAlignedMR(Stats, u, spikeCol, tMs, preBase, postLen)
%
%   The spatial-grid test compares the mean rate in a window starting when the
%   ball crosses a grid cell with the trial's pre-stimulus baseline. To show
%   exactly that comparison, each returned trial row is the concatenation of
%     columns 1 : preBase             the pre-trial baseline [onset - preBase, onset)
%     columns preBase+1 : preBase+postLen   [crossing, crossing + postLen)
%   for the trials of the unit's preferred direction that cross its preferred
%   grid cell, at the speed with the unit's lowest p-value (as in AllExpAnalysis).
%   Column preBase+1 is therefore the crossing, and everything left of it is
%   true baseline, so a raster/PSTH or the FP/FN screen built on MR uses the
%   same baseline and response as the test.
%
%   INPUTS
%     Stats     grid result of StatisticsPerNeuron (SpatialGridMode=true,
%               statsVersion >= 5): Speed* fields with pvalsResponse,
%               prefDirection, prefGridCell, trialOnsets, trialDirections, crossingMs
%     u         unit index (column of goodU)
%     spikeCol  goodU(:, u) of the spike-sorting ic matrix
%     tMs       spike times p_sort.t (ms)
%     preBase   baseline length (ms)
%     postLen   post-crossing length (ms)
%
%   OUTPUTS
%     MR        [nSel x (preBase + postLen)] spike counts per 1 ms bin
%     baseRates [nAll x 1] per-trial baseline rate (spikes/s) of ALL trials of
%               that speed (stimulus-independent; for z-scoring)
%     info      struct: speedField, direction, gridCell, nTrials

spF = fieldnames(Stats);                                           % all fields of the stats struct
spF = spF(startsWith(string(spF), "Speed"));                       % Speed1, Speed2, ...
assert(~isempty(spF), 'gridAlignedMR: stats have no Speed* fields (not a moving-ball result).');
need = {'pvalsResponse','prefDirection','prefGridCell','trialOnsets','trialDirections','crossingMs'};
for k = 1:numel(spF)
    assert(all(isfield(Stats.(spF{k}), need)), ...
        ['gridAlignedMR: %s lacks grid alignment fields (needs SpatialGridMode stats with ', ...
         'statsVersion >= 5); rerun AllExpAnalysis with overwriteStats=true.'], spF{k});  % fail loudly
end

pU = cellfun(@(f) Stats.(f).pvalsResponse(u), spF);                % this unit's p-value per speed
[~, k] = min(pU);                                                  % best speed (lowest p), as in AllExpAnalysis
sub = Stats.(spF{k});                                              % that speed's grid result

dirs = unique(sub.trialDirections);                                % direction levels, same order as the grid test
dU   = dirs(sub.prefDirection(u));                                 % preferred direction value
c    = sub.prefGridCell(u);                                        % preferred grid cell index
rows = sub.trialDirections(:) == dU & ~isnan(sub.crossingMs(:, c));  % trials of that direction crossing that cell
assert(any(rows), 'gridAlignedMR: unit %d has no trials crossing its preferred cell.', u);  % invariant of the test

ons  = sub.trialOnsets(rows);                                      % trial onsets (ms)
cros = ons(:) + sub.crossingMs(rows, c);                           % crossing times (ms)
nSel = numel(ons);                                                 % selected trials

Mb = BuildBurstMatrix(spikeCol, round(tMs), round(ons(:)' - preBase), round(preBase));  % pre-trial baseline
Mr = BuildBurstMatrix(spikeCol, round(tMs), round(cros(:)'), round(postLen));            % from the crossing
MR = [reshape(Mb, nSel, []), reshape(Mr, nSel, [])];               % [nSel x (preBase + postLen)]

allOns = sub.trialOnsets(:)';                                      % every trial of this speed
Ball   = BuildBurstMatrix(spikeCol, round(tMs), round(allOns - preBase), round(preBase));  % all baselines
baseRates = mean(reshape(Ball, numel(allOns), []), 2) * 1000;      % per-trial baseline rate (spikes/s)

info = struct('speedField', spF{k}, 'direction', dU, 'gridCell', c, 'nTrials', nSel);  % what was selected
end
