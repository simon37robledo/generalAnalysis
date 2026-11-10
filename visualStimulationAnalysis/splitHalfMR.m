function [MR, baseRates, info] = splitHalfMR(sub, u, spikeCol, tMs, alignTo, preBase, postLen)
% splitHalfMR  Held-out trial spike matrices of one unit, split-half selection.
%
%   [MR, baseRates, info] = splitHalfMR(sub, u, spikeCol, tMs, alignTo, preBase, postLen)
%
%   Uses the split stored with DeltaRateCV (StatisticsPerNeuron cvInfo, or
%   cvSel/cvHalf for the MB spatial grid): the best condition chosen on half 1
%   of the repetitions is shown with the trials of half 2, and the condition
%   chosen on half 2 with the trials of half 1. Trials that chose a condition
%   are never the ones displayed for it, so a population average of these
%   PSTHs is free of the winner's-curse inflation of plotting the best
%   condition on its own trials. Average the two matrices' PSTHs (as
%   DeltaRateCV averages the two held-out estimates).
%
%   Each row is the pre-trial baseline [baseEnd - preBase, baseEnd) followed
%   by [event, event + postLen), so column preBase+1 is the alignment event.
%
%   INPUTS
%     sub       stats sub-struct of one stimulus (e.g. Stats.Speed1 for MB,
%               the RG struct, Stats.Moving for SDGm)
%     u         unit index (column of goodU)
%     spikeCol  goodU(:, u) of the spike-sorting ic matrix
%     tMs       spike times p_sort.t (ms)
%     alignTo   "onset"    trial onset (motion onset for SDGm = start of its test window)
%               "crossing" grid-tested MB only: closest approach to the selected cell
%               "offset"   onset-tested stimuli only: onset + stimulus duration
%     preBase   baseline length (ms)
%     postLen   length after the event (ms)
%
%   OUTPUTS
%     MR        {1 x 2} [nTrials_h x (preBase + postLen)] spike counts per 1 ms:
%               MR{h} = condition selected on half h, trials of the other half
%     baseRates [nAll x 1] per-trial baseline rate (spikes/s) of ALL trials
%               (stimulus-independent; for z-scoring)
%     info      struct: group and option (cell / segment) selected per half,
%               number of trials per matrix
%
%   See also: splitHalfDelta, eventAlignedMR, gridAlignedMR

alignTo = string(alignTo);                                         % accept char or string
assert(ismember(alignTo, ["onset","crossing","offset"]), ...
    'splitHalfMR: alignTo must be "onset", "crossing" or "offset".'); % valid alignment
isGrid = isfield(sub, 'prefGridCell');                             % MB spatial-grid result
MR = cell(1, 2);                                                   % held-out matrices, one per selection half
info = struct('group', [0 0], 'option', [0 0], 'nTrials', [0 0]);  % what was selected

if isGrid
    % ---- MB spatial grid: group = direction, option = grid cell ----
    assert(all(isfield(sub, {'cvSel','cvHalf','trialOnsets','trialDirections','crossingMs'})), ...
        'splitHalfMR:noCV', ['MB stats lack the stored split (cvSel/cvHalf): rerun ', ...
        'AllExpAnalysis with overwriteStats=true.']);              % fail loudly on an old cache
    assert(alignTo ~= "offset", 'splitHalfMR: "offset" is not defined for the moving ball.');  % no single offset per cell
    dirs  = unique(sub.trialDirections);                           % direction levels, same order as dirID in the test
    nGrp  = numel(dirs);                                           % groups of the split
    for h = 1:2
        sel = sub.cvSel(h, u);                                     % (direction, cell) chosen on half h
        g   = mod(sel - 1, nGrp) + 1;                              % direction index
        c   = ceil(sel / nGrp);                                    % grid cell index
        rows = sub.trialDirections(:) == dirs(g) & ~isnan(sub.crossingMs(:, c)) & ...
            sub.cvHalf(:) == 3 - h;                                % other half's trials crossing that cell in that direction
        ons  = sub.trialOnsets(rows);                              % their onsets (ms)
        if alignTo == "crossing"
            ev = ons(:) + sub.crossingMs(rows, c);                 % closest approach to the cell (ms)
        else
            ev = ons(:);                                           % trial onset
        end
        MR{h} = eventAlignedMR(spikeCol, tMs, ons, ev, preBase, postLen);  % baseline | from the event
        info.group(h) = g; info.option(h) = c; info.nTrials(h) = sum(rows);  % record the selection
    end
    baseEndAll = sub.trialOnsets(:);                               % baselines of every trial end at its onset
else
    % ---- Onset-tested stimuli: group = category, option = segment (ON/OFF, MBR segments) ----
    assert(isfield(sub, 'cvInfo') && ~isempty(sub.cvInfo.cvSel), 'splitHalfMR:noCV', ...
        ['Stats lack the stored split (cvInfo): rerun AllExpAnalysis with ', ...
        'overwriteStats=true.']);                                  % fail loudly on an old cache
    assert(alignTo ~= "crossing", 'splitHalfMR: "crossing" applies to the grid-tested moving ball only.');
    ci   = sub.cvInfo;                                             % stored split and trial timing
    nGrp = size(sub.validCats, 1);                                 % categories = groups of the split
    for h = 1:2
        sel = ci.cvSel(h, u);                                      % (category, segment) chosen on half h
        g   = mod(sel - 1, nGrp) + 1;                              % category index
        sg  = ceil(sel / nGrp);                                    % segment index (1 if not segmented)
        rows = ci.cvGrp(:) == g & ci.cvHalf(:) == 3 - h;           % other half's trials of that category
        ev = ci.cvAlignMs(rows);                                   % start of the test window (ms)
        if alignTo == "offset"
            ev = ev + ci.cvStimDur;                                % stimulus offset (ms)
        end
        MR{h} = eventAlignedMR(spikeCol, tMs, ci.cvBaseMs(rows), ev, preBase, postLen);  % baseline | from the event
        info.group(h) = g; info.option(h) = sg; info.nTrials(h) = sum(rows);  % record the selection
    end
    baseEndAll = ci.cvBaseMs(:);                                   % baselines of every trial
end

for h = 1:2
    assert(info.nTrials(h) > 0, 'splitHalfMR: unit %d has no held-out trials for selection %d.', u, h);  % invariant of the split
end
Ball = BuildBurstMatrix(spikeCol, round(tMs), round(baseEndAll' - preBase), round(preBase));  % all baselines
baseRates = mean(reshape(Ball, numel(baseEndAll), []), 2) * 1000;  % per-trial baseline rate (spikes/s)
end
