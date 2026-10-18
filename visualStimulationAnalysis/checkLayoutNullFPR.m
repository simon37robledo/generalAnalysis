function [fprTable, deltaTable, layouts] = checkLayoutNullFPR(params)
% checkLayoutNullFPR  Null calibration of the fixed-window and segmented tests on
%   the trial layouts of the other stimulus families (FFF, SDG, MBR, NI).
%
%   [fprTable, deltaTable, layouts] = checkLayoutNullFPR()
%
%   Under H0 (homogeneous firing, same rate before and during the stimulus) the
%   test only depends on the number of conditions, trials per condition and
%   windows, so each layout is read from a real recording's ResponseWindow and
%   simulated. The statistic mirrors StatisticsPerNeuron: per-trial window mean
%   minus the trial's baseline mean, maximum over conditions x windows of the
%   condition means, with the same maximum after the sign flip (segSignFlipNull).
%   A valid test gives a false-positive rate of ~alpha; the split-half delta rate
%   (splitHalfDelta) should average ~0.
%
%   Reads cached ResponseWindow only; writes nothing to disk.
%
%   See also: checkMovingWindowNullFPR, segSignFlipNull, splitHalfDelta

arguments
    params.ratesHz double = [0.5 1 2 5 10 20]   % simulated firing rates (Hz)
    params.nPerRate double = 240                % simulated neurons per rate
    params.nBoot double = 2000                  % sign-flip permutations
    params.alpha double = 0.05                  % responsiveness threshold
    params.winMs double = 300                   % fixed response window length (ms)
    params.nSegments double = 9                 % equal segments per stimulus for segmented layouts (MBR)
    params.baseMs double = 300                  % baseline length (ms)
    params.sdgWinMs double = [1000 1923]        % [static drifting] grating windows (ms), as figureDesign sdgWindows; [] = winMs
    params.chunk double = 40                    % neurons simulated at once (memory)
    params.seed double = 1                      % RNG seed
end

rng(params.seed);                                                  % reproducible simulation

% =========================================================================
% 1. Real layouts: conditions (nCats), trials per condition (tc), windows (nWin)
% =========================================================================
% {label, recording, constructor, sub-field ('' = flat), segmented?, NI lizard subset?}
spec = {'FFF fixed',        52, @fullFieldFlashAnalysis,        '',       false, false; ...
        'SDGs fixed',       52, @StaticDriftingGratingAnalysis, 'Static', false, false; ...
        'SDGm fixed',       52, @StaticDriftingGratingAnalysis, 'Moving', false, false; ...
        'MBR 9 segments',   70, @linearlyMovingBarAnalysis,     'Speed1', true,  false; ...
        'NI fixed (lizards)', 92, @imageAnalysis,               '',       false, true};
layouts = struct('name', {}, 'nCats', {}, 'tc', {}, 'nWin', {}, 'winLen', {}, 'baseLen', {});   % filled below
for k = 1:size(spec, 1)
    evalc('NP = loadNPclassFromTable(spec{k,2}); o = spec{k,3}(NP); rw = o.ResponseWindow;');  % cached, read-only
    if isempty(spec{k,4}), sub = rw; else, sub = rw.(spec{k,4}); end   % condition matrix holder
    if isfield(sub, 'C'), C = sub.C; else, C = rw.C; end            % SDG phases share the flat C
    stimDur = sub.stimDur;                                           % this phase's stimulus duration (ms)
    nCats   = size(sub.NeuronVals, 2);                               % conditions, as in StatisticsPerNeuron
    tc      = round(size(C, 1) / nCats);                             % trials per condition (balanced design)
    if spec{k,6}                                                     % NI keeps trials [1:15 61:75] (lizard images)
        nKeep = 30;                                                  % 2 images x 15 repetitions
        nCats = round(nKeep / tc);                                   % conditions after the subset
    end
    nWin = 1; winLen = params.winMs; baseLen = params.baseMs;        % fixed window vs full baseline
    if startsWith(spec{k,1}, 'SDG') && ~isempty(params.sdgWinMs)    % gratings: per-phase window (StatisticsPerNeuron SDGWindows)
        winLen  = min(params.sdgWinMs(1 + strcmp(spec{k,4}, 'Moving')), floor(stimDur));   % capped at the phase duration
        baseLen = min(winLen, rw.stimInter - 100);                   % same-length blank baseline, capped by the ITI
    end
    if spec{k,5}, nWin = params.nSegments; winLen = round(stimDur / nWin); baseLen = winLen; end  % equal segments, matched-length baseline (StatisticsPerNeuron)
    layouts(k) = struct('name', spec{k,1}, 'nCats', nCats, 'tc', tc, 'nWin', nWin, 'winLen', winLen, 'baseLen', baseLen);  % this family's layout
    fprintf('%-20s exp %d: %d conditions x %d trials, %d window(s) of %d ms\n', spec{k,1}, spec{k,2}, nCats, tc, nWin, winLen);
end

% =========================================================================
% 2. Simulate no-response neurons and run the test on each layout
% =========================================================================
fprTable = table(); deltaTable = table();                          % accumulators
for k = 1:numel(layouts)
    L   = layouts(k);
    nTr = L.nCats * L.tc;                                          % trials (contiguous condition blocks)
    grp  = repelem((1:L.nCats)', L.tc);                            % condition of each trial
    half = repmat(2 - mod((1:L.tc)', 2), L.nCats, 1);              % odd / even repetition
    for r = params.ratesHz
        lam = r / 1000;                                            % spike probability per 1 ms bin
        p = []; d = [];                                            % p-values and split-half deltas
        for c = 1:ceil(params.nPerRate / params.chunk)
            nU = min(params.chunk, params.nPerRate - (c-1)*params.chunk);  % neurons in this chunk
            base = mean(rand(nTr, nU, L.baseLen, 'single') < lam, 3);   % baseline mean rate per trial
            D = zeros(nTr, nU, L.nWin, 'single');                  % window mean minus baseline
            for w = 1:L.nWin
                D(:,:,w) = mean(rand(nTr, nU, L.winLen, 'single') < lam, 3) - base;  % independent window, shared baseline
            end
            p = [p, segSignFlipNull(double(D), L.tc, L.nCats, params.nBoot)];  %#ok<AGROW> max over conditions x windows
            d = [d, splitHalfDelta(D, grp, half)];                 %#ok<AGROW> cross-validated delta rate
        end
        nSig = sum(p < params.alpha);                              % false positives
        [phat, ci] = binofit(nSig, numel(p));                      % Clopper-Pearson 95% CI
        fprTable = [fprTable; table(string(L.name), r, numel(p), nSig, phat, ci(1), ci(2), ...
            'VariableNames', {'test','rateHz','n','nSig','FPR','ciLo','ciHi'})];   %#ok<AGROW>
        deltaTable = [deltaTable; table(string(L.name), r, numel(d), mean(d, 'omitnan'), ...
            std(d, 'omitnan') / sqrt(sum(~isnan(d))), ...
            'VariableNames', {'layout','rateHz','n','meanDeltaHz','seDeltaHz'})];  %#ok<AGROW>
    end
end
disp(fprTable); disp(deltaTable);                                  % print both summaries
end
