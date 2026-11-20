function [Tr, info] = eyeTrialMetrics(vsObj, E, calib, params)
% EYETRIALMETRICS  Eye movement in every trial of one stimulus block.
%
%   [Tr, info] = eyeTrialMetrics(vsObj, E, calib)
%     vsObj   stimulus object (linearlyMovingBallAnalysis, rectGridAnalysis, ...); trial times from
%             vsObj.getSyncedDiodeTriggers (same clock as the camera TTLs: sync_events of the catgt run)
%     E, calib  from getEyePosition
%
%   Windows per trial: stimulus (onset -> offset), the first 300 ms after onset (same length for
%   every stimulus), and the 300 ms before onset (spontaneous reference).
%   Per window: fraction of valid frames, maximum displacement from the window's starting position,
%   saccade-like events. Per stimulus window also: mean position offset from the block median and
%   the net velocity (end - start position / duration).
%   Degrees = asind(px / R) with R = calib.RPx; pupil-diameter units = px / calib.pupilDiamPx.
%
%   Saccade-like event: >= minEventFrames consecutive frames with speed above
%   max(noise median + 5 MAD, minSpeedDegS) and a net amplitude >= minEventDeg.
%
%   Smoothing (user, 2026-11-19): camera vibration (49.4 Hz in rec 65, sampled at 99.5 Hz, so it
%   alternates frame to frame: lag-1 autocorrelation of dx -0.80; jitter up to 0.46 deg/frame)
%   would otherwise dominate the maximum displacement. In 100 Hz videos x and y first pass a
%   2-frame mean (gain ~0.02 at 49.5 Hz; a median alone keeps alternating values), then a
%   centred 50 ms moving median (smoothMs), which keeps saccade steps. At 50 Hz (PV97) 49 Hz
%   vibration aliases to ~1 Hz and is not removed; PV97 jitter is small (0.025 deg/frame).
%   Raw and smoothed jitter are reported (info.rawNoisePx, info.smoothNoisePx).

arguments
    vsObj                                                                          % stimulus object
    E table                                                                        % eye position (getEyePosition)
    calib struct                                                                   % calibration (getEyePosition)
    params.windowMs double = 300                                                   % fixed post-onset and pre-onset window [ms]
    params.edgeMs double = 50                                                      % start/end position = median of the first/last edgeMs
    params.minValidFrac double = 0.8                                               % trial window kept if >= this fraction of frames is valid
    params.minEventFrames double = 2                                               % saccade-like event: consecutive fast frames
    params.minSpeedDegS double = 20                                                % saccade-like event: speed floor [deg/s]
    params.minEventDeg double = 0.5                                                % saccade-like event: net amplitude [deg]
    params.smoothMs double = 50                                                    % moving-median window [ms] applied to x and y first (0 = raw)
    params.matchMs double = NaN                                                    % extra window 0 -> matchMs after onset (e.g. the SB duration, for MB); NaN = none
end

S = vsObj.getSyncedDiodeTriggers;                                                  % trial onsets/offsets [ms]
on = S.stimOnFlipTimes(:); off = S.stimOffFlipTimes(:);                            % columns
assert(numel(on) == numel(off) && all(off > on), 'eyeTrialMetrics:trials', ...
    '%s: %d onsets, %d offsets, or an offset before its onset.', class(vsObj), numel(on), numel(off));   % paired trials
assert(numel(on) == vsObj.VST.nTotTrials, 'eyeTrialMetrics:nTrials', ...
    '%s: %d synced onsets but %d trials in VST.', class(vsObj), numel(on), vsObj.VST.nTotTrials);   % trial order = VST order

t = E.tMs; fin = isfinite(t);                                                      % frame times and frames that have one
assert(all(fin(1:find(fin, 1, 'last'))), 'eyeTrialMetrics:timeGaps', 'Frames without a time inside the video.');   % NaNs only at the end
t = t(fin); x = E.xPx(fin); y = E.yPx(fin); v = E.valid(fin);                      % frames with a time
assert(all(diff(t) > 0), 'eyeTrialMetrics:timeOrder', 'Frame times are not strictly increasing.');   % needed for the index lookup
fp = calib.framePeriodMs;                                                          % frame period [ms]
R = calib.RPx;                                                                     % eye radius [px] (NaN = no degrees)
rawStep = hypot(diff(x), diff(y));                                                 % frame-to-frame change before smoothing [px]
rawPairOK = v(1:end-1) & v(2:end) & diff(t) < 1.5 * fp;                            % consecutive valid frames
rawNoisePx = median(rawStep(rawPairOK));                                           % raw jitter (camera vibration + fit noise) [px/frame]
if params.smoothMs > 0                                                             % remove camera vibration before every metric
    x(~v) = NaN; y(~v) = NaN;                                                      % invalid frames do not enter the filters
    if 1000 / fp > 60                                                              % ~100 Hz video: vibration (~49.5 Hz) alternates frame to frame
        x = movmean(x, 2, 'omitnan'); y = movmean(y, 2, 'omitnan');                % 2-frame mean: gain ~0.02 at 49.5 Hz (a median keeps alternation)
    end
    x = movmedian(x, params.smoothMs, 'omitnan', 'SamplePoints', t);               % centred moving median over smoothMs of ephys time
    y = movmedian(y, params.smoothMs, 'omitnan', 'SamplePoints', t);               % (window in ms, so dropped frames do not widen it)
end
toDeg = @(px) asind(min(px ./ R, 1)) + 0 * R;                                      % px shift -> rotation [deg]; NaN when R is NaN
                                                                                   % (+0*R: min(NaN, 1) returns 1 in MATLAB, which gave 90 deg; fixed 2026-11-19)
nOutside = nnz(on < t(1) | off > t(end));                                          % trials outside the video
if nOutside > 0                                                                    % video started late / stopped early
    warning('eyeTrialMetrics:outsideVideo', '%s: %d of %d trials lie (partly) outside the video.', class(vsObj), nOutside, numel(on));
end

% Saccade threshold from the frame-to-frame speed of the whole recording
step = hypot(diff(x), diff(y));                                                    % px per frame
pairOK = v(1:end-1) & v(2:end) & diff(t) < 1.5 * fp;                               % consecutive valid frames
noiseThr = median(step(pairOK)) + 5 * 1.4826 * mad(step(pairOK), 1);               % noise threshold [px/frame]
floorThr = params.minSpeedDegS * fp / 1000;                                        % speed floor [deg/frame]
floorPx = R * sind(floorThr);                                                      % ... in px/frame
fastThr = max(noiseThr, floorPx);                                                  % event speed threshold [px/frame] (noise only when R is NaN)
fast = [pairOK & step > fastThr; false];                                           % frame k fast = step k->k+1 above threshold

blockOn = on(1) - params.windowMs; blockOff = off(end);                            % block span incl. the first baseline
inBlock = t >= blockOn & t <= blockOff & v;                                        % valid frames of the block
xMed = median(x(inBlock)); yMed = median(y(inBlock));                              % block median position [px]

nT = numel(on);                                                                    % trials
Tr = table((1:nT)', on, off, off - on, 'VariableNames', {'trial', 'onMs', 'offMs', 'durMs'});   % output
stimW = windowStats(on, off, true);                                                % whole stimulus window
postW = windowStats(on, on + params.windowMs, false);                              % fixed 300 ms after onset
baseW = windowStats(on - params.windowMs, on, false);                              % 300 ms before onset
Tr = [Tr, prefix(stimW, 'stim_'), prefix(postW, 'post_'), prefix(baseW, 'base_')]; % one row per trial
if isfinite(params.matchMs)                                                        % duration-matched window (MB vs SB comparison)
    Tr = [Tr, prefix(windowStats(on, on + params.matchMs, false), 'match_')];      % 0 -> matchMs after onset
end

if isprop(vsObj, 'VST') && isfield(vsObj.VST, 'directions')                        % moving-ball trials: ball velocity
    Tr.direction = vsObj.VST.directions(:);                                        % ball direction [rad]
    Tr.speedPxS = vsObj.VST.speeds(:);                                             % ball speed [px/s on screen]
end

info = struct('noiseThrPx', noiseThr, 'fastThrPx', fastThr, 'xMedPx', xMed, 'yMedPx', yMed, ...
    'rawNoisePx', rawNoisePx, 'smoothNoisePx', median(step(pairOK)), 'smoothMs', params.smoothMs, ...   % jitter before / after smoothing
    'nTrials', nT, 'nKept', nnz(Tr.stim_kept), 'blockStartMs', blockOn, 'blockEndMs', blockOff);   % block summary

    function W = windowStats(a, b, full)
        % Metrics of each window [a(i), b(i)) for all trials.
        n = numel(a);                                                              % trials
        W = table();                                                               % output
        W.fracValid = nan(n, 1); W.kept = false(n, 1);                             % coverage
        W.dispPx = nan(n, 1); W.nEvents = nan(n, 1);                               % displacement, events
        if full, W.meanXPx = nan(n, 1); W.meanYPx = nan(n, 1); W.velXPxS = nan(n, 1); W.velYPxS = nan(n, 1); end
        for i = 1:n
            k = find(t >= a(i) & t < b(i));                                        % frames in the window
            nExp = (b(i) - a(i)) / fp;                                             % frames expected
            if isempty(k) || numel(k) < 0.9 * nExp, continue, end                  % window not covered by the video
            W.fracValid(i) = mean(v(k));                                           % valid fraction
            if W.fracValid(i) < params.minValidFrac, continue, end                 % too few valid frames
            W.kept(i) = true;                                                      % usable window
            kv = k(v(k));                                                          % valid frames
            s0 = kv(t(kv) < t(kv(1)) + params.edgeMs);                             % first edgeMs
            s1 = kv(t(kv) > t(kv(end)) - params.edgeMs);                           % last edgeMs
            x0 = median(x(s0)); y0 = median(y(s0));                                % start position
            W.dispPx(i) = max(hypot(x(kv) - x0, y(kv) - y0));                      % max displacement from start [px]
            W.nEvents(i) = countEvents(k);                                         % saccade-like events
            if full
                W.meanXPx(i) = mean(x(kv)) - xMed; W.meanYPx(i) = mean(y(kv)) - yMed;   % offset from the block median [px]
                dur = (median(t(s1)) - median(t(s0))) / 1000;                      % time between start and end [s]
                W.velXPxS(i) = (median(x(s1)) - x0) / dur;                         % net velocity [px/s]
                W.velYPxS(i) = (median(y(s1)) - y0) / dur;
            end
        end
        W.dispDeg = toDeg(W.dispPx);                                               % [deg]
        W.dispPupil = W.dispPx / calib.pupilDiamPx;                                % [pupil diameters]
    end

    function nE = countEvents(k)
        % Runs of >= minEventFrames fast frames with net amplitude >= minEventDeg.
        f = fast(k);                                                               % fast flags in the window
        d = diff([false; f; false]);                                               % run edges
        st = find(d == 1); en = find(d == -1) - 1;                                 % run starts / ends (window-relative)
        nE = 0;                                                                    % count
        for j = 1:numel(st)
            if en(j) - st(j) + 1 < params.minEventFrames, continue, end            % too short
            i0 = k(st(j)); i1 = min(k(en(j)) + 1, numel(x));                       % frames before and after the run
            amp = hypot(x(i1) - x(i0), y(i1) - y(i0));                             % net amplitude [px]
            if toDeg(amp) >= params.minEventDeg, nE = nE + 1; end                  % large enough
        end
        if isnan(R), nE = NaN; end                                                 % no calibration: no amplitude criterion, so no event count
    end
end

function W = prefix(W, p)
% Prefix every column name.
W.Properties.VariableNames = strcat(p, W.Properties.VariableNames);                % e.g. stim_dispDeg
end
