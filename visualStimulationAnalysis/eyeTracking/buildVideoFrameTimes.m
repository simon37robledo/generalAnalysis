function F = buildVideoFrameTimes(experimentNum, params)
% BUILDVIDEOFRAMETIMES  Ephys time of every saved eye-video frame, from the camera TTLs and the
% camera's own frame timestamps.
%
%   F = buildVideoFrameTimes(experimentNum)              % load the cache or build it
%   F = buildVideoFrameTimes(experimentNum, overwrite=true)
%
%   The camera is hardware-triggered by the TTLs on digital bit 5 (sync_events\out_5.txt, which
%   getTrigger returns as T{7}). Every trigger produces an exposure, but not every exposure is saved:
%     - frames dropped by the acquisition PC show up as gaps in the camera timestamp CSV;
%     - PV97 saved every second exposure (100 Hz trigger, 50 Hz video).
%   Saved frame i therefore sits on TTL  1 + ratio * (i - 1 + dropped frames before i)  of the saved
%   segment (the TTLs between the two >1 s pauses that bracket the video). This is exact when
%   saved-segment TTLs == ratio * (frames + dropped). A mismatch of up to maxClosureFrames frames
%   is tolerated with a warning: frame times then carry an unknown constant offset of at most that
%   many frames (F.timeUncertaintyMs). A larger mismatch is an error.
%
%   Replaces videoTimeStampsSynced.mat (video_analysis.m, NeuroPix_Vstim_Analysis): that file removed
%   TTLs at the wrong indices for dropped frames, so frames slipped by up to ~2 s after the first drop
%   (checked 2026-11-19). The old file is left untouched; this writes <catgt>\videoFrameTimes.mat.
%
%   Output F (struct): frameTimeMs [nFrames x 1] (NaN = no TTL, e.g. video ran past the end of ephys),
%   counts, closure and residual checks, syncVersion.

arguments
    experimentNum (1,1) double                                                  % Excel row of the insertion
    params.manifest string = "W:\Large_scale_mapping_NP\EyeTracking\manifest.csv"   % file paths per insertion
    params.overwrite logical = false                                            % rebuild even if cached
    params.pauseMs double = 1000                                                % TTL pause that brackets the saved video
    params.maxClosureFrames double = 2                                          % tolerated TTL/frame count mismatch, in frames (PV97 49, 53: 2)
end

syncVersion = 2;                                                                % bump when the mapping or its outputs change (2: closure tolerance, uncertainty field)

M = readtable(params.manifest, 'TextType', 'string');                           % manifest (makeEyeTrackingManifest)
r = M(M.experimentNum == experimentNum, :);                                     % this insertion
assert(height(r) == 1, 'buildVideoFrameTimes:row', 'experimentNum %d: %d manifest rows (need 1).', experimentNum, height(r));
outF = fullfile(r.rec_dir, "videoFrameTimes.mat");                              % cache next to the old sync file

if isfile(outF) && ~params.overwrite                                            % reuse the cache
    S = load(outF, 'F');                                                        % cached result
    if isfield(S.F, 'syncVersion') && S.F.syncVersion == syncVersion            % current version
        F = S.F; return                                                         % done
    end
    warning('buildVideoFrameTimes:staleCache', 'experimentNum %d: cache version differs, rebuilding.', experimentNum);
end

%% Camera TTLs, read exactly as getTrigger does
sd = fullfile(r.rec_dir, "sync_events");                                        % sync event folder of the catgt run
names = string({dir(sd).name});                                                 % its files, in dir() order (as getTrigger)
names = names(contains(names, 'out'));                                          % the files getTrigger reads
assert(numel(names) >= 7 && names(7) == "out_5.txt", 'buildVideoFrameTimes:bit', ...
    'experimentNum %d: getTrigger T{7} is %s, not out_5.txt (camera bit 5).', experimentNum, names(min(7, end)));   % camera trigger file
ttl = 1000 * importdata(fullfile(sd, "out_5.txt"));                             % rising edges of the camera trigger [ms], as getTrigger
ttl = ttl(:);                                                                   % column
pauses = find(diff(ttl) > params.pauseMs);                                      % pauses in the trigger train
assert(~isempty(pauses) && numel(pauses) <= 2, 'buildVideoFrameTimes:pauses', ...
    'experimentNum %d: %d trigger pauses > %g ms (expected 1 or 2).', experimentNum, numel(pauses), params.pauseMs);   % unknown layout
if numel(pauses) == 2                                                           % video started and stopped inside ephys
    seg = ttl(pauses(1)+1 : pauses(2));                                         % TTLs of the saved video
    endsWithEphys = false;                                                      % video stopped first
else                                                                            % ephys stopped while the video was still saving
    seg = ttl(pauses(1)+1 : end);                                               % TTLs up to the end of ephys
    endsWithEphys = true;                                                       % trailing frames have no TTL
end

%% Camera frame timestamps (PC clock, one row per saved frame)
ts = readmatrix(r.camera_ts_csv);                                               % timestamps [s]
ts = (ts(:) - ts(1)) * 1000;                                                    % elapsed [ms]
nF = numel(ts);                                                                 % saved frames
assert(all(~isnan(ts)), 'buildVideoFrameTimes:tsNaN', 'experimentNum %d: NaN camera timestamps.', experimentNum);   % unreadable rows
assert(nF == r.n_dlc_frames, 'buildVideoFrameTimes:dlcRows', ...
    'experimentNum %d: %d camera timestamps but %d DLC rows.', experimentNum, nF, r.n_dlc_frames);   % DLC analysed every saved frame
V = VideoReader(r.dlc_video);                                                   % the video DLC ran on
assert(V.NumFrames == nF, 'buildVideoFrameTimes:mp4Frames', ...
    'experimentNum %d: %d mp4 frames but %d camera timestamps.', experimentNum, V.NumFrames, nF);   % video and timestamps agree

%% TTLs per saved frame and dropped frames
framePeriod = median(diff(ts));                                                 % saved-frame period [ms]
ttlPeriod = median(diff(seg));                                                  % trigger period [ms]
ratio = round(framePeriod / ttlPeriod);                                         % TTLs per saved frame (1, or 2 for PV97)
assert(ismember(ratio, [1 2]) && abs(framePeriod / ttlPeriod - ratio) < 0.01 * ratio, 'buildVideoFrameTimes:ratio', ...
    'experimentNum %d: frame/TTL period ratio %.4f is not 1 or 2.', experimentNum, framePeriod / ttlPeriod);   % unexpected rates
assert(all(diff(seg) < 1.5 * ttlPeriod), 'buildVideoFrameTimes:ttlGap', ...
    'experimentNum %d: missing camera TTLs inside the saved segment.', experimentNum);   % TTLs must be regular for the mapping
step = diff(ts) / framePeriod;                                                  % intervals in frame periods
assert(all(step > 0.5), 'buildVideoFrameTimes:tsOrder', ...
    'experimentNum %d: %d camera timestamp intervals shorter than half a frame.', experimentNum, nnz(step <= 0.5));   % duplicated/out-of-order frames
dropped = round(step) - 1;                                                      % frames dropped before each next saved frame
nDropped = sum(dropped);                                                        % total dropped frames

%% Closure: saved TTLs == ratio x (saved + dropped frames)
closure = numel(seg) - ratio * (nF + nDropped);                                 % 0 = exact
if endsWithEphys                                                                % video outlasts ephys: fewer TTLs is expected
    assert(closure <= ratio, 'buildVideoFrameTimes:closureEnd', ...
        'experimentNum %d: %d more TTLs than frames although ephys ended first.', experimentNum, closure);   % impossible layout
elseif closure ~= 0                                                             % video bracketed by pauses: must close
    assert(abs(closure) <= ratio * params.maxClosureFrames, 'buildVideoFrameTimes:closure', ...
        'experimentNum %d: saved TTLs %d vs %d x (frames %d + dropped %d): off by %d TTLs.', ...
        experimentNum, numel(seg), ratio, nF, nDropped, closure);               % more than maxClosureFrames: mapping unknown
    warning('buildVideoFrameTimes:closureOff', ...
        'experimentNum %d: closure off by %d TTL(s); frame times uncertain by up to %d frame(s).', ...
        experimentNum, closure, ceil(abs(closure) / ratio));                    % tolerated, reported
end
uncertaintyMs = (~endsWithEphys) * ceil(abs(closure) / ratio) * framePeriod;    % worst-case constant offset of all frame times [ms]

%% Frame -> TTL
k = 1 + ratio * ((0:nF-1)' + [0; cumsum(dropped)]);                             % TTL index of each saved frame
has = k <= numel(seg);                                                          % frames that have a TTL
frameTimeMs = nan(nF, 1);                                                       % NaN = no ephys time
frameTimeMs(has) = seg(k(has));                                                 % ephys time of each frame [ms]
if ~endsWithEphys                                                               % every frame must have a TTL
    assert(all(has) || closure < 0, 'buildVideoFrameTimes:unmapped', ...
        'experimentNum %d: %d frames without a TTL.', experimentNum, nnz(~has));   % only allowed by a -1 closure
end

%% Residual check: TTL time vs camera clock (slow drift removed)
d = frameTimeMs(has) - ts(has);                                                 % ephys minus camera clock [ms]
res = d - movmedian(d, 2001);                                                   % remove clock offset and drift (~20-40 s window)
maxRes = max(abs(res));                                                         % worst frame
assert(maxRes < 0.5 * framePeriod, 'buildVideoFrameTimes:residual', ...
    'experimentNum %d: frame-time residual %.1f ms >= half a frame (%.1f ms): a frame slipped.', ...
    experimentNum, maxRes, 0.5 * framePeriod);                                  % a slip moves the residual by a whole frame

%% Compare with the old mapping (information only)
oldMax = NaN;                                                                   % old file missing
oldF = fullfile(r.rec_dir, "videoTimeStampsSynced.mat");                        % video_analysis.m output
if isfile(oldF)
    old = load(oldF, 'vidTTLsR'); old = old.vidTTLsR(:);                        % old frame times [ms]
    m = min(numel(old), nF);                                                    % row-matched part
    oldMax = max(abs(old(1:m) - frameTimeMs(1:m)), [], 'omitnan');              % largest disagreement [ms]
end

%% Save
F = struct('experimentNum', experimentNum, 'animal', r.animal, 'insertion', r.insertion, ...
    'frameTimeMs', frameTimeMs, 'nFrames', nF, 'nDropped', nDropped, 'nSegmentTTL', numel(seg), ...
    'ttlPerFrame', ratio, 'framePeriodMs', framePeriod, 'closure', closure, 'endsWithEphys', endsWithEphys, ...
    'timeUncertaintyMs', uncertaintyMs, ...
    'nFramesNoTTL', nnz(~has), 'maxResidualMs', maxRes, 'sdResidualMs', std(res), ...
    'oldMappingMaxDiffMs', oldMax, 'syncVersion', syncVersion, 'created', datetime('now'));   % result and provenance
save(outF, 'F');                                                                % cache
fprintf('experimentNum %d (%s ins %d): %d frames, %d dropped, %d TTL/frame, closure %d, %d frames without TTL, residual max %.1f ms; old mapping off by up to %.1f ms.\n', ...
    experimentNum, r.animal, r.insertion, nF, nDropped, ratio, closure, nnz(~has), maxRes, oldMax);   % summary line
end
