function [E, calib] = getEyePosition(experimentNum, params)
% GETEYEPOSITION  Pupil position of one insertion on the ephys clock, with a per-recording
% pixel -> degree calibration.
%
%   [E, calib] = getEyePosition(experimentNum)
%
%   E (table, one row per saved video frame):
%     frame   video frame index (0-based, as in the pupil table)
%     tMs     ephys time of the frame [ms] (buildVideoFrameTimes; NaN = no camera TTL)
%     xPx, yPx  pupil-ellipse centre in camera pixels (Elena's table, no flips)
%     diamPx  pupil major diameter [px]
%     valid   fit_valid in the pupil table AND a finite ephys time
%
%   calib (struct): eye radius R [px] = half the distance between the DLC eye-corner points
%   (Rostral_edge, Caudal_edge); median pupil diameter; frame-to-frame noise floor.
%   Degrees: rotation = asind(shift_px / R). mm per px from the 3 mm retractor (Rostral_scale,
%   Caudal_scale) is stored for information only: it gave eye widths of 5.7-16.7 mm, also within
%   one animal, so it is not used.
%   Camera obliquity shortens shifts (underestimates rotation); the eye-opening width is at most
%   the globe diameter (R too small -> overestimates rotation). Report a +-20% R sensitivity.
%
%   Pupil table: fitPupilFromDLC output by default (EyeTracking\tables); Elena's pupil_metrics tables
%   give identical fits. Only the ellipse columns are read: her vid_ttls_r / ttl_to_frame columns copy
%   the old videoTimeStampsSynced.mat mapping, which slips after dropped frames, and are ignored.

arguments
    experimentNum (1,1) double                                                     % Excel row of the insertion
    params.manifest string = "W:\Large_scale_mapping_NP\EyeTracking\manifest.csv"  % file paths per insertion
    params.tableDir string = "W:\Large_scale_mapping_NP\EyeTracking\tables"        % pupil tables from fitPupilFromDLC (identical to Elena's
                                                                                   % pupil_metrics tables on 49, 57, 64, 65, 68, 69: 0.000 px, same valid flags)
    params.cacheDir string = "W:\Large_scale_mapping_NP\EyeTracking\cache"         % cached results
    params.overwrite logical = false                                               % recompute
    params.maxMissingRows double = 10                                              % table shorter than the video by more -> incomplete
    params.retractorMm double = 3                                                  % retractor length (user, 2026-11-19)
    params.maxLandmarkIQR double = 0.1                                             % landmark length accepted if IQR/median below this
    params.dlcSampleFrames double = 20000                                          % frames read from each of 3 places in the DLC CSV
end

eyeVersion = 1;                                                                    % bump when the output changes

M = readtable(params.manifest, 'TextType', 'string');                              % manifest
r = M(M.experimentNum == experimentNum, :);                                        % this insertion
assert(height(r) == 1, 'getEyePosition:row', 'experimentNum %d: %d manifest rows (need 1).', experimentNum, height(r));
tabF = fullfile(params.tableDir, r.out_table);                                     % Elena's table for this insertion
assert(isfile(tabF), 'getEyePosition:noTable', 'experimentNum %d: no pupil table %s.', experimentNum, tabF);   % not returned yet
tabInfo = dir(tabF);                                                               % table file date
tmp = dir(fullfile(params.tableDir, "." + extractBefore(r.out_table, ".csv") + "*"));   % hidden temp copy written while saving
if ~isempty(tmp) && max([tmp.datenum]) >= tabInfo.datenum                          % temp newer than the table: still being written
    error('getEyePosition:incompleteTable', 'experimentNum %d: %s is still being written (%s).', experimentNum, r.out_table, tmp(1).name);
elseif ~isempty(tmp)                                                               % older leftover of an earlier save
    warning('getEyePosition:staleTemp', 'experimentNum %d: ignoring older temp file %s.', experimentNum, tmp(1).name);
end

F = buildVideoFrameTimes(experimentNum);                                           % corrected frame times (cached)

cacheF = fullfile(params.cacheDir, sprintf('eyePos_%d.mat', experimentNum));       % cache file
if isfile(cacheF) && ~params.overwrite                                             % reuse
    S = load(cacheF, 'E', 'calib');                                                % cached
    if S.calib.eyeVersion == eyeVersion && S.calib.tableDatenum == tabInfo.datenum && S.calib.syncVersion == F.syncVersion   % same code, table and sync
        E = S.E; calib = S.calib; return                                           % done
    end
end

%% Pupil table
opts = detectImportOptions(tabF);                                                  % column layout
need = {'experiment_num','animal','insertion','frame_index','center_x_px','center_y_px','major_diameter_px','fit_valid','pupil_likelihood_threshold'};   % columns used
assert(all(ismember(need, opts.VariableNames)), 'getEyePosition:columns', ...
    'experimentNum %d: pupil table lacks %s.', experimentNum, strjoin(setdiff(need, opts.VariableNames), ', '));   % format changed
opts.SelectedVariableNames = need;                                                 % read only these
T = readtable(tabF, opts);                                                         % pupil table
assert(all(T.experiment_num == experimentNum) && all(string(T.animal) == r.animal) && all(T.insertion == r.insertion), ...
    'getEyePosition:key', 'experimentNum %d: table key does not match the manifest row.', experimentNum);   % right insertion
assert(T.frame_index(1) == 0 && all(diff(T.frame_index) == 1), 'getEyePosition:frames', ...
    'experimentNum %d: frame_index is not contiguous from 0.', experimentNum);    % row i = video frame i-1
nRows = height(T);                                                                 % frames in the table
assert(nRows <= F.nFrames, 'getEyePosition:tooManyRows', 'experimentNum %d: %d table rows > %d video frames.', experimentNum, nRows, F.nFrames);
assert(F.nFrames - nRows <= params.maxMissingRows, 'getEyePosition:incompleteTable', ...
    'experimentNum %d: table has %d of %d frames (incomplete).', experimentNum, nRows, F.nFrames);   % partial file
if nRows < F.nFrames                                                               % a few trailing frames missing
    warning('getEyePosition:shortTable', 'experimentNum %d: table lacks the last %d frames.', experimentNum, F.nFrames - nRows);
end
fitValid = strcmpi(string(T.fit_valid), "true") | string(T.fit_valid) == "1";      % parse True/False text

E = table();                                                                       % output
E.frame = T.frame_index;                                                           % 0-based video frame
E.tMs = F.frameTimeMs(E.frame + 1);                                                % ephys time [ms] (frame_index is 0-based)
E.xPx = T.center_x_px; E.yPx = T.center_y_px;                                      % pupil centre [px]
E.diamPx = T.major_diameter_px;                                                    % pupil major diameter [px]
E.valid = fitValid & isfinite(E.tMs) & isfinite(E.xPx) & isfinite(E.yPx);         % usable frames

%% Calibration from DLC landmarks
% The landmarks are static, so their DLC likelihood is often low (0.2-0.5) while the position is
% stable to < 1%. Accept a landmark length when it is stable across the sampled frames (IQR/median
% below maxLandmarkIQR) instead of filtering by likelihood (checked 2026-11-19 on 49, 57, 64, 69).
thr = T.pupil_likelihood_threshold(1);                                             % likelihood threshold Elena used for the pupil
D = readDlcSample(r.dlc_csv, params.dlcSampleFrames);                              % sampled DLC rows (bodypart x/y/likelihood)
[scaleLen, scaleLik] = landmarkDistance(D, "Rostral_scale", "Caudal_scale");       % retractor length per frame [px], median likelihood
[eyeWidth, edgeLik] = landmarkDistance(D, "Rostral_edge", "Caudal_edge");          % eye-opening width per frame [px], median likelihood

calib = struct();                                                                  % calibration and quality
calib.experimentNum = experimentNum; calib.animal = r.animal; calib.insertion = r.insertion;   % key
calib.eyeWidthPx = median(eyeWidth); calib.eyeWidthIQR = iqr(eyeWidth) / calib.eyeWidthPx;   % eye-opening width [px] and its relative spread
calib.retractorPx = median(scaleLen); calib.retractorIQR = iqr(scaleLen) / calib.retractorPx;   % retractor length [px] and its relative spread
calib.edgeLikelihood = edgeLik; calib.scaleLikelihood = scaleLik;                  % median DLC likelihoods (information)
if ~(calib.eyeWidthIQR < params.maxLandmarkIQR)                                    % eye corners not located consistently
    warning('getEyePosition:noEdges', 'experimentNum %d: eye-corner distance varies (IQR/median %.3f); no degrees.', experimentNum, calib.eyeWidthIQR);
    calib.eyeWidthPx = NaN;                                                        % degrees not available
end
if ~(calib.retractorIQR < params.maxLandmarkIQR)                                   % retractor ends not located consistently
    warning('getEyePosition:noScale', 'experimentNum %d: retractor length varies (IQR/median %.3f); no mm.', experimentNum, calib.retractorIQR);
    calib.retractorPx = NaN;                                                       % mm not available
end
calib.RPx = calib.eyeWidthPx / 2;                                                  % eye radius used for degrees [px]
calib.mmPerPx = params.retractorMm / calib.retractorPx;                            % mm per pixel (information only: varied 5.7-16.7 mm eye widths, 2026-11-19)
calib.eyeWidthMm = calib.eyeWidthPx * calib.mmPerPx;                               % eye-opening width [mm] (information only)
calib.pupilDiamPx = median(E.diamPx(E.valid), 'omitnan');                          % pupil diameter [px]
calib.pupilToEyeWidth = calib.pupilDiamPx / calib.eyeWidthPx;                      % sanity check: should be similar across insertions
calib.fracValid = mean(E.valid);                                                   % fraction of usable frames
valid2 = E.valid(1:end-1) & E.valid(2:end);                                        % consecutive valid frame pairs
step = hypot(diff(E.xPx), diff(E.yPx));                                            % frame-to-frame centre change [px]
calib.noiseFloorPx = median(step(valid2));                                         % typical frame-to-frame jitter [px]
calib.framePeriodMs = F.framePeriodMs;                                             % saved-frame period [ms]
calib.likelihoodThreshold = thr;                                                   % pupil likelihood threshold of the table
calib.eyeVersion = eyeVersion; calib.syncVersion = F.syncVersion;                  % provenance
calib.tableDatenum = tabInfo.datenum; calib.tableFile = tabF;                      % table used

if ~isfolder(params.cacheDir), mkdir(params.cacheDir); end                         % cache folder
save(cacheF, 'E', 'calib');                                                        % cache
fprintf('experimentNum %d (%s ins %d): %.0f%% valid frames, pupil %.1f px, eye width %.1f px (%.2f mm), R %.1f px, noise %.3f px/frame.\n', ...
    experimentNum, r.animal, r.insertion, 100 * calib.fracValid, calib.pupilDiamPx, calib.eyeWidthPx, calib.eyeWidthMm, calib.RPx, calib.noiseFloorPx);
end

function D = readDlcSample(f, nPer)
% Read nPer data rows from the start, middle and 3/4 of a DLC CSV (3 header lines).
fid = fopen(f, 'r');                                                               % open
assert(fid > 0, 'getEyePosition:open', 'Cannot open %s', f);                       % unreadable
cleaner = onCleanup(@() fclose(fid));                                              % always close
fgetl(fid);                                                                        % line 1: scorer
parts = strtrim(split(string(fgetl(fid)), ","));                                   % line 2: bodyparts (trim a trailing \r of CRLF files)
coords = strtrim(split(string(fgetl(fid)), ","));                                  % line 3: x / y / likelihood
dataStart = ftell(fid);                                                            % byte where the data begin
fmt = repmat('%f', 1, numel(parts));                                               % all columns numeric
fseek(fid, 0, 'eof'); fileEnd = ftell(fid);                                        % file size
rows = [];                                                                         % sampled rows
for frac = [0 0.5 0.75]                                                            % three places in the session
    fseek(fid, dataStart + round(frac * (fileEnd - dataStart)), 'bof');            % jump
    if frac > 0, fgetl(fid); end                                                   % skip the partial line
    C = textscan(fid, fmt, nPer, 'Delimiter', ',', 'CollectOutput', true);        % read nPer rows
    rows = [rows; C{1}]; %#ok<AGROW>                                               % append
end
D = struct('parts', parts, 'coords', coords, 'data', rows);                        % column names + values
end

function [L, lik] = landmarkDistance(D, a, b)
% Per-frame distance between two DLC bodyparts, and their median likelihood.
col = @(p, c) find(D.parts == p & D.coords == c);                                  % column of bodypart p, coordinate c
need = [col(a,"x") col(a,"y") col(a,"likelihood") col(b,"x") col(b,"y") col(b,"likelihood")];   % the six columns
assert(numel(need) == 6, 'getEyePosition:bodypart', 'DLC file lacks %s or %s.', a, b);   % model without these points
X = D.data(:, need);                                                               % their values
L = hypot(X(:,1) - X(:,4), X(:,2) - X(:,5));                                       % distance [px], every sampled frame
lik = median(min(X(:,3), X(:,6)));                                                 % median of the weaker point's likelihood
end
