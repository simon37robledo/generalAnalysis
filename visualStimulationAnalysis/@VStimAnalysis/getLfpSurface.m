function S = getLfpSurface(obj, params)
% getLfpSurface  Brain-surface position (probe y, um from tip) for this insertion from the LFP.
%
%   S = obj.getLfpSurface()
%
%   Stimulus recordings span at most ~2.1 mm above the tip, below the brain surface, so the
%   LFP boundary (estimateLfpSurface) is measured on the animal's whole-shank / depth-test
%   recordings (folders matching 'wholeShank' or 'depthTest' under the animal folder).
%   Each test is matched to an insertion: by the insertion number in its folder name
%   ('<n>_depthTest'), otherwise by the nearest recording in time. Per test,
%   delta = LFP surface - coor_Z of the matched insertion. Per animal,
%   offset = mean delta over valid tests (boundary not at the column edge, bands agree,
%   coor_Z known). For this insertion: surfaceY = coor_Z + offset (source "lfp-offset").
%   An animal without a valid test keeps coor_Z (source "coorZ-noTest", warning).
%
%   Output S: exp (Excel row), surfaceY, coorZ, offsetUm, source, tests (table).
%   Each test's estimate is cached as lfpSurface.mat in the test's folder.

arguments
    obj
    params.excelPath char = '\\sil3\data\Large_scale_mapping_NP\Experiment_Excel.xlsx'   % master table (coor_Z)
    params.overwrite (1,1) logical = false          % recompute the cached test estimates
    params.maxDeltaUm (1,1) double = 300            % |LFP - coor_Z| above this -> warning
    params.verbose (1,1) logical = true             % print the test table
end
lfpVersion = 1;                                     % bump when estimateLfpSurface changes results (invalidates caches)

%% This insertion's row in the master table
X = readtable(params.excelPath);                                       % master table
rowName = string(X.Exp_name) + "_" + string(X.Insertion);              % '<Exp_name>_<Insertion>' per row (as recordingName)
ex = find(rowName == string(obj.dataObj.recordingName));               % row(s) of this recording
assert(isscalar(ex), 'getLfpSurface:row', '%d Excel rows match recording %s (need exactly 1).', ...
    numel(ex), obj.dataObj.recordingName);                             % alignment invariant
animal = string(X.Animal_ID(ex));                                      % animal ID
rows = find(string(X.Animal_ID) == animal);                            % all insertions of this animal

%% The animal's test recordings
recDir = string(obj.dataObj.recordingDir);                             % ...\<animal>\<experiment>\Insertion<n>\catgt_...
animalDir = fileparts(fileparts(fileparts(recDir)));                   % ...\<animal>
f = dir(fullfile(animalDir, '**', '*.ap.meta'));                       % every AP meta file of the animal
isTest = ~cellfun(@isempty, regexpi({f.folder}, 'wholeshank|depth_?test', 'once'));   % test folders only
isCatGT = contains({f.name}, '_tcat.');                                % CatGT copies of the same recordings (time-shifted data)
f = f(isTest & ~isCatGT);                                              % raw SpikeGLX test meta files only (as the notebook)

%% Start/end times of the animal's insertions (catgt meta)
insStart = NaT(numel(rows), 1); insEnd = NaT(numel(rows), 1);          % recording interval per insertion
for k = 1:numel(rows)
    r = rows(k);                                                       % Excel row
    d = fullfile(animalDir, string(X.Exp_name(r)), "Insertion" + X.Insertion(r), ...
        "catgt_" + string(X.Exp_name(r)) + "_" + X.Insertion(r) + "_g0");   % catgt folder (as loadNPclassFromTable)
    mf = dir(fullfile(d, '*.ap.meta'));                                % its AP meta
    if isempty(mf), continue, end                                      % insertion without data: cannot be matched
    [insStart(k), insEnd(k)] = metaInterval(fullfile(mf(1).folder, mf(1).name));   % start / end time
end

%% Estimate (or load) each test and match it to an insertion
T = table();                                                           % one row per test
for k = 1:numel(f)
    metaPath = string(fullfile(f(k).folder, f(k).name));               % test meta file
    cacheF = fullfile(f(k).folder, 'lfpSurface.mat');                  % cached estimate
    R = [];                                                            % estimate of this test
    if ~params.overwrite && isfile(cacheF)                             % reuse a cache of the current version
        C = load(cacheF, 'R', 'lfpVersion');                           % cached estimate and its version
        if isfield(C, 'lfpVersion') && C.lfpVersion == lfpVersion, R = C.R; end   % stale version -> recompute
    end
    if isempty(R)
        try
            R = estimateLfpSurface(metaPath);                          % LFP boundary of this test
        catch ME
            if ~strcmp(ME.identifier, 'estimateLfpSurface:fewSites'), rethrow(ME); end   % only "too few sites" is tolerated
            warning('getLfpSurface:testSkipped', 'Test %s skipped: %s', f(k).name, ME.message);   % explicit
            continue
        end
        save(cacheF, 'R', 'lfpVersion');                               % cache next to the test recording
    end

    [tS, tE] = metaInterval(metaPath);                                 % test recording interval
    gapH = max(max(hours(insStart - tE), hours(tS - insEnd)), 0);      % hours between test and each insertion (0 = overlap)
    [gap, kNear] = min(gapH);                                          % nearest insertion in time (NaT rows -> NaN, ignored)
    folderRest = erase(string(f(k).folder), unique(string(X.Exp_name(rows))));   % drop the experiment name (its date ends in digits)
    nTok = regexpi(folderRest, '_(\d+)_depth_?test', 'tokens', 'once');  % insertion number in the folder name, if any
    if ~isempty(nTok)                                                  % named test: use the name
        kM = find(X.Insertion(rows) == str2double(nTok{1}));           % that insertion of this animal
        assert(isscalar(kM), 'getLfpSurface:testName', 'Test %s names insertion %s, not found once.', f(k).name, nTok{1});
    else
        kM = kNear;                                                    % unnamed test: nearest in time
    end
    if ~isempty(nTok) && X.coor_Z(rows(kM)) ~= X.coor_Z(rows(kNear))   % name and time point to different penetrations
        warning('getLfpSurface:matchConflict', 'Test %s: name -> insertion %d, time -> insertion %d (different coor_Z).', ...
            f(k).name, X.Insertion(rows(kM)), X.Insertion(rows(kNear)));
    end
    coorZ = X.coor_Z(rows(kM));                                        % manipulator depth of the matched insertion
    delta = R.surfaceMeanUm - coorZ;                                   % LFP surface minus coor_Z (um)
    valid = ~R.nearEdge && ~R.bandsDisagree && isfinite(coorZ);        % usable for the animal offset
    if valid && abs(delta) > params.maxDeltaUm                         % large disagreement: report
        warning('getLfpSurface:largeDelta', 'Test %s: LFP %.0f vs coor_Z %.0f um (delta %.0f).', ...
            f(k).name, R.surfaceMeanUm, coorZ, delta);
    end
    T = [T; table(string(f(k).name), rows(kM), X.Insertion(rows(kM)), gap, R.shank, R.yUm(end), ...
        R.surfaceMeanUm, R.surfaceMedianUm, R.surfaceStdUm, coorZ, delta, R.nearEdge, R.bandsDisagree, valid, ...
        'VariableNames', {'test','matchedExp','matchedInsertion','gapHours','shank','topSiteUm', ...
        'lfpMeanUm','lfpMedianUm','lfpSdUm','coorZ','deltaUm','nearEdge','bandsDisagree','valid'})];   %#ok<AGROW>
end

%% Per-animal offset and this insertion's surface
S = struct('exp', ex, 'coorZ', X.coor_Z(ex), 'tests', T);              % output
if ~isempty(T) && any(T.valid)
    S.offsetUm = mean(T.deltaUm(T.valid));                             % mean LFP - coor_Z over valid tests
    S.source = "lfp-offset";                                           % calibrated depth
else
    S.offsetUm = 0;                                                    % no calibration available
    S.source = "coorZ-noTest";                                         % falls back to the manipulator depth
    warning('getLfpSurface:noTest', '%s (exp %d): no valid depth test; using coor_Z.', animal, ex);
end
S.surfaceY = S.coorZ + S.offsetUm;                                     % surface position on the probe (um from tip)
assert(isfinite(S.surfaceY), 'getLfpSurface:noCoorZ', 'Exp %d: coor_Z missing, no surface position.', ex);
if params.verbose
    fprintf('%s exp %d: coor_Z %.0f, offset %.1f um (%s) -> surface %.0f um\n', ...
        animal, ex, S.coorZ, S.offsetUm, S.source, S.surfaceY);        % result
    disp(T);                                                           % per-test table
end
end

% =========================================================================
function [t0, t1] = metaInterval(metaPath)
% Start and end time of a SpikeGLX recording (original creation time if catgt rewrote it).
txt = fileread(metaPath);                                              % meta file text
c = regexp(txt, 'fileCreateTime_original=([^\r\n]+)', 'tokens', 'once');   % catgt keeps the original time here
if isempty(c), c = regexp(txt, 'fileCreateTime=([^\r\n]+)', 'tokens', 'once'); end   % raw recordings
dur = regexp(txt, 'fileTimeSecs=([^\r\n]+)', 'tokens', 'once');        % recording length (s)
assert(~isempty(c) && ~isempty(dur), 'getLfpSurface:meta', 'No creation time / length in %s', metaPath);
t0 = datetime(c{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');         % start time
t1 = t0 + seconds(str2double(dur{1}));                                 % end time
end
