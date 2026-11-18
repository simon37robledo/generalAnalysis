function M = makeEyeTrackingManifest(params)
% MAKEEYETRACKINGMANIFEST  One row per insertion with an eye video: which DLC output to analyse,
% which video it was run on, and the frame counts needed to sync eye frames to ephys time.
%
%   M = makeEyeTrackingManifest()                  % build and print, write nothing
%   M = makeEyeTrackingManifest(write=true)        % also write <outDir>\manifest.csv
%
%   Key: experimentNum (row of Experiment_Excel.xlsx, 1-based as in MATLAB) + animal + insertion.
%   The same experimentNum is the insertion code used by hierBootMatchFreq. Never key by a
%   0-based pandas row: pandas row k is MATLAB row k+1.
%
%   DLC choice rule (one raw DLC CSV per insertion):
%     1. only CSVs written by DLC for this insertion's video ("<video stem>DLC_<scorer>.csv"),
%        excluding re-analyses of labeled videos and derived files (_filtered, _inverted,
%        -FilteredJitter, ellipse / good-frame tables);
%     2. newest network first: IR_Eye_movementsApr29 > EyeTracking3Oct3 > any other;
%     3. highest snapshot;
%     4. tie between shuffles -> the shuffle that already has an ellipse fit (it was inspected),
%        else the lowest shuffle; the tie is written to `note` so it can be reviewed.
%
%   Status columns:
%     status_elena  "ready" when the chosen DLC CSV exists and has one row per saved video frame
%                   (DLC rows == camera-timestamp rows); otherwise the reason (needs_dlc, ...).
%     status_sync   "ok" when videoTimeStampsSynced.mat (video_analysis.m, dropped frames removed)
%                   has exactly one TTL per saved frame; otherwise ttl_mismatch / no_sync /
%                   camera_ts_csv x0. Only our side: it does not stop Elena from analysing the video.

arguments
    params.excelFile string = "W:\Large_scale_mapping_NP\Experiment_Excel.xlsx"   % experiment sheet
    params.outDir string = "W:\Large_scale_mapping_NP\EyeTracking"                % handoff folder
    params.rows double = []                                                        % Excel rows to include (default: every row with an eye-video folder)
    params.write logical = false                                                   % write manifest.csv to outDir
end

X = readtable(params.excelFile, 'Format', 'auto');                     % experiment sheet, one row per insertion
eyeDir = string(X.Eye_video_dir);                                      % eye-video folder per row
hasDir = ~ismissing(eyeDir) & strlength(eyeDir) > 0 & eyeDir ~= "NaN"; % rows that list a folder
rows = find(hasDir)';                                                  % candidate rows
if ~isempty(params.rows), rows = intersect(rows, params.rows, 'stable'); end   % user subset
assert(~isempty(rows), 'makeEyeTrackingManifest:noRows', 'No Excel rows with an eye-video folder selected.');

modelRank = dictionary(["IR_Eye_movementsApr29" "EyeTracking3Oct3"], [2 1]);   % newer DLC network ranks higher

C = cell(numel(rows), 1);                                              % one table row per insertion
for k = 1:numel(rows)
    ex = rows(k);                                                      % Excel row = experimentNum
    ins = X.Insertion(ex);                                             % insertion number within the animal
    R = struct();                                                      % this insertion's manifest row
    R.experimentNum = ex;                                              % globally unique insertion code
    R.animal = string(X.Animal_ID(ex));                                % animal ID
    R.insertion = ins;                                                 % insertion number
    R.video_dir = toW(eyeDir(ex));                                     % folder with video and DLC files (W: path)
    R.dlc_csv = ""; R.dlc_video = ""; R.camera_ts_csv = "";            % chosen files, filled below
    R.dlc_model = ""; R.dlc_shuffle = NaN; R.dlc_snapshot = NaN;       % chosen network
    R.n_dlc_frames = NaN; R.n_camera_frames = NaN; R.n_synced_ttl = NaN;   % frame counts to compare
    R.frame_w = NaN; R.frame_h = NaN;                                  % pixel size of the analysed video
    R.had_inverted_file = false;                                       % an older "_inverted" copy exists (mirror/upside-down setup)
    R.n_dlc_candidates = 0; R.note = "";                               % review fields
    R.status_elena = ""; R.status_sync = ""; R.ttl_minus_frames = NaN; % input readiness for Elena / ephys-sync readiness for us (TTLs - saved frames)
    R.rec_dir = toW(fullfile(string(X.Base_path(ex)), string(X.Exp_name(ex)), "Insertion" + ins, ...
        "catgt_" + string(X.Exp_name(ex)) + "_" + ins + "_g0"));       % catgt folder (as loadNPclassFromTable)

    problems = strings(0);                                             % reasons Elena cannot analyse this row yet
    syncProblems = strings(0);                                         % reasons we cannot put its frames on ephys time yet
    notes = strings(0);                                                % information, not blocking
    if ~isfolder(R.video_dir)                                          % folder listed but not on disk
        warning('makeEyeTrackingManifest:noDir', 'Row %d (%s ins %d): folder missing: %s', ex, R.animal, ins, R.video_dir);
        R.status_elena = "missing_dir"; R.status_sync = "unchecked"; C{k} = R; continue   % nothing else can be checked
    end
    L = dir(R.video_dir); L = L(~[L.isdir]);                           % files in the folder
    names = string({L.name});                                          % file names
    tok = regexp(names, '(?i)(?:_|insertion)(\d+)_Camera', 'tokens', 'once');   % insertion number written just before "_Camera"
    insOf = nan(size(names));                                          % insertion number per file (NaN = none)
    for j = 1:numel(tok), if ~isempty(tok{j}), insOf(j) = str2double(tok{j}{1}); end, end
    mine = insOf == ins;                                               % files of this insertion only

    % Camera timestamps written by the acquisition software (plain "<...>Camera1_<date>-<time>.csv")
    isCam = mine & ~cellfun(@isempty, regexp(names, 'Camera\d_\d{8}-\d{6}\.csv$', 'once'));   % plain timestamp CSV
    if nnz(isCam) == 1                                                 % exactly one, as required
        R.camera_ts_csv = fullfile(R.video_dir, names(isCam));         % full path
        hdr = string(readlines(R.camera_ts_csv, 'EmptyLineRule', 'skip', 'WhitespaceRule', 'trim'));   % (small file) header + rows
        assert(contains(hdr(1), "timestamp"), 'makeEyeTrackingManifest:tsColumn', ...
            'Row %d: camera CSV has no "timestamp" column: %s', ex, R.camera_ts_csv);   % unexpected format
        R.n_camera_frames = numel(hdr) - 1;                            % frames saved in the video (rows minus header)
    else
        syncProblems(end+1) = sprintf("camera_ts_csv x%d", nnz(isCam)); % none or ambiguous: frame count cannot be verified
    end

    % DLC candidates: "<video stem>DLC_<scorer>.csv", raw only
    isDlc = mine & endsWith(names, ".csv") & contains(names, "DLC_");  % DLC-written or derived CSVs of this insertion
    R.had_inverted_file = any(isDlc & contains(names, "_inverted"));   % hint that coordinates were flipped before
    derived = contains(names, ["_filtered" "_inverted" "FilteredJitter" "ellipse" "good_frames" "GazeMapper" "labeled"], 'IgnoreCase', true) ...
        | startsWith(names, "Diff vibration");                         % not raw DLC output (the last: a stray figure file named by video_analysis.m's print(sprintf(...)))
    cand = find(isDlc & ~derived);                                     % raw DLC CSVs
    R.n_dlc_candidates = numel(cand);                                  % how many raw outputs exist
    if isempty(cand)                                                   % no DLC run for this insertion
        problems(end+1) = "needs_dlc";                                 % Elena has nothing to analyse
    else
        P = regexp(names(cand), '^(?<stem>.*)DLC_Resnet\d+_(?<model>.+?)shuffle(?<shuffle>\d+)_snapshot_(?<snap>\d+)\.csv$', 'names');   % parse scorer
        if iscell(P), P = [P{:}]; end                                  % struct array, one per candidate (regexp returns a cell only for >1 name)
        assert(numel(P) == numel(cand), 'makeEyeTrackingManifest:parse', ...
            'Row %d: could not parse a DLC file name in %s', ex, strjoin(names(cand), ', '));   % new naming scheme
        rank = zeros(1, numel(P));                                     % network rank per candidate
        for j = 1:numel(P), if isKey(modelRank, P(j).model), rank(j) = modelRank(P(j).model); end, end
        snap = str2double([P.snap]); shuf = str2double([P.shuffle]);   % snapshot and shuffle numbers (string array; a cell of strings would give NaN)
        assert(all(~isnan([snap shuf])), 'makeEyeTrackingManifest:parseNum', 'Row %d: unreadable shuffle/snapshot number.', ex);   % parse check
        best = find(rank == max(rank));                                % newest network
        best = best(snap(best) == max(snap(best)));                    % highest snapshot of it
        if numel(best) > 1                                             % shuffles tie
            fitted = arrayfun(@(j) any(contains(names, P(j).stem + "DLC_") & contains(names, "shuffle" + P(j).shuffle + "_snapshot_" + P(j).snap) ...
                & contains(names, "ellipse", 'IgnoreCase', true)), best);   % shuffle already used for an ellipse fit
            tieMsg = "tie shuffles " + strjoin(string(shuf(best)), "/");    % for review
            if nnz(fitted) == 1, best = best(fitted); tieMsg = tieMsg + " -> " + shuf(best) + " (has earlier ellipse fit)";   % inspected before
            else, [~, m] = min(shuf(best)); best = best(m); tieMsg = tieMsg + " -> " + shuf(best) + " (lowest)";   % arbitrary but fixed
            end
            notes(end+1) = tieMsg;                                     % flag for review
        end
        R.dlc_csv = fullfile(R.video_dir, names(cand(best)));          % chosen DLC output
        R.dlc_model = P(best).model; R.dlc_shuffle = shuf(best); R.dlc_snapshot = snap(best);   % chosen network
        R.n_dlc_frames = countDataRows(R.dlc_csv);                     % frames DLC analysed
        vid = fullfile(R.video_dir, P(best).stem + ".mp4");            % the video DLC was run on (its pixel frame)
        if isfile(vid)                                                 % present
            R.dlc_video = vid;                                         % record it
            V = VideoReader(vid);                                      % header only
            R.frame_w = V.Width; R.frame_h = V.Height;                 % pixel size of the DLC coordinates
        else
            problems(end+1) = "dlc_video_missing";                     % cannot check pixel frame
        end
        if isempty(regexp(P(best).stem, '^(?i)(pv|sa)\d+_.*Camera\d_\d{8}-\d{6}$', 'once'))   % DLC run on a processed copy (cropped/enhanced), not the raw camera video
            notes(end+1) = "processed video (cropped/enhanced): pixels are in that frame";   % fine for Elena, matters for degrees
        end
        if ~isnan(R.n_camera_frames) && R.n_dlc_frames ~= R.n_camera_frames   % DLC did not analyse every saved frame
            problems(end+1) = "dlc_incomplete";                        % rows would not map to frames
            warning('makeEyeTrackingManifest:dlcIncomplete', 'Row %d (%s ins %d): DLC %d rows, video %d frames.', ...
                ex, R.animal, ins, R.n_dlc_frames, R.n_camera_frames); % loud
        end
    end

    % Synced camera TTLs (one ephys time per saved frame)
    syncF = fullfile(R.rec_dir, "videoTimeStampsSynced.mat");          % written by video_analysis.m
    if isfile(syncF)
        S = load(syncF, 'vidTTLsR');                                   % TTL times of saved frames (ms)
        R.n_synced_ttl = numel(S.vidTTLsR);                            % number of synced frames
        R.ttl_minus_frames = R.n_synced_ttl - R.n_camera_frames;       % 0 = one TTL per saved frame
        if R.ttl_minus_frames ~= 0                                     % includes NaN camera count (comparison false) -> not flagged here
            syncProblems(end+1) = "ttl_mismatch";                      % index-based frame->time mapping would be off
        end
    else
        syncProblems(end+1) = "no_sync";                               % video_analysis.m not run for this insertion
    end

    if isempty(problems), R.status_elena = "ready"; else, R.status_elena = strjoin(problems, ";"); end   % Elena's input
    if isempty(syncProblems), R.status_sync = "ok"; else, R.status_sync = strjoin(syncProblems, ";"); end   % our side
    R.note = strjoin(notes, "; ");                                     % review notes
    C{k} = R;                                                          % store
end

M = struct2table([C{:}]', 'AsArray', true);                            % manifest table
assert(numel(unique(M.experimentNum)) == height(M), 'makeEyeTrackingManifest:dupKey', 'Duplicate experimentNum in manifest.');   % unique key
M.out_table = "eye_" + M.experimentNum + "_" + M.animal + "_ins" + M.insertion + ".csv";   % file name Elena returns per row

disp(M(:, {'experimentNum' 'animal' 'insertion' 'dlc_model' 'dlc_shuffle' 'dlc_snapshot' 'n_dlc_candidates' ...
    'n_dlc_frames' 'n_camera_frames' 'ttl_minus_frames' 'frame_w' 'frame_h' 'had_inverted_file' 'status_elena' 'status_sync' 'note'}));   % review print
fprintf('%d insertions: %d ready for Elena, %d with ephys sync ok.\n', height(M), nnz(M.status_elena == "ready"), nnz(M.status_sync == "ok"));   % summary

if params.write                                                        % handoff file
    if ~isfolder(params.outDir), mkdir(params.outDir); end             % create handoff folder
    if ~isfolder(fullfile(params.outDir, "tables")), mkdir(fullfile(params.outDir, "tables")); end   % where Elena writes her tables
    outF = fullfile(params.outDir, "manifest.csv");                    % manifest path
    writetable(M, outF);                                               % write
    fprintf('Manifest written: %s\n', outF);                           % confirm
end
end

function p = toW(p)
% Map the \\sil3\data share to the W: drive letter (the same storage).
p = string(p);                                                         % work with strings
p = regexprep(p, '^\\\\sil3\\data\\', 'W:\\', 'ignorecase');           % \\sil3\data\... -> W:\...
end

function n = countDataRows(f)
% Number of data rows in a DLC CSV (3 header lines: scorer, bodyparts, coords).
fid = fopen(f, 'r');                                                   % open read-only
assert(fid > 0, 'makeEyeTrackingManifest:open', 'Cannot open %s', f);  % unreadable file
cleaner = onCleanup(@() fclose(fid));                                  % always close
nl = 0; last = uint8(10);                                              % newline count, last byte read
while true
    b = fread(fid, 2^24, '*uint8');                                    % 16 MB chunk
    if isempty(b), break, end                                          % end of file
    nl = nl + nnz(b == 10);                                            % count line feeds
    last = b(end);                                                     % remember final byte
end
nLines = nl + double(last ~= 10);                                      % a last line without trailing newline still counts
n = nLines - 3;                                                        % minus the 3 DLC header lines
end
