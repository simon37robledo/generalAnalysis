function R = checkVideoStimTiming(experimentNum, vsObj, params)
% CHECKVIDEOSTIMTIMING  Independent check of the eye-video clock: does the video brightness change
% at the diode-measured stimulus onsets?
%
%   R = checkVideoStimTiming(experimentNum, vsObj)
%
%   Light from the screen reaches the eye camera (corneal reflection, eyelid, retractor), so if
%   the frame times from buildVideoFrameTimes are right, the mean frame brightness averaged over
%   many trials should step within about one frame of the stimulus onset
%   (vsObj.getSyncedDiodeTriggers, measured with the photodiode). A constant offset of the video
%   clock shifts that step by the same amount. If the camera does not see the screen (IR filter),
%   no step appears and the check is inconclusive (R.stepFound = false).
%
%   Output R: time axis [ms], mean brightness (baseline-subtracted), step latency [ms], frames used.

arguments
    experimentNum (1,1) double                                                   % Excel row
    vsObj                                                                        % stimulus object of that recording
    params.nTrials double = 300                                                  % trials sampled
    params.winMs double = [-200 400]                                             % window around onset [ms]
    params.manifest string = "W:\Large_scale_mapping_NP\EyeTracking\manifest.csv"   % video path
    params.plot logical = true                                                   % draw the average
end

F = buildVideoFrameTimes(experimentNum);                                         % frame times [ms]
M = readtable(params.manifest, 'TextType', 'string');                            % manifest
r = M(M.experimentNum == experimentNum, :);                                      % this insertion
V = VideoReader(r.dlc_video);                                                    % the video DLC ran on
S = vsObj.getSyncedDiodeTriggers;                                                % diode-measured onsets [ms]
on = S.stimOnFlipTimes(:);                                                       % onsets
inVid = on + params.winMs(1) > min(F.frameTimeMs) & on + params.winMs(2) < max(F.frameTimeMs);   % onsets inside the video
on = on(inVid);                                                                  % usable
pick = round(linspace(1, numel(on), min(params.nTrials, numel(on))));            % evenly spread sample of trials
fp = F.framePeriodMs;                                                            % frame period [ms]
edges = params.winMs(1):fp:params.winMs(2);                                      % one bin per frame period
sumB = zeros(1, numel(edges) - 1); nB = zeros(1, numel(edges) - 1);              % accumulators
for i = pick
    k = find(F.frameTimeMs >= on(i) + params.winMs(1) & F.frameTimeMs < on(i) + params.winMs(2));   % frames in the window (1-based)
    if isempty(k), continue, end                                                 % none
    fr = read(V, [k(1) k(end)]);                                                 % contiguous frames (H x W x 3 x n)
    b = squeeze(mean(fr, [1 2 3]))';                                             % mean brightness per frame
    rel = F.frameTimeMs(k(1):k(end))' - on(i);                                   % frame time relative to onset [ms]
    bin = discretize(rel, edges);                                                % bin of each frame
    okb = ~isnan(bin);                                                           % inside the window
    sumB = sumB + accumarray(bin(okb)', b(okb)', [numel(nB) 1])';                % add brightness
    nB = nB + accumarray(bin(okb)', 1, [numel(nB) 1])';                          % add counts
end
tc = edges(1:end-1) + fp / 2;                                                    % bin centres [ms]
m = sumB ./ nB;                                                                  % mean brightness per bin
pre = tc < -fp;                                                                  % bins clearly before onset
m = m - mean(m(pre));                                                            % baseline-subtracted
sdPre = std(m(pre));                                                             % baseline noise
thr = 5 * sdPre;                                                                 % step threshold
post = find(tc > -2 * fp & abs(m) > thr, 1);                                     % first bin past threshold near/after onset
R = struct('experimentNum', experimentNum, 'stim', string(class(vsObj)), 't', tc, 'brightness', m, ...
    'nTrials', numel(pick), 'sdPre', sdPre, 'stepFound', ~isempty(post), 'latencyMs', NaN, 'framePeriodMs', fp);   % result
if ~isempty(post), R.latencyMs = tc(post); end                                   % first frame-bin with the step (bin centre)
if params.plot
    figure('Position', [100 100 500 320]); stairs(edges, [m m(end)], 'k'); hold on;   % bins as steps between edges
    xline(0, 'r'); yline([-thr thr], 'k:');                                      % onset, threshold
    xlabel('Time from diode onset (ms)'); ylabel('Mean frame brightness (a.u., baseline 0)');
    title(sprintf('%d %s: step at %.0f ms (frame %.0f ms)', experimentNum, R.stim, R.latencyMs, fp), 'Interpreter', 'none');
end
fprintf('experimentNum %d %s: %d trials, step %s at %.1f ms (frame period %.1f ms, timing uncertainty %.0f ms).\n', ...
    experimentNum, R.stim, R.nTrials, string(R.stepFound), R.latencyMs, fp, F.timeUncertaintyMs);
end
