%% SYNC — trial timing from the trial TTL + the recording's measured TTL->diode lag (run before the figures)
% Where the per-trial diode pairing is unreliable, each trial onset is its trial TTL plus the
% recording's median TTL->diode lag, measured on the trials where the onset is visible in the
% diode (getSyncedDiodeTriggers 'trialOnsetFromTTL'; >= 95% of those within one frame, asserted).
%   - Gratings (SDG), every recording: the static onset is a large diode drop ~80-95 ms after the
%     TTL, visible in about every other trial (when the previous trial ended on a bright frame).
%     The older per-trial pairing put recordings 1 and 8 one static period early (all trials),
%     placed ~30/80 trials of PV97 / SA9 at the TTL, and its caches could not be regenerated.
%   - SA5 (rows 28-32, 2023 protocol): FFF, MB, RG and SDG (MB: diode active before ~40% of
%     onsets; RG/SDG: missing or extra flips near some onsets).
% Existing caches are backed up to <analysis folder>\backup_preTTLsync_20261105 (never overwritten).
% Output: one row per recording x stimulus in Paper_figs\Stats summary tables\SyncTTLLag_lags.csv.

d = figureDesign();                                               % recording lists and output folder
sdgRows = [d.exps, 56];                                           % gratings: every Pogona recording + SA9 ins 2
jobs = [num2cell(sdgRows'), repmat({{'SDG'}}, numel(sdgRows), 1); ...
        num2cell((28:32)'), repmat({{'FFF', 'MB', 'RG', 'SDG'}}, 5, 1)];   % SA5: all four stimuli
ctor = struct('FFF', @fullFieldFlashAnalysis, 'MB', @linearlyMovingBallAnalysis, ...
    'RG', @rectGridAnalysis, 'SDG', @StaticDriftingGratingAnalysis);   % analysis class per stimulus
bkName = 'backup_preTTLsync_20261105';                            % backup subfolder
rows = {};                                                        % lag table
for j = 1:size(jobs, 1)
    ex = jobs{j, 1}; NP = loadNPclassFromTable(ex);               % recording
    for s = jobs{j, 2}
        try
            o = ctor.(s{1})(NP);                                  % stimulus object (errors if not recorded)
        catch
            continue                                              % stimulus not recorded in this recording
        end
        if isempty(o.VST), continue, end                          % no stimulus file
        f = o.visualStimAnalysisFolder;                           % cache folder
        bk = fullfile(f, bkName); if ~isfolder(bk), mkdir(bk); end
        for c = {'getSyncedDiodeTriggers.mat', 'getDiodeTriggers.mat'}
            if isfile(fullfile(f, c{1})) && ~isfile(fullfile(bk, c{1}))
                copyfile(fullfile(f, c{1}), bk);                  % keep the previous cache
            end
        end
        o.getSessionTime;                                         % session limits (cached)
        o.getSyncedDiodeTriggers('overwrite', true, 'trialOnsetFromTTL', true);   % rebuild the onsets
        sy = load(fullfile(f, 'getSyncedDiodeTriggers.mat'));     % new cache
        assert(numel(sy.stimOnFlipTimes) == o.VST.nTotTrials, 'row %d %s: %d onsets, expected %d trials', ...
            ex, s{1}, numel(sy.stimOnFlipTimes), o.VST.nTotTrials);   % one onset per trial
        L = sy.ttlSync.lagOn(sy.ttlSync.cleanTrials);             % measured lags (ms)
        rows(end+1, :) = {ex, string(NP.recordingName), string(s{1}), o.VST.nTotTrials, numel(sy.ttlSync.cleanTrials), ...
            sy.ttlSync.medianOnLag, prctile(L, 5), prctile(L, 95), 100*sy.ttlSync.fracCleanWithin1Frame, sy.ttlSync.medianOffLag}; %#ok<SAGROW>
        fprintf('SYNC row %3d %-4s %3d trials | lag %.1f ms (%d measured, %.0f%% within 1 frame)\n', ex, s{1}, ...
            o.VST.nTotTrials, sy.ttlSync.medianOnLag, numel(sy.ttlSync.cleanTrials), 100*sy.ttlSync.fracCleanWithin1Frame);
    end
end
T = cell2table(rows, 'VariableNames', {'exp', 'recording', 'stimulus', 'nTrials', 'nMeasured', 'lagMedianMs', ...
    'lagP5Ms', 'lagP95Ms', 'pctWithin1Frame', 'offLagMedianMs'});
writetable(T, fullfile(d.paperFigs, 'Stats summary tables', 'SyncTTLLag_lags.csv'));   % methods table
disp(T);
