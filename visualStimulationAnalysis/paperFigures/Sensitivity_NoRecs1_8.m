%% SENSITIVITY — every cross-animal comparison without recordings 1 and 8
% Recordings 1 (PV103, insertion 1) and 8 (PV27, insertion 2) use the old protocol:
% single-trajectory moving ball (8 directions, black background), 4x4 static-ball
% grid, gratings at 20% contrast with an unsaved temporal frequency. They are kept
% in the figures (MB everywhere except Figure 2); this script reports every
% comparison with them removed. Each is its animal's only insertion, so removing
% them = removing PV103 and PV27. Statistics are re-bootstrapped from the saved
% figure structs (buildStatsTable, the routine behind the figures and the LOAO
% table); per-recording statistics are not recomputed. Run after Figs 1, 3, 4, 5.

d = figureDesign();                                               % shared design (recordings, output folder)
drop = ["PV103", "PV27"];                                         % animals of recordings 1 and 8
nBoot = 1e4;                                                      % bootstrap replicates (as the LOAO table)
outCsv = fullfile(d.paperFigs, 'Stats summary tables', 'Sensitivity_without_recordings_1_8.csv');   % output table

%% Cross-animal comparisons (Figures 1, 3, 4)
spec = { ...                                                      % {label, struct, pairs}, as in the figure scripts
    'Figure 1',        d.structPath(d.exps, {'MB','RG'}),    "all"; ...
    'Figure 3',        d.structPath(d.exps, {'SDGm','SDGs'}), "all"; ...
    'Figure 4 static', d.structPath(d.exps, {'RG','SDGs'}),  "all"; ...
    'Figure 4 MB-SDGm', d.structPath(d.exps, {'MB','SDGm'}), "all"; ...
    'Figure 4 MB-MBR', d.structPath(d.exps, {'MB','MBR'}),   "all"};
metrics = {'Z-score', 'DeltaRate', 'FracResponsive'};             % the reported metrics
R = table();                                                      % one row per comparison x metric
for k = 1:size(spec, 1)
    L = load(spec{k, 2}, 'TableStimComp');                        % animals present in this struct
    present = intersect(drop, string(unique(L.TableStimComp.animal)));   % dropAnimals errors on absent names
    Tall = buildStatsTable(spec(k, :), metrics=metrics, nBootPval=nBoot, pvalSeed=0, write=false);   % all recordings
    if isempty(present)                                           % neither animal in this comparison
        Tw = Tall;                                                % identical by construction
    else
        Tw = buildStatsTable(spec(k, :), metrics=metrics, nBootPval=nBoot, pvalSeed=0, write=false, ...
            excludeAnimals=present);                              % without recordings 1 and 8
    end
    assert(height(Tall) == height(Tw) && isequal(Tall.Metric, Tw.Metric), ...
        'Sensitivity: %s rows differ between runs.', spec{k, 1}); % rows must pair one to one
    R = [R; table(string(spec{k, 1}) + repmat("", height(Tall), 1), string(Tall.Metric), ...
        paperLabel(Tall.StimA) + " vs " + paperLabel(Tall.StimB), ...   % paper labels (SB, SDGm, SDGs)
        Tall.n, Tall.nAnimals, Tall.diffMean, Tall.diffCIlo, Tall.diffCIhi, string(Tall.pFormatted), ...
        Tw.n, Tw.nAnimals, Tw.diffMean, Tw.diffCIlo, Tw.diffCIhi, string(Tw.pFormatted), ...
        repmat(strjoin(present, '+'), height(Tall), 1), ...
        'VariableNames', {'Figure','Metric','Comparison','nAll','animalsAll','effectAll','ciLoAll','ciHiAll','pAll', ...
        'nWithout','animalsWithout','effectWithout','ciLoWithout','ciHiWithout','pWithout','removed'})];   %#ok<AGROW>
end

%% Figure 5: fraction of direction-selective units vs 5% chance
load(fullfile(d.paperFigs, 'Fig5_directionTuning_units.mat'), 'U'); % per-unit tuning table (Figure 5 cache)
A = U(U.responsive & U.present & ~isnan(U.DSI), :);               % analysable units, as in Figure 5
rng(0);                                                           % same seed as Figure 5
for s = ["MB", "SDGm", "MBR"]
    a = A(A.stim == s, :);                                        % this stimulus, all recordings
    w = a(~ismember(a.exp, [1 8]), :);                            % without recordings 1 and 8
    [fa, ca, pa, na] = fracDS(a, nBoot);                          % all recordings
    [fw, cw, pw, nw] = fracDS(w, nBoot);                          % without 1 and 8
    R = [R; table("Figure 5", "Fraction direction-selective", s + " vs 5% chance", height(a), na, fa, ca(1), ca(2), ...
        pOneTailed(pa, nBoot), height(w), nw, fw, cw(1), cw(2), pOneTailed(pw, nBoot), ...
        strjoin(intersect(drop, unique(a.animal)), '+'), 'VariableNames', R.Properties.VariableNames)];   %#ok<AGROW>
end

disp(R);                                                          % print the table
writetable(R, outCsv);                                            % supplementary CSV
fprintf('Wrote %s\n', outCsv);

%% Local functions
function [f, ci, p, nAn] = fracDS(a, nBoot)
% Per-insertion fraction of direction-selective units, hierarchical bootstrap
% (animal, then insertion) and one-tailed p vs 0.05, exactly as in Figure 5.
[g, ex, an] = findgroups(a.exp, a.animal);                        % insertions
fr = splitapply(@(p) mean(p < 0.05), a.pDSI, g);                  % fraction DS per insertion
b = hierBootMatchFreq(fr, nBoot, double(categorical(an)), ex);    % (data, nrep, ANIMAL, INSERTION)
f = mean(b); ci = prctile(b, [2.5 97.5]);                         % bootstrap mean and 95% CI
p = mean(b <= 0.05);                                              % one-tailed P(fraction <= chance)
nAn = numel(unique(an));                                          % animals contributing
end

function s = pOneTailed(p, nBoot)
% One-tailed bootstrap p as text; 0 replicates reaching chance -> upper bound 1/nBoot.
if p == 0, s = string(sprintf('<%g (one-tailed)', 1 / nBoot));    % no replicate at or below chance
else,      s = string(sprintf('%.4g (one-tailed)', p)); end       % exact bootstrap proportion
end

function s = paperLabel(x)
% Paper labels for the stimulus codes stored in the structs (as makeStatsSummaryTables).
s = string(x);                                                    % default: unchanged
s(s == "RG") = "SB"; s(s == "MG") = "SDGm"; s(s == "SG") = "SDGs";   % static ball, drifting / static grating
end
