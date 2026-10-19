function outPaths = makeStatsSummaryTables(params)
% makeStatsSummaryTables  Titled CSV tables for the supplementary results.
%
%   outPaths = makeStatsSummaryTables()
%   outPaths = makeStatsSummaryTables(runNull = false)
%
%   Writes one CSV per table into <Paper_figs>\Stats summary tables\. Each file
%   starts with a title row and a one-line description, then the column
%   header and data, then footnotes (definitions, tests, sample sizes), so it
%   can be pasted into the supplement as-is. Files are UTF-8 with a byte-order
%   mark (Excel then renders Δ and − correctly) and every field is quoted.
%
%   Tables
%     Fig1_MB_vs_SB_summary.csv                 group and paired-difference statistics
%     Fig1_units_per_animal.csv                 recorded / responsive units per animal
%     Fig1_delta_rate_per_animal.csv            per-animal medians of the delta rate
%     Fig1_LOAO_delta_rate.csv                  leave-one-animal-out, delta rate
%     Fig1_LOAO_zscore.csv                      leave-one-animal-out, z-score
%     Fig1_LOAO_fraction_responsive.csv         leave-one-animal-out, fraction responsive
%     Null_calibration_false_positive_rate.csv  responsiveness tests under H0
%     Null_calibration_delta_rate_bias.csv      split-half delta rate under H0
%     Null_calibration_other_stimuli.csv        FFF, SDG, MBR, NI designs under H0
%
%   Every number is read from the saved Figure 1 struct, the leave-one-animal-out
%   spreadsheet (leaveOneAnimalOut) or recomputed by checkMovingWindowNullFPR; nothing is
%   typed in by hand.
%
%   See also: AllExpAnalysis, leaveOneAnimalOut, checkMovingWindowNullFPR

arguments
    params.structPath char = 'W:\Large_scale_mapping_NP\Paper_figs\Ex_1-97_Combined_MB-RG_analysisStruct.mat'
    params.loaoPath   char = 'W:\Large_scale_mapping_NP\Paper_figs\LeaveOneAnimalOut.xlsx'  % xlsx: the CSV twin has free text with commas that readtable cannot delimit
    params.outDir     char = 'W:\Large_scale_mapping_NP\Paper_figs\Stats summary tables'
    params.runNull    logical = true      % recompute the null calibration (~1 min)
    params.nullExp    double  = 97        % recording whose trial layout the simulation uses
    params.nullNPerRate double = 240      % simulated neurons per firing rate
    params.nullBoot   double  = 2000      % sign-flip permutations in the simulation
end

assert(isfile(params.structPath), 'Figure 1 struct not found:\n  %s', params.structPath);  % fail loudly
assert(isfile(params.loaoPath), ['Leave-one-animal-out table not found:\n  %s\n', ...
    'Run leaveOneAnimalOut on the Figure 1 struct first.'], params.loaoPath);              % fail loudly
if ~isfolder(params.outDir), mkdir(params.outDir); end                                     % create output folder

S = load(params.structPath);                                       % Figure 1 analysis struct
assert(isfield(S, 'summaryStats') && isfield(S.summaryStats, 'DeltaRate') && ...
    ~isempty(S.summaryStats.DeltaRate), ...
    'Struct lacks DeltaRate statistics (statsVersion < 4); rerun the Figure 1 call.');      % outdated struct
pRes  = 2 / S.params.nBootPval;                                    % resolution of the figure p-values
disp_ = @(s) replace(string(s), "RG", "SB");                       % figure labels: RG is shown as SB (static ball)
outPaths = {};                                                     % files written

% Shared methods text for the Figure 1 footnotes
nInsUnits = numel(unique(string(S.TableStimComp.animal) + "|" + string(S.TableStimComp.insertion)));  % insertions contributing responsive units
isMBrow = string(S.TableRespNeurs.stimulus) == "MB";             % one row per insertion
mUnits = sprintf(['Units: somatic units (Kilosort4 + BombCell) present throughout both the MB and SB blocks ', ...
    '(per-block presence ratio >= 0.7 over 10 bins) and responsive to MB or SB (permutation p < %.2g; no FDR): ', ...
    '%d units from %d insertions in %d animals (recorded: %d somatic units, %d present in both blocks; %d insertions). ', ...
    'Fraction responsive uses all %d insertions, over the units present in both blocks.'], ...
    S.params.threshold, height(S.TableStimComp)/2, nInsUnits, ...
    numel(unique(string(S.TableStimComp.animal))), ...
    sum(S.TableRespNeurs.recordedSomaticN(isMBrow)), sum(S.TableRespNeurs.totalSomaticN(isMBrow)), ...
    numel(unique(string(S.TableRespNeurs.insertion))), numel(unique(string(S.TableRespNeurs.insertion))));
mTests = ['Responsiveness: max-statistic sign-flip permutation test (10000 permutations; right-tailed). ', ...
    'MB: spatial grid (9 x 9 cells); mean rate in a 300 ms window from each ball crossing; maximum over cells x directions. ', ...
    'SB: mean rate in a 300 ms ON window from onset and a 300 ms OFF window from stimulus offset (per-recording stimulus duration); ', ...
    'maximum over positions x {ON, OFF}. Maxima taken after the sign flip in both tests. Both against a 300 ms pre-stimulus baseline; ', ...
    'white stimuli only (luminosity 255).'];
mStats = sprintf(['Statistics: hierarchical bootstrap (animals -> insertions -> units; Saravanan et al. 2020); ', ...
    '95%% percentile CIs from 10000 replicates; two-tailed p-values from %g replicates ', ...
    '(p below %s reported as an upper bound).'], S.params.nBootPval, formatPValue(0, pRes));

% =========================================================================
% 1. Figure 1 summary: groups and paired differences, three metrics
% =========================================================================
metrics = {'Zscore', 'Z-score'; ...                                % summaryStats field, display name
           'DeltaRate', 'Δ firing rate (spikes/s)'; ...
           'FracResponsive', 'Fraction responsive'};
R = strings(0, 13);                                                % accumulated rows
for k = 1:size(metrics, 1)
    st   = S.summaryStats.(metrics{k,1});                          % stats block of this metric
    unit = "units"; if k == 3, unit = "insertions"; end            % what n counts
    for g = 1:height(st.groups)                                    % one row per group (MB, SB)
        G = st.groups(g, :);
        R(end+1, :) = [metrics{k,2}, disp_(G.group), G.n, unit, G.nAnimals, G.nInsertions, ...
            f(G.mean,3), f(G.median,3), f(G.ciLo,3), f(G.ciHi,3), "", "", ""];   %#ok<AGROW>
    end
    D = st.diffs(1, :);                                            % the single MB - SB pair
    [pt, po] = formatPValue(D.p, pRes);                            % exact p, or upper bound if 0
    R(end+1, :) = [metrics{k,2}, replace(disp_(D.group), "-", " − "), D.n, unit, D.nAnimals, D.nInsertions, ...
        f(D.mean,3), f(D.median,3), f(D.ciLo,3), f(D.ciHi,3), ...
        sprintf('%s [%s, %s]', f(D.effect,2), f(D.effectLo,2), f(D.effectHi,2)), ...
        sprintf('%s [%s, %s]', f(D.propPos,2), f(D.propPosLo,2), f(D.propPosHi,2)), ...
        pText(po, pt)];                                  %#ok<AGROW>
end
hdr = ["Metric","Group","n","n counts","Animals","Insertions","Mean","Median", ...
       "95% CI lower","95% CI upper","Effect size d_z [95% CI]","P(MB > SB) [95% CI]","p (MB vs SB)"];
outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, 'Fig1_MB_vs_SB_summary.csv'), ...
    'Figure 1 — Moving ball (MB) vs static ball (SB): response magnitude and fraction of responsive units', ...
    'Group rows: values per stimulus. MB − SB rows: paired within-unit (or within-insertion) differences.', ...
    hdr, R, {mUnits, mTests, ...
    ['Z-score: (observed statistic − mean of its permutation null) / SD of per-trial baseline rates. ', ...
     'Δ firing rate: split-half cross-validated peak response − baseline (preferred cell/position chosen on odd repetitions; measured on even; and vice versa).'], ...
    'Fraction responsive: responsive / somatic units per insertion. d_z: mean / SD of the paired differences.', mStats});

% =========================================================================
% 2. Units per animal
% =========================================================================
TR  = S.TableRespNeurs; TR.animal = string(TR.animal); TR.stimulus = string(TR.stimulus);  % per insertion x stimulus
TSm = S.TableStimComp(string(S.TableStimComp.stimulus) == "MB", :);  % one row per analysed unit
animals = unique(TR.animal);                                       % animals in the comparison
R = strings(0, 9);
for a = animals'
    mb = TR(TR.animal == a & TR.stimulus == "MB", :);              % MB rows of this animal
    sb = TR(TR.animal == a & TR.stimulus == "RG", :);              % SB rows of this animal
    assert(isequal(sort(string(mb.insertion)), sort(string(sb.insertion))), ...
        'Insertions differ between MB and SB rows for %s.', a);    % same insertions for both stimuli
    nSom = sum(mb.totalSomaticN);                                  % units present in both blocks (analysable)
    R(end+1, :) = [a, height(mb), sum(mb.recordedSomaticN), nSom, sum(mb.respNeur), sum(sb.respNeur), ...
        sum(string(TSm.animal) == a), f(sum(mb.respNeur)/nSom,3), f(sum(sb.respNeur)/nSom,3)];  %#ok<AGROW>
end
nAll = sum(TR.totalSomaticN(TR.stimulus=="MB"));                   % analysable units, all animals
R(end+1, :) = ["Total", height(TR)/2, sum(TR.recordedSomaticN(TR.stimulus=="MB")), nAll, ...
    sum(TR.respNeur(TR.stimulus=="MB")), sum(TR.respNeur(TR.stimulus=="RG")), height(TSm), ...
    f(sum(TR.respNeur(TR.stimulus=="MB"))/nAll,3), f(sum(TR.respNeur(TR.stimulus=="RG"))/nAll,3)];
hdr = ["Animal","Insertions","Somatic units recorded","Present in both blocks","MB-responsive","SB-responsive", ...
       "Responsive to MB or SB","Fraction MB-responsive","Fraction SB-responsive"];
outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, 'Fig1_units_per_animal.csv'), ...
    'Figure 1 — Recorded and responsive units per animal', ...
    'Pooled fractions per animal (responsive / units present in both blocks); the statistics in the summary table use per-insertion fractions.', ...
    hdr, R, {mUnits, mTests});

% =========================================================================
% 3. Delta firing rate per animal
% =========================================================================
TS = S.TableStimComp; TS.animal = string(TS.animal); TS.stimulus = string(TS.stimulus);
W  = unstack(TS(:, {'animal','realInsertion','NeurID','stimulus','DeltaRate'}), 'DeltaRate', 'stimulus');  % MB, RG per unit
assert(all(ismember({'MB','RG'}, W.Properties.VariableNames)) && ~any(isnan([W.MB; W.RG])), ...
    'Every analysed unit must have a delta rate for both MB and SB.');  % pairing invariant
R = strings(0, 6);
for a = unique(W.animal)'
    w = W(W.animal == a, :);                                       % units of this animal
    R(end+1, :) = [a, height(w), f(median(w.MB),2), f(median(w.RG),2), ...
        f(median(w.MB - w.RG),2), f(mean(w.MB - w.RG),2)];         %#ok<AGROW>
end
hdr = ["Animal","Units","Median Δ rate MB (spikes/s)","Median Δ rate SB (spikes/s)", ...
       "Median paired difference MB − SB (spikes/s)","Mean paired difference MB − SB (spikes/s)"];
outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, 'Fig1_delta_rate_per_animal.csv'), ...
    'Figure 1 — Cross-validated Δ firing rate per animal', ...
    'Descriptive per-animal summary of the units in the MB vs SB comparison (no per-animal inference).', ...
    hdr, R, {mUnits, ['Δ firing rate: split-half cross-validated peak response − baseline; ', ...
    'a unit without a consistent preferred cell/position can be negative.']});

% =========================================================================
% 4-6. Leave-one-animal-out, one file per metric
% =========================================================================
L = readtable(params.loaoPath, 'TextType', 'string');              % written by leaveOneAnimalOut
assert(all(ismember({'droppedAnimal','Metric','diffMean','pValue'}, L.Properties.VariableNames)), ...
    'Unexpected leave-one-animal-out table layout.');              % fail loudly on a wrong file
loaoRes = 2 / 1e4;                                                 % leaveOneAnimalOut default nBootPval = 1e4
loao = {'DeltaRate', 'Δ firing rate (spikes/s)', 'Fig1_LOAO_delta_rate.csv', "units"; ...
        'Z-score', 'Z-score', 'Fig1_LOAO_zscore.csv', "units"; ...
        'FracResponsive', 'Fraction responsive', 'Fig1_LOAO_fraction_responsive.csv', "insertions"};
for k = 1:size(loao, 1)
    Lk = L(L.Metric == loao{k,1}, :);                              % rows of this metric
    assert(height(Lk) == numel(unique(Lk.droppedAnimal)) && any(Lk.droppedAnimal == "none"), ...
        'LOAO %s: expected one row per dropped animal plus "none".', loao{k,1});  % one pair per run
    R = strings(height(Lk), 8);
    for i = 1:height(Lk)
        [pt, po] = formatPValue(Lk.pValue(i), loaoRes);            % exact p, or upper bound if 0
        lbl = Lk.droppedAnimal(i); if lbl == "none", lbl = "None (all animals)"; end
        R(i, :) = [lbl, Lk.n(i), Lk.nAnimals(i), f(Lk.diffMean(i),3), ...
            sprintf('[%s, %s]', f(Lk.diffCIlo(i),3), f(Lk.diffCIhi(i),3)), f(Lk.CohensD(i),2), ...
            f(Lk.PropAgtB(i),2), pText(po, pt)];
    end
    hdr = ["Animal removed", "n (" + loao{k,4} + ")", "Animals", "Mean difference MB − SB", ...
           "95% CI", "Effect size d_z", "P(MB > SB)", "p"];
    outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, loao{k,3}), ...
        sprintf('Figure 1 — Leave-one-animal-out sensitivity: %s (MB − SB)', loao{k,2}), ...
        'Each row repeats the MB vs SB comparison with one animal removed; the first row uses all animals.', ...
        hdr, R, {sprintf(['Same hierarchical bootstrap as the summary table; p-values from 10000 replicates ', ...
        '(resolution %s) so they differ slightly from the summary table by Monte Carlo error.'], ...
        formatPValue(loaoRes, loaoRes)), mUnits});                 %#ok<AGROW>
end

% =========================================================================
% 7-8. Null calibration (simulated neurons with no stimulus effect)
% =========================================================================
if params.runNull
    [F, D, info] = checkMovingWindowNullFPR(params.nullExp, nPerRate = params.nullNPerRate, ...
        nBoot = params.nullBoot);                                  % recomputed, never copied
    show = ["MB_slide", "MB_grid", "RG", "RG_onoff"];              % tests reported (legacy pooled-null variants omitted)
    name = ["MB sliding window (best 250 ms window per trial vs one baseline window) — not used", ...
            "MB spatial grid (300 ms window per crossing) — used", ...
            "SB ON window only (first 300 ms) — not used", ...
            "SB ON + OFF windows (300 ms each) — used"];
    R = strings(0, 6);
    for t = 1:numel(show)
        Ft = F(F.test == show(t), :);
        for i = 1:height(Ft)
            R(end+1, :) = [name(t), Ft.rateHz(i), Ft.n(i), Ft.nSig(i), f(Ft.FPR(i),3), ...
                sprintf('[%s, %s]', f(Ft.ciLo(i),3), f(Ft.ciHi(i),3))];  %#ok<AGROW>
        end
    end
    mSim = sprintf(['Simulation: homogeneous Bernoulli spiking (1 ms bins) at the same rate before and during the stimulus; ', ...
        'trial layout of one recording (MB: %d categories x %d trials, %d ms; SB: %d positions x %d trials); ', ...
        '%d neurons per rate in simulated recordings mixing all rates; %d sign-flip permutations; alpha = 0.05.'], ...
        info.nCatMB, info.tcMB, info.stimDurMB, info.nCatRG, info.tcRG, params.nullNPerRate, params.nullBoot);
    hdr = ["Test","Firing rate (Hz)","Simulated neurons","False positives","False-positive rate","95% CI (Clopper-Pearson)"];
    outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, 'Null_calibration_false_positive_rate.csv'), ...
        'Null calibration of the responsiveness tests', ...
        'False-positive rate at alpha = 0.05 for simulated units with no stimulus response; a valid test gives about 0.05.', ...
        hdr, R, {mSim});

    labD = containers.Map({'MB_slide','MB_grid','RG','RG_onoff'}, ...
        {'MB, 250 ms windows over the full trajectory', 'MB spatial grid (cell x direction)', ...
         'SB ON window (position)', 'SB ON + OFF windows (position x window)'});
    R = strings(height(D), 5);
    for i = 1:height(D)
        R(i, :) = [labD(char(D.layout(i))), D.rateHz(i), D.n(i), f(D.meanDeltaHz(i),3), f(D.seDeltaHz(i),3)];
    end
    hdr = ["Selection space","Firing rate (Hz)","Simulated neurons","Mean Δ rate (spikes/s)","SEM (spikes/s)"];
    outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, 'Null_calibration_delta_rate_bias.csv'), ...
        'Null calibration of the split-half Δ firing rate', ...
        'Mean cross-validated Δ rate for simulated units with no stimulus response; an unbiased estimator averages 0.', ...
        hdr, R, {mSim});

    % Designs of the other stimulus families (Figures 3-4), on their real layouts
    [F2, ~, L2] = checkLayoutNullFPR(nPerRate = params.nullNPerRate, nBoot = params.nullBoot);  % recomputed
    R = strings(height(F2), 6);
    for i = 1:height(F2)
        R(i, :) = [F2.test(i), F2.rateHz(i), F2.n(i), F2.nSig(i), f(F2.FPR(i),3), ...
            sprintf('[%s, %s]', f(F2.ciLo(i),3), f(F2.ciHi(i),3))];
    end
    lay = strjoin(arrayfun(@(x) sprintf('%s: %d conditions x %d trials, %d window(s) of %d ms vs %d ms baseline', ...
        x.name, x.nCats, x.tc, x.nWin, x.winLen, x.baseLen), L2, 'UniformOutput', false), '; ');  % layouts used
    hdr = ["Test","Firing rate (Hz)","Simulated neurons","False positives","False-positive rate","95% CI (Clopper-Pearson)"];
    outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, 'Null_calibration_other_stimuli.csv'), ...
        'Null calibration of the responsiveness tests for FFF, SB grating, drifting grating, moving bar and natural images', ...
        'False-positive rate at alpha = 0.05 for simulated units with no stimulus response; a valid test gives about 0.05.', ...
        hdr, R, {['Layouts (real recordings): ' lay], ...
        'Simulation: homogeneous Bernoulli spiking (1 ms bins); maximum over conditions x windows after the sign flip.'});
end

% =========================================================================
% 9. Leave-one-animal-out summary across every comparison (supplementary)
% =========================================================================
figDir = fileparts(params.loaoPath);                               % LOAO files live next to the structs
R = strings(0, 10);                                                % accumulated rows
% Figures 1, 3, 4: leaveOneAnimalOut spreadsheets (paired hierarchical bootstrap, two-tailed)
xl = {'LeaveOneAnimalOut.xlsx', 'LeaveOneAnimalOut_SDGm-SDGs.xlsx', 'LeaveOneAnimalOut_RG-SDGs.xlsx', ...
      'LeaveOneAnimalOut_MB-MBR.xlsx', 'LeaveOneAnimalOut_MB-SDGm.xlsx'};
for k = 1:numel(xl)
    fx = fullfile(figDir, xl{k});
    if ~isfile(fx), warning('makeStatsSummaryTables:noLOAO', 'Missing %s - skipped.', fx); continue, end   % explicit
    Lx = readtable(fx, 'TextType', 'string');                       % one row per run x metric x pair
    [G, gid] = findgroups(Lx(:, {'Figure','Metric','StimA','StimB'}));   % one group per comparison
    for g = 1:height(gid)
        Lg = Lx(G == g, :);
        full = Lg(Lg.droppedAnimal == "none", :); drops = Lg(Lg.droppedAnimal ~= "none", :);
        assert(height(full) == 1 && height(drops) >= 1, 'LOAO %s: expected one full-data row and >= 1 drop.', xl{k});
        [~, iw] = max(drops.pValue);                               % run whose removal weakens the result most
        R(end+1, :) = [full.Figure, stimLabel(gid.StimA(g)) + " vs " + stimLabel(gid.StimB(g)), gid.Metric(g), full.n, full.nAnimals, ...
            sprintf('%s [%s, %s]', f(full.diffMean,3), f(full.diffCIlo,3), f(full.diffCIhi,3)), loaoP(full.pValue, loaoRes), ...
            loaoP(min(drops.pValue), loaoRes) + " – " + loaoP(max(drops.pValue), loaoRes), drops.droppedAnimal(iw), ...
            sprintf('%d / %d', sum(drops.pValue < 0.05), height(drops))];   %#ok<AGROW>
    end
end
% Figure 2: spatial tuning index, MB vs SB (unpaired, one-tailed P(SB >= MB))
f2 = fullfile(figDir, 'LeaveOneAnimalOut_Fig2_SpatialTuning.csv');
if isfile(f2)
    L2 = readtable(f2, 'TextType', 'string'); full = L2(L2.leftOut == "none", :); drops = L2(L2.leftOut ~= "none", :);
    [~, iw] = max(drops.pOneTailed);
    R(end+1, :) = ["Figure 2", "MB vs SB (unpaired)", "Spatial tuning index", full.nMB + " / " + full.nSB, "", ...
        sprintf('%s (MB %s, SB %s)', f(full.meanMB - full.meanSB, 3), f(full.meanMB, 3), f(full.meanSB, 3)), ...
        loaoP(full.pOneTailed, loaoRes) + " (one-tailed)", loaoP(min(drops.pOneTailed), loaoRes) + " – " + loaoP(max(drops.pOneTailed), loaoRes), ...
        drops.leftOut(iw), sprintf('%d / %d', sum(drops.pOneTailed < 0.05), height(drops))];
else
    warning('makeStatsSummaryTables:noLOAO', 'Missing %s - skipped.', f2);
end
% Figure 5: fraction of direction-selective units vs 5% chance (one-tailed)
f5a = fullfile(figDir, 'Fig5_fraction_direction_selective.csv'); f5b = fullfile(figDir, 'Fig5_LeaveOneAnimalOut_fractionDS.csv');
if isfile(f5a) && isfile(f5b)
    A5 = readtable(f5a, 'TextType', 'string'); B5 = readtable(f5b, 'TextType', 'string');
    for i = 1:height(A5)
        drops = B5(B5.stim == A5.stim(i), :); [~, iw] = max(drops.pVsChance);
        R(end+1, :) = ["Figure 5", A5.stim(i) + " vs 5% chance", "Fraction direction-selective", A5.nUnits(i), A5.nAnimals(i), ...
            sprintf('%s [%s, %s]', f(A5.fracDS(i),3), f(A5.ciLo(i),3), f(A5.ciHi(i),3)), loaoP(A5.pVsChance(i), loaoRes) + " (one-tailed)", ...
            loaoP(min(drops.pVsChance), loaoRes) + " – " + loaoP(max(drops.pVsChance), loaoRes), drops.dropped(iw), ...
            sprintf('%d / %d', sum(drops.pVsChance < 0.05), height(drops))];   %#ok<AGROW>
    end
else
    warning('makeStatsSummaryTables:noLOAO', 'Missing Figure 5 tables - skipped.');
end
hdr = ["Figure", "Comparison", "Metric", "n (all animals)", "Animals", "Effect (all animals) [95% CI]", ...
       "p (all animals)", "p range, one animal removed", "Removal giving the largest p", "Runs with p < 0.05"];
outPaths{end+1} = writeTitledCsv(fullfile(params.outDir, 'Leave_one_animal_out_summary.csv'), ...
    'Leave-one-animal-out sensitivity of every cross-animal comparison', ...
    'Each comparison repeated with one animal removed at a time; the p range shows how much any single animal drives the result.', ...
    hdr, R, {['Figures 1, 3, 4: paired hierarchical bootstrap of the per-unit (or per-insertion, fraction responsive) ', ...
    'difference, two-tailed. Figure 2: two-sample hierarchical bootstrap, one-tailed P(SB >= MB). Figure 5: per-insertion ', ...
    'fraction of direction-selective units vs 0.05, one-tailed. All: animals -> insertions -> units (hierBootMatchFreq), ', ...
    '10000 replicates.']});

fprintf('Wrote %d tables to %s\n', numel(outPaths), params.outDir);  % report
end

function s = stimLabel(x)
% Paper labels for stimulus codes (RG is shown as SB; MG/SG are the drifting/static gratings).
map = struct('RG', "SB", 'MG', "SDGm", 'SG', "SDGs");                % code -> paper label
x = string(x); s = x;                                              % default: unchanged
if isfield(map, char(x)), s = map.(char(x)); end                  % exact-code mapping only
end

function s = loaoP(p, res)
% p as text: exact, or an upper bound when no replicate reached it.
[pt, po] = formatPValue(p, res); s = pText(po, pt);
end

% =========================================================================
% Local functions
% =========================================================================
function s = f(x, nDec)
% Fixed-point text with nDec decimals; empty for NaN.
if isnan(x), s = ""; else, s = string(sprintf('%.*f', nDec, x)); end
end

function s = pText(op, txt)
% Exact p as a plain number ("0.00056"); an upper bound keeps its "<" ("<0.0001").
if strcmp(op, '<'), s = "<" + string(txt); else, s = string(txt); end
end

function p = writeTitledCsv(p, title, desc, hdr, R, notes)
% Write title, description, blank line, header, rows, blank line, notes.
% UTF-8 with BOM; every field quoted (internal quotes doubled).
q   = @(c) """" + replace(string(c), """", """""") + """";       % quote one field
row = @(c) strjoin(q(c), ",");                                     % quoted CSV line
fid = fopen(p, 'w', 'n', 'UTF-8');                                 % UTF-8 text file
assert(fid > 0, 'Cannot open %s for writing.', p);                 % fail loudly
cleaner = onCleanup(@() fclose(fid));                              % close on any exit
fprintf(fid, '%s', char(65279));                                   % byte-order mark (Excel UTF-8 detection)
fprintf(fid, '%s\n', row(title));                                  % title row
fprintf(fid, '%s\n', row(desc));                                   % one-line description
fprintf(fid, '\n');                                                % spacer
fprintf(fid, '%s\n', row(hdr));                                    % column header
for i = 1:size(R, 1)
    fprintf(fid, '%s\n', row(R(i, :)));                            % data rows
end
fprintf(fid, '\n');                                                % spacer
for i = 1:numel(notes)
    fprintf(fid, '%s\n', row(notes{i}));                           % footnotes
end
end
