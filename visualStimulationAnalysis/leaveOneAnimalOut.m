function [T, outPaths] = leaveOneAnimalOut(structPath, params)
% leaveOneAnimalOut  Sensitivity of a comparison to each single animal.
%
%   [T, outPaths] = leaveOneAnimalOut(structPath)
%   [T, outPaths] = leaveOneAnimalOut(structPath, metrics = {'DeltaRate'}, ...)
%
%   Recomputes the supplementary statistics (mean difference, hierarchical
%   bootstrap 95% CI, Cohen's d, P(A>B), exact p) once with every animal and
%   once with each animal removed. Every row comes from buildStatsTable (via
%   its excludeAnimals option), i.e. from the same routines as the figures, so
%   the "all animals" row equals the published one up to nBootPval resolution.
%
%   INPUTS
%     structPath  *_analysisStruct.mat written by AllExpAnalysis (PaperFig=true)
%     params      Name-value:
%       figLabel  - label for the Figure column (default 'Figure 1')
%       pairs     - "all" or an Mx2 cell of stimulus pairs (default "all")
%       metrics   - metrics to test (default Z-score, DeltaRate, FracResponsive;
%                   SpkR is left out because z_mean is inflated by the number of
%                   categories a stimulus offers)
%       nBootPval - replicates per p-value (default 1e4: p resolution 2e-4,
%                   sufficient for a sensitivity analysis and 10x faster)
%       pvalSeed  - RNG seed applied before every p-value (default 0)
%       write     - write .xlsx/.csv next to the struct (default true)
%       outName   - base filename (default 'LeaveOneAnimalOut')
%
%   OUTPUT
%     T         one row per (dropped animal x metric x pair); droppedAnimal
%               "none" is the full data set
%     outPaths  files written
%
%   See also: buildStatsTable, makeSupplementaryStatsTable

arguments
    structPath   (1,:) char
    params.figLabel  char   = 'Figure 1'
    params.pairs            = "all"
    params.metrics   cell   = {'Z-score','DeltaRate','FracResponsive'}
    params.nBootPval double = 1e4
    params.pvalSeed  double = 0
    params.write     logical = true
    params.outName   char   = 'LeaveOneAnimalOut'
end

assert(isfile(structPath), 'leaveOneAnimalOut: struct not found:\n  %s', structPath);  % fail loudly
S       = load(structPath, 'TableStimComp');                     % only the animal list is needed here
animals = unique(string(S.TableStimComp.animal))';               % every animal in the comparison
drops   = ["", animals];                                         % "" = keep all animals, then drop each in turn
spec    = {params.figLabel, structPath, params.pairs};           % one comparison, buildStatsTable format

T = table();                                                     % accumulator
for k = 1:numel(drops)
    excl = drops(k);                                             % animal removed in this run
    if excl == "", excl = strings(1,0); end                      % empty list = full data set
    fprintf('\n=== leave-one-animal-out %d/%d: dropping %s ===\n', k, numel(drops), ...
        string(strjoin(excl, ',')));                             % progress
    Tk = buildStatsTable(spec, metrics = params.metrics, nBootPval = params.nBootPval, ...
        pvalSeed = params.pvalSeed, write = false, excludeAnimals = excl);  % same routines as the figures
    Tk.droppedAnimal = repmat(drops(k), height(Tk), 1);          % which animal was left out
    T = [T; Tk];                                                 %#ok<AGROW>
end
T.droppedAnimal(T.droppedAnimal == "") = "none";                 % readable label for the full data set
T = movevars(T, 'droppedAnimal', 'Before', 1);                   % key column first

% ---- Compact console summary --------------------------------------------
disp(T(:, {'droppedAnimal','Metric','StimA','StimB','n','nAnimals', ...
    'diffMean','diffCIlo','diffCIhi','CohensD','pFormatted'}));   % one line per run x metric x pair

% ---- Write ----------------------------------------------------------------
outPaths = {};
if params.write
    outDir   = fileparts(structPath);                            % next to the analysed struct
    xlsxPath = fullfile(outDir, [params.outName '.xlsx']);
    csvPath  = fullfile(outDir, [params.outName '.csv']);
    writetable(T, xlsxPath);                                     % spreadsheet for the supplement
    writetable(T, csvPath);                                      % plain-text copy
    outPaths = {xlsxPath, csvPath};
    fprintf('\nWrote %d rows to:\n  %s\n  %s\n', height(T), xlsxPath, csvPath);
end
end
