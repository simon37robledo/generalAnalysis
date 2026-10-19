%% FIGURE 2 — Spatial tuning of moving-ball (MB) vs static-ball (SB, rect grid) responses
% Units: responsive (permutation p < 0.05) AND present in both the MB and SB blocks, read
% from the Figure 1 statistics caches (run Fig1_MBvsRG.m first). MB uses, per unit,
% the speed with the lowest p, as in Figure 1.
% Index: L_amplitude_diff on 9x9 grid maps of the receptive field, split half
% (MB preferred direction chosen on odd repetitions and the index measured on
% even, and vice versa; SB averaged over the same halves).
% Stripes: anisotropy on the 9x9 stimulus grid vs a pixel-permutation null; MB
% takes the max over directions in the observed score and in every surrogate.
% Group comparisons in analyzeStripeNeurons that select on the same RFs are
% descriptive (CIs only, no p-value).
% Run the cells in order; rasters must be drawn in an interactive session.

d    = figureDesign();                                            % shared design (recordings, output folder)
paperFigPrefix('fig2');                                           % saved panels are named fig2_<panel>
exps = intersect(d.exps, [40:43, 49:54, 64:66, 68:85, 87:97]);    % recordings with MB receptive-field maps
tuneArgs = {'indexType', "L_amplitude_diff", 'topPercent', 5, 'useRF', true, 'prefDir', true, ...
    'onOff', 1, 'sizeIdx', 1, 'lumIdx', 1, 'unionResponsive', false, 'allResponsive', true, ...
    'splitHalf', true};                                           % index design (shared by both calls)

%% Receptive fields on each half of the repetitions (slow: run once, cached per recording)
% Computes and caches ...-odd/-even RF files (and speed 2 where a unit's Figure-1
% speed is 2). SpatialTuningIndex would do this on demand; this cell only makes
% the long computation explicit and resumable.
for ex = exps
    NP = loadNPclassFromTable(ex);                                % recording
    for half = ["all", "odd", "even"]
        rectGridAnalysis(NP).CalculateReceptiveFields(AllResponsiveNeurons=false, trialHalf=half);   % SB maps
        mb = linearlyMovingBallAnalysis(NP);                      % MB analysis object
        for sp = 1:numel(unique(mb.VST.speed))                    % every recorded speed
            mb.CalculateReceptiveFields(AllResponsiveNeurons=false, trialHalf=half, speed=sp);       % MB maps
        end
    end
end

%% Spatial tuning index and stripe detection (swarm MB vs SB, hierarchical bootstrap)
results = SpatialTuningIndex(exps, tuneArgs{:}, overwrite=true, PaperFig=true, ...
    plotRFs=false, plotRFunion=false, detectStripe=true, testStripeDetection=false);   % cache name does not encode every param: always overwrite

%% Stripe units vs responsive non-stripe units (depth, z, OSI/DSI, direction)
stripeArgs ={'indexType', "L_amplitude_diff", 'onOff', 1, 'sizeIdx', 1, 'lumIdx', 1, 'useRF', true, ...
    'prefDir', true, 'unionResponsive', false, 'allResponsive', true, 'splitHalf', true};   % must match the call above
resStripe = analyzeStripeNeurons(exps, stripeArgs{:}, PaperFig=true);

%% Leave-one-animal-out: MB vs SB index (two-sample hierarchical bootstrap, as in SpatialTuningIndex)
S   = results;                                                    % SpatialTuningIndex output
tbl = S.tbl(S.tbl.onOff == 1 & S.tbl.sizeIdx == 1 & S.tbl.lumIdx == 1, :);   % plotted condition
animals = categories(removecats(tbl.animal));                     % animals in the figure
loao = table('Size', [numel(animals) + 1, 6], 'VariableTypes', {'string','double','double','double','double','double'}, ...
    'VariableNames', {'leftOut','nMB','nSB','meanMB','meanSB','pOneTailed'});   % one row per run
for k = 0:numel(animals)
    keep = true(height(tbl), 1);                                  % k = 0: all animals
    if k > 0, keep = tbl.animal ~= animals{k}; end                % drop one animal
    mb = tbl(keep & tbl.stimulus == "linearlyMovingBall", :);     % MB units
    sb = tbl(keep & tbl.stimulus == "rectGrid", :);               % SB units
    bMB = hierBootMatchFreq(mb.L_amplitude_diff, 10000, double(mb.animal), double(mb.experimentNum));   % animal, then insertion
    bSB = hierBootMatchFreq(sb.L_amplitude_diff, 10000, double(sb.animal), double(sb.experimentNum));
    name = "none"; if k > 0, name = string(animals{k}); end       % label of the left-out animal
    loao(k + 1, :) = {name, height(mb), height(sb), mean(mb.L_amplitude_diff), mean(sb.L_amplitude_diff), ...
        mean(bSB >= bMB)};                                        % P(SB >= MB), same convention as SpatialTuningIndex
end
disp(loao);
writetable(loao, fullfile(d.paperFigs, 'LeaveOneAnimalOut_Fig2_SpatialTuning.csv'));   % supplementary table

%% Example neurons: rasters and receptive fields (interactive session)
NP = loadNPclassFromTable(81);                                    % SB example recording
vsRe = rectGridAnalysis(NP);
vsRe.plotRaster(MergeNtrials=1, overwrite=true, AllResponsiveNeurons=false, exNeuronsPhyID=288, exNeurons=[], selectedLum=255, oneTrial=true, PaperFig=true)
colorbarLimsRG = vsRe.PlotReceptiveFields(exNeurons=14, allStimParamsCombined=false, PaperFig=true, overwrite=true);

NP = loadNPclassFromTable(70);                                    % MB example recording
vs = linearlyMovingBallAnalysis(NP);
vs.plotRaster('exNeuronsPhyID', 180, 'overwrite', true, 'MergeNtrials', 2, 'PaperFig', true, OneLuminosity='white', OneDirection='down', bin=20)
colorbarLims = vs.PlotReceptiveFields('exNeurons', 19, 'overwrite', true, 'OneDirection', 'down', 'OneLuminosity', 'white', 'PaperFig', true);
