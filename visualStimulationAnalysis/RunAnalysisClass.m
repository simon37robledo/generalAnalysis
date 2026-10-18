cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);



%% AllExpAnalysis
%[49:54 57:81] MBR all experiments 'NV','NI'
%[44:56,64:88] All experiments
%[28:32,44,45,47,48,56,98] All SA experiments
%Check triggers 45, SA82 44,45,47:54,56,64:88 
% All stim: 'FFF','SDG','MBR','MB','RG','NI','NV'
%[49:54,64:97] %All PV good experiments [49:54,64:85 87:97]
% %%[89,90,92,93,95,96,97]  %Al NV and NI experiments 
%[49:54,84:90,92:96] %All SDG experiments
%solve MBR
%bootsrapRespBase

%% FIGURE 1 (MB vs SB) and the dataset summary (responsive union) now live in
% paperFigures/Fig1_MBvsRG.m and paperFigures/DatasetSummary_ResponsiveUnion.m;
% every figure takes its statistics design from paperFigures/figureDesign.m.

%% FIGURE 2: SPATIAL TUNING ANS STRIPE ANALYSIS
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%%  Calculate spatial tuning
for tp = 5
results = SpatialTuningIndex([40:43,49:54,64:66, 68:85 87:97], indexType =  "L_amplitude_diff" ,overwrite=true...
    , topPercent = tp,useRF=true,onOff=1,unionResponsive = false,allResponsive=true, PaperFig=true, plotRFs=false, plotRFunion=false, detectStripe = true, testStripeDetection= true);

end

%%
results = analyzeStripeNeurons([40:43,49:54,64:66, 68:85 87:97], ...
    indexType   = "L_amplitude_diff", ...
    onOff       = 1, ...
    sizeIdx     = 1, ...
    lumIdx      = 1, ...
    useRF       = true, ...
    prefDir     = true, ...
    unionResponsive = false,...
    allResponsive=true, PaperFig=true);   % match what you used in SpatialTuningIndex

%% Plot examples neurons MB and SB
%%% Rect Grid

for ex = [81] %84:91
    NP = loadNPclassFromTable(ex); %73 81
    vsRe = rectGridAnalysis(NP);
    vsRe.plotRaster(MergeNtrials=1,overwrite=true,AllResponsiveNeurons=false,exNeuronsPhyID = [288], exNeurons= [], selectedLum=255,oneTrial = true,PaperFig = true) 
    [colorbarLimsRG] = vsRe.PlotReceptiveFields(exNeurons=14,allStimParamsCombined=false,PaperFig=true,overwrite=true);
end

%%% Moving ball

for ex =[70]%
    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBallAnalysis(NP);
    vs.plotRaster('exNeuronsPhyID',[180],'overwrite',true,'MergeNtrials',2,'PaperFig',true,OneLuminosity= 'white',OneDirection='down',bin=20)
    colorbarLims=vs.PlotReceptiveFields('exNeurons',19,'overwrite',true,'OneDirection','down','OneLuminosity','white','PaperFig',true);
end

%% FIGURE 3 MOVING VS STATIC COMPARISON
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%% Figure 3 main comparison (SDGm vs SDGs), PSTH and raster: see
% paperFigures/Fig3_MovingVsStatic.m. Category cells below are pending revalidation.
% DO NOT RUN them before the category round: their BaseRespWindow / maxCategory
% differ from figureDesign.m, so they would overwrite the SDG and MB statistics
% caches that Figures 1, 3 and 4 and the union are built from.

%% %% Compare SDGm vs SDGs, across directions
[tempTableMW] = AllExpAnalysis([40:43 49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'SDGm','SDGs'},CompareCategory={"angles","angles"},CompareLevels={[0,90,180,270],[0,90,180,270]},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true,useFDR=false,maxCategory=false,BaseRespWindow=1000);


%% SDGm spatial frequency
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'SDGm'},CompareCategory={"tempFrequency"},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 1000);


%% NI / NV: NV is not in the paper; the NI lizard-posture comparison (imgeOrder) is
% rerun in the category round (StatisticsPerNeuronPerCategory fixes pending).

%% FIGURE 4 SIZES AND LOCALITY COMPARISON
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
 
% %% Compares MB across different sizes
% [tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=true,ComparePairs={'MB'},CompareCategory="sizes",PaperFig=true,...
%     overwriteResponse=false,overwriteStats=true, BaseRespWindow = 500,CategoryMaximized = "offsets");
% 
% 
% %% Compares MB and SDGm across different categories, z-scores are computed with a moving window and responsive units are selected on the per category p value.
% [tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'MB','SDGm'},CompareCategory={"directions","angles"},CompareLevels={[0,1.57,3.14, 4.71],[0,90,180,270]},PaperFig=true,...
%     overwriteResponse=false,overwriteStats=true, BaseRespWindow = 1000,useGeneralFilter=true,SpatialGridMode = false,maxCategory = true);
%% Figure 4 main comparisons (RG-SDGs-FFF, MB-MBR, MB-SDGm, MBR-SDGm in PV97), their
% PSTHs and rasters: see paperFigures/Fig4_Locality.m (shared design: figureDesign.m).


%% SDGm spatial frequency
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'SDGs'},CompareCategory="spatFrequency",PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 1000);

%% Plot different spatial freuqncies
plotPSTH_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, zScore=true,TakeTopPercentTrials=[], PaperFig=true, byDepth=false, smooth=50, postStim= 2000, stimTypes={"SDGm"},splitBy="spatFrequency",useCategoryPvals = true); %stimTypes=["linearlyMovingBall"]


%% Plot MB raster sorted per size
plotRaster_MultiExp([49:54,64:66,68:85 87:97], overwrite=false, sortBy="preferredCategory", ...
    splitCategory="Sizes", stimTypes="MB",zScore=true,PaperFig=true)

%% Plot SDGm raster sorted per size
plotRaster_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, sortBy="preferredCategory", ...
    splitCategory="spatFrequency", stimTypes="SDGm",zScore=true,PaperFig=true)

%% %% FIGURE 5 DIRECTION TUNING COMPARISON
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%%  Compares MB across differen directions
[tempTableMW] = AllExpAnalysis([40:43, 49:54,64:66,68:85 87:97], overwrite=true,ComparePairs={'MB'},CompareCategory="directions",PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 1500);
%%
plotPSTH_MultiExp([40:43, 49:54,64:66,68:85 87:97], overwrite=false, zScore=true,TakeTopPercentTrials=[], PaperFig=true,postStim =2000, byDepth=false, smooth=50, stimTypes={"MB"},splitBy="directions"); %stimTypes=["linearlyMovingBall"]

%%  Compares MB across differen directions
[tempTableMW] = AllExpAnalysis([40:43, 49:54,64:66,68:85 87:97], overwrite=true,ComparePairs={'SDGm'},CompareCategory="angles",CompareLevels = {[0,90,180,270]},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 1000);
%%
plotPSTH_MultiExp([40:43, 49:54,64:66,68:85 87:97], overwrite=true, zScore=true,TakeTopPercentTrials=[], PaperFig=true,postStim =1000, byDepth=false, smooth=50, stimTypes={"SDGm"},splitBy="angles",splitLevels=[0,90,180,270]); %stimTypes=["linearlyMovingBall"]

%% PLot OSI and DSI

AllExpDirectionTuning([40:43, 49:54,64:66,68:85 87:97], overwrite=true, stimuli = {'MB'},PaperFig=true)


%% Plot SDGm raster sorted per direction
plotRaster_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, sortBy="preferredCategory", ...
    splitCategory="angles",   splitLevels ={[0,90,180,270]}, stimTypes="SDGm",zScore=true,PaperFig=true)

%% Plot SDGm raster sorted per direction
plotRaster_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, sortBy="preferredCategory", ...
    splitCategory="directions", stimTypes="MB",zScore=true,PaperFig=true)

%% END OF FIGURE CALLS
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%%
%% Rect Grid
j = 1;
for ex = [69] %84:91
    NP = loadNPclassFromTable(ex); %73 81
    vsRe = rectGridAnalysis(NP);
    phyIds = T.PhyID(j);
    % vsRe.getSessionTime("overwrite",true);
    % %vsRe.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    % vsRe.getDiodeTriggers('overwrite',true);
    % vsRe.getSyncedDiodeTriggers("overwrite",true);
    % % vsRe.plotSpatialTuningSpikes;
    % % vsRe.plotSpatialTuningLFP;
    % vsRe.ResponseWindow('overwrite',true)
    % results = vsRe.ShufflingAnalysis('overwrite',true);
    %vsRe.plotRaster(MergeNtrials=1,overwrite=true,AllResponsiveNeurons = true, selectedLum=[],oneTrial = true,PaperFig = true) %43
    vsRe.plotRaster(MergeNtrials=1,overwrite=true,AllResponsiveNeurons=false,exNeuronsPhyID = [], exNeurons= [21], selectedLum=255,oneTrial = true,PaperFig = true) %43
    vsRe.CalculateReceptiveFields('overwrite',true,AllResponsiveNeurons=false)
    % % [colorbarLimsRG] = vsRe.PlotReceptiveFields(exNeurons=13,allStimParamsCombined=false,PaperFig=true,overwrite=true);
    % %result = vsRe.BootstrapPerNeuron('overwrite',true);
    % result = vsRe.StatisticsPerNeuron;

    j = j+1;

end
% vsRe.CalculateReceptiveFields
%vsRe.PlotReceptiveFields("exNeurons",18)

%% Moving ball
j =1;
for ex =[97]% 
    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBallAnalysis(NP);
    phyIds = T.PhyID(j);
    % vs.getSessionTime("overwrite",true);
    % vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    % % % %vs.plotDiodeTriggers
    % vs.getSyncedDiodeTriggers("overwrite",true);
    % % % % %vs.plotSpatialTuningSpikes;
    % r = vs.ResponseWindow('overwrite',true);
    % % % results = vs.ShufflingAnalysis('overwrite',true);
    % % % % vs.plotRaster('AllSomaticNeurons',true,'overwrite',true,'MergeNtrials',3)
    % vs.plotRaster('AllResponsiveNeurons',true,'overwrite',true,'MergeNtrials',1,'bin',50,'GaussianLength',30,'MaxVal_1', false, oneLuminosity = "white", OneDirection="left", ...
    %     sortingOrder=["size","direction","luminosity","offset","speed"])
    %vs.plotRaster('AllSomaticNeurons',true,'overwrite',true,'MergeNtrials',3,PaperFig=true)
    vs.plotRaster('exNeuronsPhyID',357,'overwrite',true,'MergeNtrials',1,'PaperFig',true,OneLuminosity= 'white',OneDirection='up',bin=20)
    % % % % %vs.plotCorrSpikePattern
    % vs.plotRaster('exNeurons',82,'overwrite',true,'OneDirection','up','OneLuminosity','white','MergeNtrials',1,'PaperFig',false)
    % %vs.plotRaster('AllSomaticNeurons',true,'overwrite',true,'MergeNtrials',3,MaxVal_1=false,PaperFig=false)
    %vs.CalculateReceptiveFields('overwrite',true,testConvolution=false,AllResponsiveNeurons=false);   
    colorbarLims=vs.PlotReceptiveFields('exNeurons',19,'overwrite',true,'OneDirection','up','OneLuminosity','white','PaperFig',true);
    % %result = vs.BootstrapPerNeuron('overwrite',true);%('overwrite',true);
    % % pvals0_6Filter =result.Speed2.pvalsResponse';
    % % compare = [pvals,pvalsNoFilt,pvals0_6Filter];
    % result = vs.StatisticsPerNeuron('overwrite',true);
    % result = vs.StatisticsPerNeuronPerCategory('compareCategory','sizes','overwrite',true);
    j = j+1;
end


%% Raster for all experiment
plotRaster_MultiExp([49:54,64:66,68:85 87:97], sortBy = "spatialTuning",overwrite=true,TakeTopPercentTrials=[],PaperFig=true)

%% Calculate spatial tuning
results= SpatialTuningIndex([49:54,64:66, 68:85 87:97], indexType =  "L_amplitude_diff" ,overwrite=true...
    , topPercent = 30,useRF=true,onOff=1,unionResponsive = false,allResponsive=true, PaperFig=true, plotRFs=false, plotRFunion=false);

%% Get neuron depths
getNeuronDepths([49:54,64:97]) %[49:54,64:72,84:97] %% PV140 missing depth coordinates
%% Gratings

for ex = [40]
    NP = loadNPclassFromTable(ex); %73 81
    vs = StaticDriftingGratingAnalysis(NP);
    vs.getSessionTime("overwrite",true);
    vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    dT = vs.getDiodeTriggers;
    % vs.plotDiodeTriggers
    vs.getSyncedDiodeTriggers("overwrite",true);
    %r = vs.ResponseWindow('overwrite',true);
    % results = vs.ShufflingAnalysis('overwrite',true);
    % result = vs.BootstrapPerNeuron('overwrite',true);
    % vs.StatisticsPerNeuron(overwrite=true)
    vs.plotRaster(MaxVal_1=true,OneAngle=270,exNeurons=28,AllResponsiveNeurons=false,PaperFig=false) %0.5208 %2.0833
    vs.plotRaster(MaxVal_1=false)
    close all
end

%% movie

for ex =  [92:97]
    NP = loadNPclassFromTable(ex); %73 81
    vs = movieAnalysis(NP);
    % vs.getSessionTime("overwrite",true);
    % vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    % dT = vs.getDiodeTriggers;
    % vs.plotDiodeTriggers
    %vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    %results = vs.ShufflingAnalysis('overwrite',true);
    result = vs.StatisticsPerNeuron('overwrite',true);
    vs.plotRaster(AllResponsiveNeurons=true)
    close all
end


%% image

for ex = [97]
    NP = loadNPclassFromTable(ex); %73 81
    vs = imageAnalysis(NP);
    %vs.getSessionTime("overwrite",true);
    %vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %dT = vs.getDiodeTriggers;
    % vs.plotDiodeTriggers
    %vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    %results = vs.ShufflingAnalysis('overwrite',true);
    vs.plotRaster('exNeurons',13,MergeNtrials=1,overwrite=true, selectCats =[], PaperFig=true)
    close all
    %result = vs.StatisticsPerNeuron('overwrite',true);

end


%% Moving bar
for ex = 81
    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBarAnalysis(NP);
    vs.getSessionTime("overwrite",true);
    vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %vs.plotDiodeTriggers
    vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    results = vs.ShufflingAnalysis('overwrite',true);
    if ~any(results.Speed1.pvalsResponse<0.05)
         fprintf('%d-No responsive neurons.\n',ex)
         continue
    end
    vs.CalculateReceptiveFields('overwrite',true,'nShuffle',20);
    vs.PlotReceptiveFields('overwrite',true,'RFsDivision',{'Directions','','Luminosities'},meanAllNeurons=true)
end

%% FFF
for ex = 56
    NP = loadNPclassFromTable(ex); %73 81
    vs = fullFieldFlashAnalysis(NP);
    vs.getSessionTime("overwrite",true);
    vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %vs.plotDiodeTriggers
    vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    results = vs.ShufflingAnalysis('overwrite',true);
    if ~any(results.Speed1.pvalsResponse<0.05)
         fprintf('%d-No responsive neurons.\n',ex)
         continue
    end
    vs.CalculateReceptiveFields('overwrite',true,'nShuffle',20);
    vs.PlotReceptiveFields('overwrite',true,'RFsDivision',{'Directions','','Luminosities'},meanAllNeurons=true)
end


%% Run for all
for ex = 85:88
    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBallAnalysis(NP);
    vs.getSessionTime("overwrite",true);
    vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %vs.plotDiodeTriggers
    vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    results = vs.ShufflingAnalysis('overwrite',true);
    if ~any(results.Speed1.pvalsResponse<0.05)
         fprintf('%d-No responsive neurons.\n',ex)
         continue
    end
    vs.CalculateReceptiveFields('overwrite',true,'nShuffle',20);
    vs.PlotReceptiveFields('overwrite',true,'RFsDivision',{'Directions','','Luminosities'},meanAllNeurons=true)
end

%% Check experiments in timseseries viewer
timeSeriesViewer(NP)
t=NP.getTrigger;
data.VS_ordered(ex)

stimOn = t{3};
stimOff = t{4};

MBRtOn = stimOn(stimOn > t{1}(1) & stimOn < t{2}(1));
MBRtOff = stimOff(stimOff > t{1}(1) & stimOff < t{2}(1));

MBtOn = stimOn(stimOn > t{1}(2) & stimOn < t{2}(2));
MBtOff = stimOff(stimOff > t{1}(2) & stimOff < t{2}(2));

RGtOn = stimOn(stimOn > t{1}(3) & stimOn < t{2}(3));
RGtOff = stimOff(stimOff > t{1}(3) & stimOff < t{2}(3));

NGtOn = stimOn(stimOn > t{1}(4) & stimOn < t{2}(4));
NGtOff = stimOff(stimOff > t{1}(4) & stimOff < t{2}(4));

DtOn = stimOn(stimOn > t{1}(5) & stimOn < t{2}(5));
DtOff = stimOff(stimOff > t{1}(5) & stimOff < t{2}(5));

MovingBallTriggersDiode = d3.stimOnFlipTimes;



%% %% check neural data sync and analog data sync 

allTimes = [stimOn(:); stimOff(:); onSync(:); offSync(:)];  % concatenate as column

% Sort from earliest to latest
sortedTimesDiodeOldMethod = sort(allTimes);
