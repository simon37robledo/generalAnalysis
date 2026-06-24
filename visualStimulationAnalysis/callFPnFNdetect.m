D = findFalseNegAndPos('W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_1-97_n44_Raster_MB-RG_union.mat');

Dfn = sortrows(D(D.FN_for~="",:), 'Severity', 'descend');   % strongest missed responses first
plotFalseNegPosNeurons(Dfn(1:15,:));     % test batch — check the MB call and the filenames
% plotFalseNegPosNeurons(D);          % full run once the test looks right