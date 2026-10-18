function pres = statPresence(Stats, stimType, fieldName)
% Unit presence during this stimulus block (StatisticsPerNeuron 'present'):
% MB/MBR must be present during every speed's block, combined SDG during both
% phases; errors if the cache predates presence (statsVersion < 9).
if ismember(stimType, ["MB","MBR"])
    f = fieldnames(Stats); f = f(startsWith(string(f), "Speed"));  % all speed sub-fields
elseif stimType == "SDG"
    f = {'Static', 'Moving'};                                        % both grating phases
elseif isfield(Stats, char(fieldName))
    f = {char(fieldName)};                                           % stimulus sub-field
else
    f = {};                                                          % flat struct (RG, FFF, NI, ...)
end
if isempty(f)
    assert(isfield(Stats, 'present'), ...
        'Stats lack unit presence (statsVersion < 9); rerun AllExpAnalysis with overwriteStats=true.');
    pres = logical(Stats.present(:));                                % flat presence vector
    return
end
pres = true;                                                         % AND over sub-fields
for k = 1:numel(f)
    assert(isfield(Stats.(f{k}), 'present'), ...
        'Stats.%s lacks unit presence (statsVersion < 9); rerun AllExpAnalysis with overwriteStats=true.', f{k});
    pres = pres & logical(Stats.(f{k}).present(:));                  % present in this block too
end
end
