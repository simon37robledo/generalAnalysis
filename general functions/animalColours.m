function cmap = animalColours(names)
% animalColours  Fixed colour per animal, identical in every figure and comparison.
%
%   cmap = animalColours(names)
%
%   names - cellstr / string / categorical of animal IDs (e.g. categories(tbl.animal)).
%   cmap  - numel(names) x 3 RGB matrix, row k = colour of names(k).
%
%   Before 2026-11-17 each panel coloured animals by lines(nAnimals) over the animals present
%   in that panel, so the same animal changed colour between comparisons (e.g. SA9 yellow in
%   MB vs SB, orange-red in the grating panels), and lines(8) repeats its first colour (two
%   Pogona animals were both blue). Colours now come from this fixed registry.
%   Pogona and Laudakia are never pooled in one comparison, so the two species reuse the
%   same palette slots; a mixed call warns.

names = cellstr(string(names));                                     % accept cellstr, string or categorical
registry = {'PV103','PV132','PV140','PV233','PV27','PV35','PV51','PV97', ...   % Pogona, alphabetical (keeps the old Fig 1 colours)
            'SA5','SA8','SA9','SA11'};                                         % Laudakia (keeps the old Fig 7 MB-SB colours)
slot     = [1:8, 1:4];                                              % palette slot of each registry entry (species reuse slots)
palette  = [lines(7); 0.90 0.45 0.70];                              % 8 distinct colours: MATLAB's 7 + pink (lines(8) would repeat blue)
unknown  = [0.5 0.5 0.5];                                           % grey for an animal missing from the registry

cmap = zeros(numel(names), 3);                                      % output, one row per name
for k = 1:numel(names)
    idx = find(strcmp(registry, names{k}));                         % position in the registry
    if isscalar(idx)
        cmap(k, :) = palette(slot(idx), :);                         % fixed colour of this animal
    else
        cmap(k, :) = unknown;                                       % not registered: grey, never a silent reuse
        warning('animalColours:unknownAnimal', ...
            'Animal "%s" is not in the animalColours registry: drawn grey. Add it to the registry.', names{k});
    end
end

if any(startsWith(names, 'PV')) && any(startsWith(names, 'SA'))    % both species in one call
    warning('animalColours:mixedSpecies', ...
        'Pogona and Laudakia animals in one call: they share palette slots, so colours repeat.');
end
end
