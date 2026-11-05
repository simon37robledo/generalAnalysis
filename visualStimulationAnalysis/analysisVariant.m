function v = analysisVariant(tag)
% analysisVariant  Get or set the sensitivity-analysis variant ('' = main analysis).
%
%   analysisVariant('noPR')   later stimulus objects read labels from <sorting parent>\noPR_<sorting folder>,
%                             keep their caches in <stim>_Analysis_noPR (seeded with the main sync caches),
%                             and pooled structs/figures go to Paper_figs\Sensitivity_noPR
%   analysisVariant('')       back to the main analysis (default)
%   v = analysisVariant()     current variant, '' if none
%
%   A variant never reads or writes the main labels, caches or figures, so a sensitivity
%   run cannot change published results. Set it before constructing stimulus objects
%   (the folders are fixed in VStimAnalysis.initialize) and reset it when done.

if nargin == 1                                                   % set mode
    assert(isempty(tag) || isvarname(char(tag)), 'analysisVariant: tag must be a valid name, e.g. ''noPR''.');  % used in folder names
    setappdata(groot, 'analysisVariant', char(tag));             % '' clears the variant
end
v = '';                                                          % default: main analysis
if isappdata(groot, 'analysisVariant')
    v = getappdata(groot, 'analysisVariant');                    % current variant
end
end
