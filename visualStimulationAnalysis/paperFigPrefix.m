function p = paperFigPrefix(tag)
% paperFigPrefix  Get or set the figure-number prefix of saved paper panels.
%
%   paperFigPrefix('fig1')   later PaperFig saves are named 'fig1_<name>'
%   paperFigPrefix('')       no prefix (ad-hoc plots)
%   p = paperFigPrefix()     current prefix, '' if none
%
%   Each paperFigures/FigN_*.m script sets its own prefix at the top, so every
%   panel it saves (printFig, exportgraphics) carries the figure it belongs to.
%   Data files (analysis structs, CSV tables) are not renamed: other code looks
%   them up by name.

if nargin == 1                                                   % set mode
    if isempty(tag)
        setappdata(groot, 'paperFigPrefix', '');                 % clear the prefix
    else
        setappdata(groot, 'paperFigPrefix', [char(tag) '_']);    % e.g. 'fig1_'
    end
end
p = '';                                                          % default: no prefix
if isappdata(groot, 'paperFigPrefix')
    p = getappdata(groot, 'paperFigPrefix');                     % current prefix
end
end
