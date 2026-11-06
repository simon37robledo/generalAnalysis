%% STIMULUS GEOMETRY IN DEGREES — pixel values of every stimulus converted to visual angle
% Flat 33 x 33 cm, 1080 x 1080 px stimulus field, 21.5 cm from the eye, eye facing the screen
% centre perpendicularly. Main values are AT THE SCREEN CENTRE; on a flat screen the same pixel
% distance subtends less in the periphery, so the value at the horizontal screen edge is also given.
%   positions (where):  atan(d/D)                      (pixels2eyePosition)
%   sizes (how big):    2*atan(s/(2D)), centred         (pixels2eyeDegrees)
%   speeds, spatial frequencies: local scale at the centre, dtheta/dx = 1/D rad per cm
% Output: <Paper_figs>\Stats summary tables\Stimulus_geometry_degrees.csv

d = figureDesign();                                                 % shared design (output folder)
outCsv = fullfile(d.paperFigs, 'Stats summary tables', 'Stimulus_geometry_degrees.csv');   % output table

D = 21.5;                                                           % eye-to-screen distance (cm)
px = 33 / 1080;                                                     % pixel side (cm): 33 cm over 1080 px
halfW = 1080 / 2 * px;                                              % centre-to-edge distance (cm) = 16.5

posDeg  = @(x) atand(x / D);                                        % eccentricity (deg) of a point x cm from the centre
sizeC   = @(s) 2 * atand(s / (2 * D));                              % size (deg) of an object s cm wide, centred
sizeE   = @(s) atand(halfW / D) - atand((halfW - s) / D);           % size (deg) of an object s cm wide whose outer side touches the screen edge
scaleC  = 180 / pi / D;                                             % local scale at the centre (deg per cm), d(atan(x/D))/dx at x = 0
scaleE  = 180 / pi * D / (D^2 + halfW^2);                           % local radial scale at the horizontal edge (deg per cm)

assert(abs(halfW - 16.5) < 1e-9, 'Screen half-width should be 16.5 cm');   % sanity check on the geometry constants

% Each row: {quantity, value in pixel units, pixel unit, kind}
spec = { ...
    'Screen extent (edge to edge)',          1080,   'px',     'extent'; ...
    'Screen corner eccentricity',            540*sqrt(2), 'px', 'position'; ...
    'One pixel',                             1,      'px',     'size'; ...
    'MB ball diameter',                      60,     'px',     'size'; ...
    'MB ball diameter (recordings 1 and 8)', 100,    'px',     'size'; ...
    'MB ball diameter (9x9 protocol)',       120,    'px',     'size'; ...
    'MB ball diameter',                      240,    'px',     'size'; ...
    'MB trajectory spacing / grid cell',     120,    'px',     'size'; ...
    'MB outermost trajectory offset',        480,    'px',     'position'; ...
    'SB circle diameter',                    119,    'px',     'size'; ...
    'MBR bar width',                         10,     'px',     'size'; ...
    'MBR bar width',                         120,    'px',     'size'; ...
    'MBR bar length',                        1024,   'px',     'size'; ...
    'MB/MBR speed',                          125,    'px/s',   'speed'; ...
    'MB/MBR speed',                          500,    'px/s',   'speed'; ...
    'SDG spatial frequency',                 0.0042, 'cyc/px', 'sf'; ...
    'SDG spatial frequency',                 0.0083, 'cyc/px', 'sf'; ...
    'SDG spatial frequency (recordings 1 and 8)', 0.0124, 'cyc/px', 'sf'};

n = size(spec, 1);                                                  % number of quantities
Quantity = string(spec(:, 1));                                      % quantity names
PixelValue = cell2mat(spec(:, 2));                                  % values in pixel units
PixelUnit = string(spec(:, 3));                                     % pixel units
Centre = nan(n, 1);  Edge = nan(n, 1);                              % degree values at the centre and at the edge
Unit = strings(n, 1);  Formula = strings(n, 1);                     % degree unit and formula used
for k = 1:n
    v = PixelValue(k);                                              % this quantity in pixel units
    switch spec{k, 4}
        case 'extent'                                               % whole screen: both edges as positions
            Centre(k) = 2 * posDeg(v / 2 * px);                     % edge to edge = 2*atan(16.5/21.5)
            Unit(k) = "deg";  Formula(k) = "2*atan((W/2)/D)";
        case 'position'                                             % a point v px from the centre
            Centre(k) = posDeg(v * px);                             % eccentricity of that point
            Unit(k) = "deg";  Formula(k) = "atan(d/D)";
        case 'size'                                                 % an object v px wide
            Centre(k) = sizeC(v * px);                              % centred on the line of sight
            if v * px <= halfW                                      % fits between the centre and the edge
                Edge(k) = sizeE(v * px);                            % same object at the edge of the screen
            end
            Unit(k) = "deg";  Formula(k) = "2*atan(s/(2D))";
        case 'speed'                                                % px/s -> deg/s through the local scale
            Centre(k) = v * px * scaleC;                            % at the centre
            Edge(k) = v * px * scaleE;                              % radial motion at the horizontal edge
            Unit(k) = "deg/s";  Formula(k) = "v*p/D (rad/s)";
        case 'sf'                                                   % cyc/px -> cyc/deg through the local scale
            Centre(k) = v / (px * scaleC);                          % cyc/px divided by deg/px at the centre
            Edge(k) = v / (px * scaleE);                            % radial direction at the horizontal edge
            Unit(k) = "cyc/deg";  Formula(k) = "f*D/p (per rad)";
    end
end

T = table(Quantity, PixelValue, PixelUnit, Centre, Edge, Unit, Formula, ...
    'VariableNames', {'Quantity', 'PixelValue', 'PixelUnit', 'AtScreenCentre', 'AtScreenEdge', 'Unit', 'Formula'});   % output table
T.AtScreenCentre = round(T.AtScreenCentre, 4);                      % 4 decimals is beyond display precision
T.AtScreenEdge = round(T.AtScreenEdge, 4);                          % NaN where an edge value does not apply
disp(T);                                                            % print for the log
fprintf(['Flat screen, D = %.1f cm, pixel = %.5f cm. AtScreenCentre: centred / local value at the ' ...
    'screen centre. AtScreenEdge: object touching the horizontal screen edge (16.5 cm), or local radial ' ...
    'value there. Positions: atan(d/D); sizes: 2*atan(s/(2D)).\n'], D, px);   % how to read the table
writetable(T, outCsv, 'QuoteStrings', true);                        % quoted text so readtable parses it back
fprintf('Saved %s\n', outCsv);                                      % where it went
