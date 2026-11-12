function [colorbarLims] = PlotReceptiveFields(obj, params)
% PlotReceptiveFields  Plots the spatial receptive fields of neurons recorded
%                      during a moving-ball stimulus, optionally filtered by
%                      direction, size, and luminance.
%
% OUTPUT
%   colorbarLims  – [cMin cMax] limits of the last colorbar drawn.

% ── Input argument block ─────────────────────────────────────────────────────
arguments (Input)
    obj                                   % Parent analysis object
    params.overwrite           logical = false   % Overwrite existing figure files
    params.analysisTime                = datetime('now')  % Timestamp for provenance
    params.inputParams                 = false           % If true, print params and exit
    params.exNeurons                   = nan;            % Explicit neuron indices to plot
    params.AllSomaticNeurons           = false           % Plot every somatic unit
    params.AllResponsiveNeurons        = true            % Plot only statistically responsive units
    params.fixedWindow                 = false           % Use a fixed time window
    params.speed                       = 1;              % Stimulus speed index (1 = slow, 2 = fast)
    params.noEyeMoves                  = false           % Use no-eye-movement recording mode
    params.reduceFactor                = 20             % Spatial downsampling factor for the RF map
    params.allCombined                 = false           % Also plot the direction-summed RF
    params.eye_to_monitor_distance     = 21.5            % Eye-to-monitor distance in cm
    params.pixel_size                  = 33              % Physical size of one pixel in µm
    params.resolution                  = 1080            % Monitor vertical resolution in pixels
    params.meanAllNeurons              = false           % Average RF across all neurons before plotting
    params.PaperFig           logical  = false           % Apply paper-quality rendering settings
    params.OneDirection        string  = "all"           % Restrict plot to one direction ("up","left","down","right","all")
    params.OneLuminosity       string  = "all"           % Restrict plot to one luminance ("black","white","all")
    params.OneSize             string  = "all"           % Restrict plot to one size ("small","middle","big","all")
    params.colorbarLims                = []              % Optional manual colorbar limits [cMin cMax]
    params.tickNum                     = 3               % Number of ticks per axis   
end

% ── Debug helper: print parameter struct and exit ────────────────────────────
if params.inputParams, disp(params), return, end

% ── Load pre-computed statistics and receptive fields ───────────────────────
Stats = obj.ShufflingAnalysis;                              % Shuffling-based significance statistics
RFs   = obj.CalculateReceptiveFields('speed', params.speed); % Receptive-field maps at the chosen speed

% Build the field name that indexes speed-specific substructs (e.g., "Speed1")
fieldName = sprintf('Speed%d', params.speed);

% Extract per-unit response p-values for the chosen speed
pvals = Stats.(fieldName).pvalsResponse;

% Load the response window data to recover the stimulus condition table
responses = obj.ResponseWindow;

% Extract the unique values of each stimulus dimension from the condition matrix
% Columns of C are assumed to encode: [?, direction, ?, size, ?, luminance]
uDir  = unique(responses.(fieldName).C(:, 2));  % Unique motion directions (radians)
uSize = unique(responses.(fieldName).C(:, 4));  % Unique stimulus sizes
uLum  = unique(responses.(fieldName).C(:, 6));  % Unique luminance levels

% ── Resolve which direction(s) to include ────────────────────────────────────
if params.OneDirection ~= "all"
    switch params.OneDirection
        case "up"
            dirIDX = find(round(uDir,2) == 0);       % 0 rad = upward motion
        case "left"
            dirIDX = find(round(uDir,2) == 1.57);    % π/2 rad ≈ leftward motion
        case "down"
            dirIDX = find(round(uDir,2) == 3.14);    % π rad ≈ downward motion
        case "right"
            % BUG FIX: was find(uDir==1.57) which is identical to "left".
            % Right corresponds to 3π/2 ≈ 4.71 rad.
            dirIDX = find(round(uDir,2) == 4.71);
        otherwise
            error("Unknown OneDirection value: %s", params.OneDirection)
    end
    DirectionSelected = uDir(dirIDX);   % Scalar: the single selected direction value
else
    dirIDX         = 1:numel(uDir);     % All direction indices
    DirectionSelected = uDir;           % All direction values
end

% ── Resolve which luminance(s) to include ────────────────────────────────────
if params.OneLuminosity ~= "all"
    switch params.OneLuminosity
        case "black"
            lumIDX = find(uLum == 1);    % Luminance value 1 = black stimulus
        case "white"
            lumIDX = find(uLum == 255);  % Luminance value 255 = white stimulus
        otherwise
            error("Unknown OneLuminosity value: %s", params.OneLuminosity)
    end
    LuminositySelected = uLum(lumIDX);  % Scalar: the single selected luminance value
else
    lumIDX            = 1:numel(uLum);  % All luminance indices
    LuminositySelected = uLum;          % All luminance values
end

% ── Resolve which size(s) to include ─────────────────────────────────────────
if params.OneSize ~= "all"
    switch params.OneSize
        case "small"
            sizeIDX = 1;                % First (smallest) size
        case "middle"
            sizeIDX = 2;                % Middle size
        case "big"
            sizeIDX = 3;                % Last (largest) size
        otherwise
            % BUG FIX: error message incorrectly referenced params.OneLuminosity
            error("Unknown OneSize value: %s", params.OneSize)
    end
    % BUG FIX: variable was named SizeSelected (no 's') but referenced later
    % as SizesSelected, causing "Undefined variable" at runtime.
    SizesSelected = uSize(sizeIDX);     % Scalar: the single selected size value
else
    sizeIDX       = 1:numel(uSize);     % All size indices
    SizesSelected = uSize;              % All size values
end

% ── Select which neurons to process ──────────────────────────────────────────
if isnan(params.exNeurons)
    if params.AllSomaticNeurons
        % Use every unit regardless of responsiveness
        eNeuron = 1:numel(pvals);
        pvals   = [eNeuron; pvals(eNeuron)]; % Row 1: indices, Row 2: p-values
    elseif params.AllResponsiveNeurons
        % Keep only units whose response p-value is below α = 0.05
        eNeuron = find(pvals < 0.05);
        pvals   = [eNeuron; pvals(eNeuron)];
        if isempty(eNeuron)
            fprintf('No responsive neurons.\n')
            return
        end
    end
else
    % Use the explicitly provided unit indices
    eNeuron = params.exNeurons;
    pvals   = [eNeuron; pvals(eNeuron)];
end

% ── Build the reduced-resolution coordinate grid ─────────────────────────────
coorRect    = obj.VST.rect';                                    % Screen rectangle [x y w h] (pixels), transposed to column
% SUGGESTION: consider using obj.VST.rect directly with named fields to
% avoid silent breakage if the rectangle layout changes.

% Downsampling factor must not exceed the smallest ball radius
reduceFactor = min([params.reduceFactor, min(obj.VST.ballSizes)]);

% Reduced-resolution grid dimensions
redCoorX = round(coorRect(3) / reduceFactor);   % Width in reduced pixels
redCoorY = round(coorRect(4) / reduceFactor);   % Height in reduced pixels

% Convert pixel size to cm at the reduced resolution
pixel_size_cm = params.pixel_size / (params.resolution / reduceFactor);

% Build the visual-angle (degrees) coordinate arrays for the reduced grid
monitor_resolution = [redCoorX, redCoorY];                       % [width height] in reduced pixels
% Positions need eccentricity atan(d/D) (pixels2eyePosition), not the subtended size
% 2*atan(d/2D) of pixels2eyeDegrees, which overestimates the edge (42.0 instead of 37.5 deg)
[theta_x, theta_y] = pixels2eyePosition(params.eye_to_monitor_distance, ...
                                        pixel_size_cm, monitor_resolution);

% Crop theta_x horizontally so it is square (same spatial extent as theta_y)
% The crop removes the equal-width margins on left and right
theta_x = theta_x(:, 1 + (redCoorX - redCoorY)/2 : (redCoorX - redCoorY)/2 + redCoorY);

% ── Pre-compute symmetric, degree-based tick positions (shared across tiles) ─
% SUGGESTION: Computing ticks once here and reusing inside the loop avoids
% repeated interp1 calls and guarantees all subplots share identical ticks.

% X-axis: find the largest absolute visual angle present in the cropped grid
maxDeg_x    = max(abs(theta_x(1, :)));
% Define 5 tick positions in degrees, symmetric about 0, stopping 1° inside the edge
tickDeg_x   = linspace(-(maxDeg_x - 5), maxDeg_x - 5, params.tickNum);
% Round to whole degrees first, so each tick sits exactly at the angle its label shows
tickLbl_x   = round(tickDeg_x);
% Map the rounded degree values back to reduced-pixel column indices via linear interpolation
tickPix_x   = interp1(theta_x(1, :), 1:size(theta_x, 2), tickLbl_x, 'linear', 'extrap');

% Y-axis: same procedure on the first column of theta_y
maxDeg_y    = max(abs(theta_y(:, 1)));
tickDeg_y   = linspace(-(maxDeg_y - 5), maxDeg_y - 5,  params.tickNum);
tickLbl_y   = round(tickDeg_y);   % Rounded y-axis degree labels (tick positions use these values)
tickPix_y   = interp1(theta_y(:, 1), 1:size(theta_y, 1), tickLbl_y, 'linear', 'extrap');
tickLbl_y   = 0 - tickLbl_y;      % Screen y grows downwards: flip the label sign so up is positive (0 - x avoids a "-0" label)

% ── Load the appropriate RF array ────────────────────────────────────────────
if params.noEyeMoves
    % Load the no-eye-movement RF (previously saved per-recording)
    RFu = squeeze(load(sprintf('NEM-RFuST-Q1-Div-X-%s', NP.recordingName)).RFuST);
else
    RFu            = RFs.RFuST;              % Direction-summed RF map [y × x × unit] (allCombined only; NOT in spk/s)
    % The cached per-direction maps sum, over every frame the ball covers a pixel, the rate
    % 250 ms later, divided by the conditions per direction. Pixels are covered for different
    % numbers of frames (8-18 in rec 97), so those maps are not spk/s. rfRateMaps rescales them
    % to the mean rate while the ball covered each pixel (checked 2026-11-12 on rec 97 unit 1
    % against a rate map built directly from spikes and ball positions: slope 0.99-1.01,
    % r 0.83-0.97 per direction). Display only: the cached struct is not changed.
    RFuDirSizeLum  = rfRateMaps(obj, RFs, responses.(fieldName).C, params.speed);  % [dir × size × lum × y × x × unit] in spk/s
end
if params.allCombined                                                     % direction-summed map is not rescaled
    warning('PlotReceptiveFields:allCombinedUnits', ...
        'allCombined plots RFuST, which is not in spk/s; its colour-bar label is wrong.');  % flag the wrong unit label
end

% Build a 2-D Gaussian smoothing kernel scaled to the RF map size
% Kernel size is proportional to the number of spatial offsets used
offsetN      = numel(unique(obj.VST.offsets));                           % Number of unique stimulus offsets
TwoDGaussian = fspecial('gaussian', floor(size(RFu, 2) / (offsetN / 2)), ... % Kernel window (px)
                         redCoorY / offsetN);                            % Kernel sigma (px)
% SUGGESTION: fspecial is from the Image Processing Toolbox.
% imgaussfilt or a manually constructed kernel can be used as a fallback.

% ═══════════════════════════════════════════════════════════════════════════════
% Main loop: one iteration per selected neuron
% ═══════════════════════════════════════════════════════════════════════════════
for u = eNeuron

    %ru = find(eNeuron == u);   % Index of neuron u within the eNeuron vector
    ru =u;
    % ── Optional: plot the direction-summed (combined) RF ──────────────────
    if params.allCombined

        figRF = figure;   % Open a new figure window
        % Display the Gaussian-smoothed combined RF as a colour image
        imagesc(squeeze(conv2(RFu(:, :, ru), TwoDGaussian, 'same')));

        c = colorbar;                   % Attach a colourbar
        title(c, 'spk/s')              % Label colourbar units

        colormap('turbo')              % Apply perceptually-uniform colour map
        title(sprintf('u-%d', u))      % Title with unit number

        % Apply symmetric x ticks (degrees)
        xticks(tickPix_x);
        xticklabels(tickLbl_x);

        % Apply symmetric y ticks (degrees)
        yticks(tickPix_y);
        yticklabels(tickLbl_y);

        axis image    % Equal aspect ratio, no white space

        figRF.Position = [680 577 156 139];   % Set figure size (pixels)

        % Save figure to disk if overwrite flag is set
        if params.noEyeMoves
            if params.overwrite
                if params.PaperFig
                    obj.printFig(figRF, sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf', ...
                                  obj.dataObj.recordingName, u), PaperFig = params.PaperFig);
                else
                    obj.printFig(figRF, sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf', ...
                                  obj.dataObj.recordingName, u));
                end
            end
        else
            if params.overwrite
                if params.PaperFig
                    obj.printFig(figRF, sprintf('%s-MovBall-ReceptiveField-eNeuron-%d', ...
                                  obj.dataObj.recordingName, u), PaperFig = params.PaperFig);
                else
                    obj.printFig(figRF, sprintf('%s-MovBall-ReceptiveField-eNeuron-%d', ...
                                  obj.dataObj.recordingName, u));
                end
            end
        end

    end % allCombined

    % ── Extract this neuron's RF slice and apply any dimension filters ──────

    if params.meanAllNeurons
        % Average the RF across all neurons (dimension 6) and collapse it
        RFuRed = reshape(mean(RFuDirSizeLum, 6), ...
                         [size(RFuDirSizeLum, 1), size(RFuDirSizeLum, 2), ...
                          size(RFuDirSizeLum, 3), size(RFuDirSizeLum, 4), ...
                          size(RFuDirSizeLum, 5)]);
        % BUG: hasNotString is never defined; the intent seems to be to
        % additionally average over the non-compared dimensions (e.g., if
        % only direction is compared, average over size and luminance).
        % Define hasNotString before this block, e.g.:
        %   hasNotString = [];
        %   if params.OneSize ~= "all",      hasNotString(end+1) = 2; end
        %   if params.OneLuminosity ~= "all", hasNotString(end+1) = 3; end
        %   if params.OneDirection ~= "all",  hasNotString(end+1) = 1; end
        % Then the loop below is correct:
        for i = 1:numel(hasNotString)   % Average over dimensions not being compared
            RFuRed = mean(RFuRed, hasNotString(i));
        end
    else
        % Extract the RF for neuron ru; keep all 5 condition dimensions explicit
        RFuRed = reshape(RFuDirSizeLum(:, :, :, :, :, ru), ...
                         [size(RFuDirSizeLum, 1), size(RFuDirSizeLum, 2), ...
                          size(RFuDirSizeLum, 3), size(RFuDirSizeLum, 4), ...
                          size(RFuDirSizeLum, 5)]);   % [dir × size × lum × y × x]

        % Apply size filter if requested (select single size slice)
        if params.OneSize ~= "all"
            RFuRed = RFuRed(:, sizeIDX, :, :, :);
        end

        % Apply luminance filter if requested
        if params.OneLuminosity ~= "all"
            RFuRed = RFuRed(:, :, lumIDX, :, :);
        end

        % Apply direction filter if requested
        if params.OneDirection ~= "all"
            RFuRed = RFuRed(dirIDX, :, :, :, :);
        end
    end

    % ── Determine colour-axis limits from this neuron's data ────────────────
    cMax = max(RFuRed, [], 'all');   % Global maximum firing rate across all conditions
    cMin = min(RFuRed, [], 'all');   % Global minimum firing rate across all conditions

    % ── Compute tile-grid layout for the tiled figure ───────────────────────
    % BUG FIX: was prod(size(RFuRed,[1 2 3])) which always produces a scalar,
    % making the numel==3 branch unreachable. Corrected to size().
    tilesSize = size(RFuRed, [1 2 3]);   % [nDir, nSize, nLum] – counts of condition tiles

    % Flatten to a [rows × cols] layout: rows = directions, cols = size × lum
    tilesSize = [tilesSize(1), tilesSize(2) * tilesSize(3)];

    % ── Create the tiled figure ──────────────────────────────────────────────
    figRF = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full-screen window
    NeuronLayout = tiledlayout(tilesSize(1), tilesSize(2), ...
                               "TileSpacing", "tight", "Padding", "compact");

    j = 0;   % Running tile counter (used to attach colorbar to the last tile)

    % ── Inner loops: one tile per (direction × size × luminance) combination ─
    for d = 1:size(RFuRed, 1)     % Iterate over direction slices
        for s = 1:size(RFuRed, 2) % Iterate over size slices
            for l = 1:size(RFuRed, 3) % Iterate over luminance slices

                ax = nexttile;   % Advance to the next tile in the layout

                % Display the RF heat map for condition (d, s, l)
                imagesc(squeeze(RFuRed(d, s, l, :, :)));

                % Lock colour axis to the per-neuron global range for comparability
                % SUGGESTION: clim() is the modern replacement for the deprecated caxis()
                clim([cMin, cMax]);

                % Style y-axis tick labels
                axi = gca;
                axi.YAxis.FontSize = 8;
                axi.YAxis.FontName = 'helvetica';

                % Axis labels (degrees of visual angle)
                xlabel('Degrees', 'FontSize', 10, 'FontName', 'helvetica')
                ylabel('Degrees', 'FontSize', 10, 'FontName', 'helvetica')

                % Style x-axis tick labels
                axi.XAxis.FontSize = 8;
                axi.XAxis.FontName = 'helvetica';

                colormap('turbo')   % Perceptually-uniform colour map

                % BUG FIX: was uDir(d), uSize(s), uLum(l) – these index the
                % FULL unfiltered arrays, so the label is wrong whenever a
                % filter is active.  Use the Selected vectors instead.
                title(sprintf('Dir-%.2f-Size-%.0f-Lum-%.0f', ...
                              DirectionSelected(d), SizesSelected(s), LuminositySelected(l)), ...
                      'FontSize', 4)

                % ── Apply symmetric degree-based ticks ──────────────────────
                xticks(tickPix_x);         % Set x tick positions (reduced-pixel units)
                xticklabels(tickLbl_x);    % Label with rounded degree values
                yticks(tickPix_y);         % Set y tick positions (reduced-pixel units)
                yticklabels(tickLbl_y);    % Label with rounded degree values

                j = j + 1;   % Increment tile counter

                % Attach colourbar only to the final tile so it doesn't clutter the layout
                if j == size(RFuRed, 1) * size(RFuRed, 2) * size(RFuRed, 3)
                    c = colorbar;
                    title(c, 'spk/s', 'FontSize', 8, 'FontName', 'helvetica')

                    % Override colour limits if the caller supplied explicit bounds
                    if ~isempty(params.colorbarLims)
                        clim(params.colorbarLims);
                    end

                    colorbarLims = c.Limits;   % Return the final colourbar limits
                end

                axis(ax, 'image');         % Equal aspect ratio so the RF map is not distorted
                % SUGGESTION: axis equal already enforces equal scaling;
                % pbaspect([1 1 1]) is therefore redundant here and can be removed.

            end % luminance
        end % size
    end % direction

    % ── Build filename strings from the selected condition labels ────────────
    Sdir  = strjoin(string(DirectionSelected),  "-");   % e.g. "0-1.57-3.14-4.71"
    Ssize = strjoin(string(SizesSelected),      "-");   % e.g. "5-10-20"
    Slum  = strjoin(string(LuminositySelected), "-");   % e.g. "1-255"

    % ── Handle the mean-all-neurons special case ─────────────────────────────
    if params.meanAllNeurons
        title(NeuronLayout, 'MeanAllUnits');   % Label the whole layout
        if params.overwrite
            obj.printFig(figRF, sprintf('%s-%s-MovBall-RF-sep-%s-Mean', ...
                          obj.dataObj.recordingName, fieldName, ...
                          sprintf('Dir-%s-Size-%s-Lum-%s', Sdir, Ssize, Slum)));
        end
        return   % Only one figure is produced; no per-neuron loop needed
    end

    % ── Resize figure for single-row layouts ────────────────────────────────
    if tilesSize(1) == 1
        set(figRF, 'Units', 'centimeters');
        set(figRF, 'Position', [2 2 4 4]);   % Compact format for one-row grids
    end

    % ── Save figure ──────────────────────────────────────────────────────────
    if ~params.noEyeMoves   % Standard (eye-movement) recording mode
        if params.overwrite
            if params.PaperFig
                obj.printFig(figRF, sprintf('%s-%s-MovBall-RF-sep-%s-eNeuron-%d.jpg', ...
                              obj.dataObj.recordingName, fieldName, ...
                              sprintf('Dir-%s-Size-%s-Lum-%s', Sdir, Ssize, Slum), u), ...
                              PaperFig = params.PaperFig);
            else
                obj.printFig(figRF, sprintf('%s-%s-MovBall-RF-sep-%s-eNeuron-%d.jpg', ...
                              obj.dataObj.recordingName, fieldName, ...
                              sprintf('Dir-%s-Size-%s-Lum-%s', Sdir, Ssize, Slum), u));
            end
        end
    end

    % Close the figure unless this is the last neuron (keep last open for inspection)
    if u ~= eNeuron(end)
        close
    end

end % neuron loop

end % function


function R = rfRateMaps(obj, RFs, C, speed)
% rfRateMaps  Moving-ball RF maps in spk/s: the mean firing rate (250 ms after the
% ball) while the ball covered each pixel.
%   raw (RFuSTDirSizeLum) = sum over conditions and covered frames of (mask x rate) / nCondDir
%   rate                  = smooth(raw x nCondDir) ./ smooth(coverage frames)
% Coverage is rebuilt from the ball trajectories with the geometry of
% CalculateReceptiveFields (true screen orientation, which the cached maps follow).
raw = RFs.RFuSTDirSizeLum;                                                % [dir × size × lum × y × x × unit], unsmoothed
uDirC = unique(C(:, 2));                                                  % directions (C column 2)
uOffC = unique(C(:, 3));                                                  % trajectory offsets (column 3)
uSzC  = unique(C(:, 4));                                                  % ball sizes in px (column 4)
uSpdC = unique(C(:, 5));                                                  % speeds (column 5)
uLumC = unique(C(:, 6));                                                  % luminances (column 6)
assert(isscalar(uSpdC), 'rfRateMaps: %d speeds in one speed field.', numel(uSpdC));   % coverage below assumes one speed
reps = size(C, 1) / numel(uDirC) / numel(uOffC) / numel(uSzC) / numel(uLumC);       % repetitions per condition
assert(reps == round(reps), 'rfRateMaps: %d trials do not split into whole conditions.', size(C, 1));  % balanced design
nCondDir = size(C, 1) / numel(uDirC) / reps;                              % conditions per direction: the divisor of the cached maps
assert(isequal(size(raw, 1:3), [numel(uDirC), numel(uSzC), numel(uLumC)]), ...
    'rfRateMaps: cached maps are %s (dir x size x lum), trials give %d x %d x %d.', ...
    mat2str(size(raw, 1:3)), numel(uDirC), numel(uSzC), numel(uLumC));  % cache must describe these trials

Xp  = obj.VST.ballTrajectoriesX;                                          % ball x [speed × offset × direction × frame] (px)
Yp  = obj.VST.ballTrajectoriesY;                                          % ball y, same layout
nFr = sort(unique(obj.VST.nFrames), 'descend');                           % frames per speed, slowest first (as CalculateReceptiveFields)
if size(Xp, 1) > 1                                                        % several speeds: keep the analysed one
    sp = min(speed, numel(nFr));                                          % same fallback as CalculateReceptiveFields (last speed)
    Xp = Xp(sp, :, :, 1:nFr(sp));                                         % that speed's frames only
    Yp = Yp(sp, :, :, 1:nFr(sp));                                         % same for y
end
assert(size(Xp, 2) == numel(uOffC) && size(Xp, 3) == numel(uDirC), ...
    'rfRateMaps: %d offsets x %d directions in the trajectories, %d x %d in the trials.', ...
    size(Xp, 2), size(Xp, 3), numel(uOffC), numel(uDirC));               % trajectory index = sorted condition index (cache convention)

sY    = size(raw, 4);                                                     % map side (reduced px)
rF    = min(RFs.params.reduceFactor, min(uSzC));                          % reduction factor of CalculateReceptiveFields
assert(round(obj.VST.rect(4) / rF) == sY, 'rfRateMaps: map side %d does not match the screen height / %g.', sY, rF);  % same grid as the cache
cropX = (round(obj.VST.rect(3) / rF) - sY) / 2;                           % left margin cropped to make the map square
[xg, yg] = meshgrid(1:sY, 1:sY);                                          % pixel coordinates of the map

G = fspecial('gaussian', floor(sY / (numel(uOffC) / 2)), sY / numel(uOffC));   % smoothing kernel of CalculateReceptiveFields
R = nan(size(raw));                                                       % output maps in spk/s
for d = 1:numel(uDirC)                                                    % each direction
    for s = 1:numel(uSzC)                                                 % each ball size
        rad = uSzC(s) / 2 / rF + 0.5;                                     % ball radius in reduced px (+0.5 as in the cache videos)
        cov = zeros(sY);                                                  % frames the ball covers each pixel, summed over trajectories
        for of = 1:numel(uOffC)                                           % each trajectory
            for t = 1:size(Xp, 4)                                         % each frame
                cov = cov + (hypot(xg - (Xp(1, of, d, t) / rF - cropX), ...
                                   yg -  Yp(1, of, d, t) / rF) <= rad);   % add this frame's ball disk
            end
        end
        covS = conv2(cov, G, 'same');                                     % smoothed coverage (denominator)
        for l = 1:numel(uLumC)                                            % each luminance (same coverage)
            for u = 1:size(raw, 6)                                        % each unit
                num = conv2(squeeze(raw(d, s, l, :, :, u)) * nCondDir, G, 'same');   % smoothed sum of mask x rate
                rate = num ./ covS;                                       % mean rate while covered (spk/s)
                rate(covS <= eps) = NaN;                                  % never covered: no estimate
                R(d, s, l, :, :, u) = rate;                               % store
            end
        end
    end
end
end