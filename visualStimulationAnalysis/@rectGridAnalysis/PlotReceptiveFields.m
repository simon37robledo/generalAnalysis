function [colorbarLims] = PlotReceptiveFields(obj, params)
% PlotReceptiveFields  Plots spatial receptive fields from a rectangular-grid
%                      (flash) stimulus, split by luminance, size, and
%                      On/Off response polarity.
%
% OUTPUT
%   colorbarLims  – [cMin cMax] limits of the last colorbar drawn.

% ── Input argument block ──────────────────────────────────────────────────────
arguments (Input)
    obj
    params.overwrite           logical = false   % Overwrite existing saved figures
    params.analysisTime                = datetime('now')  % Timestamp for provenance
    params.inputParams                 = false           % If true, print params and exit
    params.exNeurons                   = 1;              % Explicit neuron index/indices to plot
    params.AllSomaticNeurons           = false           % Plot every somatic unit
    params.AllResponsiveNeurons        = false           % Plot only statistically responsive units
    params.noEyeMoves                  = false           % Use no-eye-movement recording mode
    params.reduceFactor                = 20             % Spatial downsampling factor (fallback)
    params.allStimParamsCombined       = false           % Plot the fully-combined (summed) RF
    params.RFsDivision                 = {'On-Off','',''} % Which dims to show separately vs. average: {response, lum, size}
    params.eye_to_monitor_distance     = 21.5            % Eye-to-monitor distance in cm
    params.pixel_size                  = 33              % Physical pixel size in µm
    params.resolution                  = 1080            % Monitor vertical resolution in pixels
    params.meanAllNeurons              = false           % Average RF across all neurons before plotting
    params.TypeOfResponse     string   = "on"            % Which polarity to plot: "on", "off", or "both"
    params.PaperFig                    = false           % Apply paper-quality rendering settings
    params.colorbarLims                = []              % Optional manual colorbar limits [cMin cMax]
    params.tickNum                     = 3               % Number of ticks per axis
    params.smooth             logical  = false           % true: old Gaussian-smoothed pixel map (RFuFilt); false: unsmoothed 9x9 stimulus-grid map (paper, 2026-11-17)
end

% ── Debug helper: print parameter struct and exit ─────────────────────────────
if params.inputParams, disp(params), return, end

% ── Load pre-computed statistics and receptive fields ────────────────────────
Stats = obj.ShufflingAnalysis;          % Shuffling-based significance statistics
RFs   = obj.CalculateReceptiveFields;   % Receptive-field maps (no speed argument for this stimulus)

% Extract per-unit response p-values (no speed sub-struct for this stimulus)
pvals = Stats.pvalsResponse;

% Load the response window data to recover the stimulus condition table
responses = obj.ResponseWindow;

% Extract unique stimulus sizes and luminance levels from the condition matrix
uSize = unique(responses.C(:, 3));   % Unique stimulus sizes
uLum  = unique(responses.C(:, 4));   % Unique luminance levels

% ── Select which neurons to process ──────────────────────────────────────────
if params.AllSomaticNeurons
    % Use every unit regardless of responsiveness
    eNeuron = 1:numel(pvals);
    pvals   = [eNeuron; pvals(eNeuron)];   % Row 1: indices, Row 2: p-values
elseif params.AllResponsiveNeurons
    % Keep only units whose response p-value is below α = 0.05
    eNeuron = find(pvals < 0.05);
    pvals   = [eNeuron; pvals(eNeuron)];
    if isempty(eNeuron)
        fprintf('No responsive neurons.\n')
        return
    end
else
    % Use the explicitly provided unit index/indices
    eNeuron = params.exNeurons;
    pvals   = [eNeuron; pvals(eNeuron)];
end

% ── Build the reduced-resolution coordinate grid ──────────────────────────────
coorRect = obj.VST.rect';   % Screen rectangle [x y w h] in pixels, transposed to column vector

% Use the reduceFactor stored inside the RF struct (more reliable than params)
redCoorX = round(coorRect(3) / RFs.params.reduceFactor);   % Grid width in reduced pixels
redCoorY = round(coorRect(4) / RFs.params.reduceFactor);   % Grid height in reduced pixels

% Convert pixel size to cm at the reduced resolution
pixel_size_cm = params.pixel_size / (params.resolution / RFs.params.reduceFactor);

% Build the visual-angle (degrees) coordinate arrays for the reduced grid
monitor_resolution = [redCoorX, redCoorY];                        % [width height] in reduced pixels
% Positions need eccentricity atan(d/D) (pixels2eyePosition), not the subtended size
% 2*atan(d/2D) of pixels2eyeDegrees, which overestimates the edge (42.0 instead of 37.5 deg)
[theta_x, theta_y] = pixels2eyePosition(params.eye_to_monitor_distance, ...
                                        pixel_size_cm, monitor_resolution);

% Crop theta_x horizontally to be square (removes equal margins left and right)
theta_x = theta_x(:, 1 + (redCoorX - redCoorY)/2 : (redCoorX - redCoorY)/2 + redCoorY);

% ── Pre-compute symmetric, degree-based tick positions (shared across tiles) ─
% Computing ticks once here and reusing inside the loop avoids
% repeated interp1 calls and guarantees all subplots share identical ticks.

% X-axis: find the largest absolute visual angle present in the cropped grid
maxDeg_x  = max(abs(theta_x(1, :)));
% Define params.tickNum tick positions in degrees, symmetric about 0, stopping 5° inside the edge
tickDeg_x = linspace(-(maxDeg_x - 5), maxDeg_x - 5, params.tickNum);
% Round to whole degrees first, so each tick sits exactly at the angle its label shows
tickLbl_x = round(tickDeg_x);
% Map the rounded degree values back to reduced-pixel column indices via linear interpolation
tickPix_x = interp1(theta_x(1, :), 1:size(theta_x, 2), tickLbl_x, 'linear', 'extrap');

% Y-axis: same procedure on the first column of theta_y
maxDeg_y  = max(abs(theta_y(:, 1)));
tickDeg_y = linspace(-(maxDeg_y - 5), maxDeg_y - 5, params.tickNum);
tickLbl_y = round(tickDeg_y);   % Rounded y-axis degree labels (tick positions use these values)
tickPix_y = interp1(theta_y(:, 1), 1:size(theta_y, 1), tickLbl_y, 'linear', 'extrap');
tickLbl_y = 0 - tickLbl_y;      % Screen y grows downwards: flip the label sign so up is positive (0 - x avoids a "-0" label)

% Warn if requested tick range exceeds the actual RF map extent
if maxDeg_x < 5 || maxDeg_y < 5
    warning('PlotReceptiveFields: tick margin of 5° exceeds RF map extent (%.1f°, %.1f°)', ...
             maxDeg_x, maxDeg_y);
end

% ── Load the appropriate RF arrays ────────────────────────────────────────────
if params.noEyeMoves
    % No-eye-movement mode: loading not yet implemented (placeholder)
    % BUG: this branch leaves RFu and RFuFilt undefined; implement or error out.
    error('PlotReceptiveFields: noEyeMoves mode is not yet implemented for this stimulus.');
else
    % CalculateReceptiveFields averages (ball mask x rate) over ALL nCond conditions
    % (positions x sizes x luminances; mean over its dim 6), so its maps are rate / nCond,
    % not spk/s. The balls do not overlap (tiling ratio 1), so multiplying by nCond gives
    % the rate (spk/s) of the position covering each pixel (checked 2026-11-12: rec 74,
    % unit 246, peak 0.04367 x 81 = 3.537 spk/s = its best position's on-rate).
    % Display only: the cached struct (used by SpatialTuningIndex) is not changed.
    nSizeRF    = numel(unique(obj.VST.tilingRatios));                                 % sizes, counted as in CalculateReceptiveFields
    nLumsRF    = numel(unique(obj.VST.rectLuminosity(obj.VST.luminosities)));         % luminances, counted as in CalculateReceptiveFields
    trialDivRF = numel(obj.VST.pos) / numel(unique(obj.VST.pos)) / nSizeRF / nLumsRF; % repetitions per condition, as there
    nCondRF    = size(responses.C, 1) / trialDivRF;                                   % conditions the maps were averaged over
    assert(trialDivRF == round(trialDivRF) && nCondRF == round(nCondRF), ...
        'PlotReceptiveFields: %d trials do not split into whole conditions (%.3f repetitions each).', ...
        size(responses.C, 1), trialDivRF);                                            % unbalanced design: the scale would be wrong
    assert(size(RFs.RFu, 2) == nLumsRF && size(RFs.RFu, 3) == nSizeRF, ...
        'PlotReceptiveFields: RF maps have %d luminances x %d sizes, stimulus has %d x %d.', ...
        size(RFs.RFu, 2), size(RFs.RFu, 3), nLumsRF, nSizeRF);                        % cached maps must describe this stimulus
    if any(unique(obj.VST.tilingRatios) > 1)                                          % balls larger than their tile overlap
        warning('PlotReceptiveFields:overlappingBalls', ...
            'Tiling ratio > 1: where balls overlap the rescaled map sums rates, so it is not exactly spk/s there.');  % flag inexact units
    end
    RFu     = RFs.RFu     * nCondRF;   % Fully-combined RF map in spk/s [on/off × lum × size × y × x × unit]
    RFuFilt = RFs.RFuFilt * nCondRF;   % Smoothed RF map in spk/s (smoothing is linear, so the scale carries over)
end

% ── Map shown in the tiles: smoothed pixel map, or unsmoothed 9x9 grid map ────
tilePix_x = tickPix_x;  tileLbl_x = tickLbl_x;   % smoothed map: reduced-pixel ticks from above
tilePix_y = tickPix_y;  tileLbl_y = tickLbl_y;   % same for y
if params.smooth
    RFuShow = RFuFilt;                           % old display: Gaussian-smoothed pixel map
else
    % Grid map: each cell shows the rate of the position presented there = the mean of RFu*nCondRF
    % over that position's disk (exact spk/s). Block means of RFu*nCondRF are not: they scale each
    % cell by the fraction the disk covers (23-28 of 36 px, ~0.72), and the rounded disks spill
    % 1-3 px into neighbouring cells (up to 7% of the peak mixed in; checked 2026-11-18, recs 74, 81).
    nGrid = 9;                                                                         % 9x9 stimulus grid, as SpatialTuningIndex
    [masks, cellRC, sizeSlot] = rgPositionMasks(obj, responses.C, RFs.params.reduceFactor, uSize, size(RFu, 4), nGrid);  % one disk per (position, size)
    nN      = size(RFu, 6);                                                            % units in the cache
    RFflat  = reshape(RFu, [size(RFu, 1) * size(RFu, 2), size(RFu, 3), size(RFu, 4) * size(RFu, 5), nN]);  % [on/off*lum × size × pixel × unit]
    RFuShow = nan(size(RFu, 1) * size(RFu, 2), size(RFu, 3), nGrid, nGrid, nN);       % [on/off*lum × size × 9 × 9 × unit]
    filled  = false(size(RFu, 3), nGrid, nGrid);                                       % cells already given a position, per size
    for q = 1:size(masks, 3)                                                           % each (position, size)
        si = sizeSlot(q);  r = cellRC(q, 1);  c = cellRC(q, 2);                        % size slot and grid cell of this disk
        assert(~filled(si, r, c), ...
            'PlotReceptiveFields: two positions fall in grid cell (%d, %d) for size %d.', r, c, si);   % one position per cell
        filled(si, r, c) = true;                                                       % mark the cell as used
        mk = masks(:, :, q);                                                           % this disk's mask
        RFuShow(:, si, r, c, :) = mean(RFflat(:, si, mk(:), :), 3);                    % mean rate over the disk pixels (spk/s)
    end
    assert(all(filled(:)), 'PlotReceptiveFields: %d of %d grid cells have no position.', nnz(~filled), numel(filled));  % full 9x9 grid per size
    RFuShow = reshape(RFuShow, [size(RFu, 1), size(RFu, 2), size(RFu, 3), nGrid, nGrid, nN]);  % [on/off × lum × size × 9 × 9 × unit]
    % Ticks at cell centres, labelled with the angle of the centre
    bsT    = size(theta_x, 2) / nGrid;                                    % reduced pixels per cell (6)
    cenPix = ((1:nGrid) - 0.5) * bsT + 0.5;                               % cell centres in reduced-pixel coordinates
    degX   = interp1(1:size(theta_x, 2), theta_x(1, :), cenPix);          % horizontal angle at each cell centre (deg)
    degY   = interp1(1:size(theta_y, 1), theta_y(:, 1)', cenPix);         % vertical angle at each cell centre (deg)
    tilePix_x = round(linspace(1, nGrid, params.tickNum));                % tick cells: first, middle, last (for 3 ticks)
    tilePix_y = tilePix_x;                                                % same cells on y
    tileLbl_x = round(degX(tilePix_x));                                   % label = rounded angle of that cell's centre
    tileLbl_y = 0 - round(degY(tilePix_y));                               % screen y grows downwards: flip so up is positive (0 - x avoids "-0")
end

% Compute the number of spatial offset positions (used to scale the Gaussian kernel)
% SUGGESTION: sqrt(max(obj.VST.pos)) is a fragile way to recover the grid side length.
% If VST.pos contains XY positions, consider numel(unique(obj.VST.pos(:,1))) instead.
offsetN = sqrt(max(obj.VST.pos));

% Build a 2-D Gaussian smoothing kernel scaled to the RF map size
% SUGGESTION: fspecial requires the Image Processing Toolbox; imgaussfilt is more portable.
TwoDGaussian = fspecial('gaussian', ...
                         floor(size(RFu, 4) / (offsetN / 2)), ...   % Kernel window (px)
                         redCoorY / offsetN);                        % Kernel sigma (px)

% Identify which RFsDivision dimensions are EMPTY (to be averaged/combined)
% Empty cell = "combine this dimension"; non-empty = "show this dimension separately"
hasNotString = find(~cellfun(@isempty, params.RFsDivision) == 0);

% Identify which RFsDivision dimensions are NON-EMPTY (to be shown separately)
% BUG: hasString is computed but never used downstream – likely dead code or a missing feature.
hasString = find(~cellfun(@isempty, params.RFsDivision) == 1);

% ── Build response-type labels consistent with TypeOfResponse filter ──────────
% BUG FIX (original): TypeOfResponse filter was applied twice (once before
% tilesSize computation, once after creating the figure). The second pass
% would crash when TypeOfResponse is "off" because dim 1 was already size 1.
% Fixed by applying the filter only once, just before the tile loop.
% Additionally, the title used rspT{r} which always returned 'on' for "off"
% mode. Fixed by building rspLabels from the actual selected polarity.
switch params.TypeOfResponse
    case "on"
        rspLabels = {'on'};      % Single label for On-only display
    case "off"
        rspLabels = {'off'};     % Single label for Off-only display
    case "both"
        rspLabels = {'on','off'}; % Two labels when showing both polarities
    otherwise
        error('params.TypeOfResponse is not valid; options are "on", "off", "both".')
end

% ── Mean-across-neurons preprocessing ────────────────────────────────────────
if params.meanAllNeurons
    % Average the displayed map (grid or smoothed) across all neurons (dimension 6)
    RFuRed = reshape(mean(RFuShow, 6), ...
                     [size(RFuShow,1), size(RFuShow,2), ...
                      size(RFuShow,3), size(RFuShow,4), size(RFuShow,5)]);

    % Additionally average over dimensions marked as "combine" in RFsDivision
    for i = 1:numel(hasNotString)
        RFuRed = mean(RFuRed, hasNotString(i));   % Collapse dimension i by averaging
    end

    % Also collapse the summed RF used in allStimParamsCombined
    RFu = mean(sum(RFu, [1,2,3]), 6);   % Sum over condition dims, then average over neurons

    eNeuron = 1;   % Treat mean as a single "virtual" neuron
end

% ═══════════════════════════════════════════════════════════════════════════════
% Main loop: one iteration per selected neuron
% ═══════════════════════════════════════════════════════════════════════════════
for u = eNeuron

    % ── Optional: plot the fully-combined (summed across all conditions) RF ──
    if params.allStimParamsCombined

        % BUG FIX: RFu was being summed inside the loop, permanently modifying
        % it on every iteration. Use a local variable to avoid corrupting RFu
        % for subsequent neurons.
        RFu_combined = sum(RFu, [1,2,3]);   % Sum over response-type, lum, and size dimensions

        figRF = figure;   % Open a new figure window

        if params.meanAllNeurons
            % Display the pre-averaged, Gaussian-smoothed RF
            imagesc(squeeze(conv2(squeeze(RFu_combined), TwoDGaussian, 'same')));
        else
            % Display the Gaussian-smoothed combined RF for neuron u
            imagesc(squeeze(conv2(squeeze(RFu_combined(:,:,:,:,:,u)), TwoDGaussian, 'same')));
        end

        c = colorbar;              % Attach a colourbar
        title(c, 'spk/s')         % Label colourbar units
        colormap('turbo')          % Perceptually-uniform colour map
        title(sprintf('u-%d', u)) % Title with unit number

        % Apply symmetric degree-based ticks
        xticks(tickPix_x);      xticklabels(tickLbl_x);
        yticks(tickPix_y);      yticklabels(tickLbl_y);

        axis image   % Equal aspect ratio, axis box fitted tightly to image (no whitespace)

        figRF.Position = [680 577 156 139];   % Set figure size in pixels

        % Save figure to disk if overwrite flag is set
        if params.noEyeMoves
            if params.overwrite
                obj.printFig(figRF, sprintf('%s-NEM-rectGrid-ReceptiveField-eNeuron-%d.pdf', ...
                              obj.dataObj.recordingName, u));
            end
        else
            if params.overwrite
                obj.printFig(figRF, sprintf('%s-rectGrid-ReceptiveField-eNeuron-%d', ...
                              obj.dataObj.recordingName, u));
            end
        end

    end % allStimParamsCombined

    % ── Extract this neuron's RF slice and average over "combine" dimensions ─
    if ~params.meanAllNeurons
        % Extract neuron u's RF; preserve all 5 condition dimensions explicitly
        RFuRed = reshape(RFuShow(:,:,:,:,:,u), ...
                         [size(RFuShow,1), size(RFuShow,2), ...
                          size(RFuShow,3), size(RFuShow,4), size(RFuShow,5)]);   % displayed map (grid or smoothed)
        % [on/off × lum × size × y × x]

        % Average over dimensions that are marked "combine" in RFsDivision
        for i = 1:numel(hasNotString)
            RFuRed = mean(RFuRed, hasNotString(i));   % Collapse dimension i by averaging
        end
    end

    % ── Determine colour-axis limits BEFORE applying the polarity filter ──────
    % Computing cMax/cMin here (on the full On+Off data) ensures the colorbar
    % range is symmetric across polarities when TypeOfResponse is "both".
    cMax = max(RFuRed, [], 'all');   % Global maximum firing rate
    cMin = min(RFuRed, [], 'all');   % Global minimum firing rate

    % ── Apply TypeOfResponse polarity filter ─────────────────────────────────
    % BUG FIX: filter applied only once here (was applied twice in original,
    % crashing on the second pass).
    switch params.TypeOfResponse
        case "on"
            RFuRed = RFuRed(1, :, :, :, :);   % Select On-response slice (dim 1 = 1)
        case "off"
            RFuRed = RFuRed(2, :, :, :, :);   % Select Off-response slice (dim 1 = 2)
        % "both": keep RFuRed unchanged; the for-r loop handles both slices
    end

    % ── Compute tile-grid layout ──────────────────────────────────────────────
    % BUG FIX: was prod(size(RFuRed,[1 2 3])) which always produces a scalar,
    % making the numel==3 branch unreachable. Corrected to size().
    tilesSize = size(RFuRed, [1 2 3]);            % [nResponseType, nLum, nSize]
    tilesSize = [tilesSize(1), tilesSize(2) * tilesSize(3)];  % rows = polarity, cols = lum × size

    % ── Create the tiled figure ───────────────────────────────────────────────
    figRF = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full-screen window
    NeuronLayout = tiledlayout(tilesSize(1), tilesSize(2), ...
                               "TileSpacing", "tight", "Padding", "tight");

    j = 0;   % Running tile counter (used to attach colorbar to the last tile only)

    % ── Inner loops: one tile per (response-type × luminance × size) ─────────
    for r = 1:size(RFuRed, 1)     % Iterate over response-type slices (1 or 2)
        for l = 1:size(RFuRed, 2) % Iterate over luminance slices
            for s = 1:size(RFuRed, 3) % Iterate over size slices

                ax = nexttile;   % Advance to the next tile in the layout

                % Display the RF heat map for condition (r, l, s)
                imagesc(squeeze(RFuRed(r, l, s, :, :)));

                % BUG FIX: original code assigned 'xi = gca' (unused) then
                % immediately used 'axi' which was not yet defined, causing
                % "Undefined variable 'axi'" error. Fixed: use axi throughout.
                axi = gca;

                % Style axis tick labels
                axi.YAxis.FontSize = 8;
                axi.YAxis.FontName = 'helvetica';
                axi.XAxis.FontSize = 8;
                axi.XAxis.FontName = 'helvetica';

                % Axis labels (degrees of visual angle)
                xlabel('Degrees', 'FontSize', 10, 'FontName', 'helvetica')
                ylabel('Degrees', 'FontSize', 10, 'FontName', 'helvetica')

                % Lock colour axis to the per-neuron global range for comparability
                % SUGGESTION: clim() is the modern replacement for deprecated caxis()
                clim([cMin, cMax]);

                colormap('turbo')   % Perceptually-uniform colour map

                % BUG FIX: original used rspT{r} which always returned 'on'
                % when TypeOfResponse=="off" (because dim 1 is already size 1
                % after filtering). Fixed: use rspLabels built from TypeOfResponse.
                title(sprintf('respType-%s-Lum-%s-Size-%s', ...
                              rspLabels{r}, string(uLum(l)), string(uSize(s))), ...
                      'FontSize', 4)

                % Apply symmetric degree-based ticks
                xticks(tilePix_x);      xticklabels(tileLbl_x);   % x ticks in degrees (grid cells, or reduced pixels if smooth)
                yticks(tilePix_y);      yticklabels(tileLbl_y);   % y ticks in degrees (grid cells, or reduced pixels if smooth)

                j = j + 1;   % Increment tile counter

                % Attach colourbar only to the final tile
                if j == size(RFuRed,1) * size(RFuRed,2) * size(RFuRed,3)
                    c = colorbar;
                    title(c, 'spk/s', 'FontSize', 8, 'FontName', 'helvetica')

                    % Override colour limits if the caller supplied explicit bounds
                    if ~isempty(params.colorbarLims)
                        clim(params.colorbarLims);
                    end

                    colorbarLims = c.Limits;   % Return the final colourbar limits
                end

                % Equal aspect ratio, axis box fitted tightly to image (no whitespace)
                % SUGGESTION: axis image supersedes both axis equal and pbaspect([1 1 1])
                axis(ax, 'image');

            end % size
        end % luminance
    end % response type

    % ── Resize figure and build filename strings ──────────────────────────────
    set(figRF, 'Units', 'centimeters');
    set(figRF, 'Position', [2 2 4 4]);   % Compact size for single-tile or small layouts

    Slum = strjoin(string(uLum), "-");   % Luminance values joined for filename, e.g. "1-255"

    % ── Handle the mean-all-neurons special case ──────────────────────────────
    if params.meanAllNeurons
        title(NeuronLayout, 'MeanAllUnits');   % Label the whole layout
        if params.overwrite
            % BUG: fieldName is never defined in this function.
            % Removed fieldName from the format string below.
            obj.printFig(figRF, sprintf('%s-RectGrid-RF-sep-%s-Mean', ...
                          obj.dataObj.recordingName, strjoin(params.RFsDivision, '&')));
        end
        return   % Only one figure produced; exit after the mean neuron
    end

    % ── Save figure ───────────────────────────────────────────────────────────
    if ~params.noEyeMoves
        if params.overwrite
            if params.PaperFig
                obj.printFig(figRF, sprintf('%s-RectGrid-RF-lum-%s-eNeuron-%d', ...
                              obj.dataObj.recordingName, Slum, u), "PaperFig", true);
            else
                % BUG FIX: fieldName was used here but is never defined in this function.
                % Replaced with params.TypeOfResponse for a meaningful filename component.
                obj.printFig(figRF, sprintf('%s-RectGrid-RF-sep-%s-%s-eNeuron-%d', ...
                              obj.dataObj.recordingName, params.TypeOfResponse, ...
                              strjoin(params.RFsDivision, '&'), u));
            end
        end
    end

    % Close the figure unless this is the last neuron (keep last open for inspection)
    if u ~= eNeuron(end)
        close
    end

end % neuron loop

end % function


function [masks, cellRC, sizeSlot] = rgPositionMasks(obj, C, reduceFactor, uSize, sR, nGrid)
% rgPositionMasks  Disk mask of every (position, size) shown, rebuilt exactly as in
% CalculateReceptiveFields (same centre, radius, rounding and transpose), so the masks
% have the orientation of the cached RFu maps.
%   masks    [sR × sR × nPosSize] logical disk masks
%   cellRC   [nPosSize × 2] grid cell (row, column) holding each disk's centre of mass
%   sizeSlot [nPosSize × 1] size slot of the RF maps (index into uSize)
[x, y] = meshgrid(1:sR, 1:sR);                                            % pixel grid of the reduced square screen
assert(round(obj.VST.rect(4) / reduceFactor) == sR, ...
    'rgPositionMasks: map side %d does not match the screen height / %g.', sR, reduceFactor);  % same grid as the cache
rd       = obj.VST.rectData;                                              % square corners per size and position (px)
pc       = unique(C(:, 2:3), 'rows');                                     % every (position, size) shown
masks    = false(sR, sR, size(pc, 1));                                    % disk masks
cellRC   = zeros(size(pc, 1), 2);                                         % grid cell of each disk
sizeSlot = zeros(size(pc, 1), 1);                                         % size slot of each disk
bs       = sR / nGrid;                                                    % reduced pixels per grid cell
for k = 1:size(pc, 1)                                                     % each position x size
    p  = pc(k, 1);                                                        % position index
    s  = pc(k, 2);                                                        % size (index into rectData, as in CalculateReceptiveFields)
    Xc = (round((rd.X2{1, s}(p) - rd.X1{1, s}(p)) / 2) + rd.X1{1, s}(p)) / reduceFactor;   % disk centre x (reduced px)
    Yc = (round((rd.Y4{1, s}(p) - rd.Y1{1, s}(p)) / 2) + rd.Y1{1, s}(p)) / reduceFactor;   % disk centre y (reduced px)
    r  = round((rd.X2{1, s}(p) - rd.X1{1, s}(p)) / 2) / reduceFactor;    % disk radius (reduced px)
    m  = (sqrt((x - Xc).^2 + (y - Yc).^2) <= r)';                         % disk mask, transposed as the cached VideoScreen
    assert(any(m(:)), 'rgPositionMasks: position %d (size %d) covers no pixel.', p, s);   % every position must be on the map
    [ii, jj] = find(m);                                                   % map row/column of the disk pixels
    masks(:, :, k) = m;                                                   % store the mask
    cellRC(k, :)   = [ceil(mean(ii) / bs), ceil(mean(jj) / bs)];          % cell containing the disk's centre of mass
    sizeSlot(k)    = find(uSize == s);                                    % size slot of the RF maps
end
end