function R = estimateLfpSurface(apMetaPath, params)
% estimateLfpSurface  Brain-surface position (probe y, um from tip) from the LFP depth profile.
%
%   R = estimateLfpSurface(apMetaPath)
%   R = estimateLfpSurface(apMetaPath, plot=true, saveDir="...")
%
%   MATLAB port of the collaborator's electrophysiological_profiling.ipynb (PV140):
%     1. one probe column (shank, x = 27 um), sites sorted by y (um from the tip);
%     2. notch 50 Hz (Q = 30) + Butterworth low-pass 300 Hz (zero-phase), resample to 1 kHz;
%     3. a 2 s window (default 99-101 s from recording start);
%     4. C-factor: per-channel median removed, gain = median over time of
%        (depth-smoothed signal / signal), Gaussian across channels sigma = 3;
%     5. Welch power (scipy defaults: Hann 256, 50% overlap, nfft 4096, mean detrend),
%        mean power in 5 bands, in dB;
%     6. per band, a continuous 6-segment piecewise-linear fit of power vs depth
%        (as pwlf.fit(6)); boundary = topmost interior breakpoint;
%     7. surface = mean of the band boundaries (median and SD also returned).
%   Only the window plus padding is loaded (the notebook filters the whole recording;
%   with zero-phase filters and >= 5 s padding the 2 s window is the same).
%
%   The method finds the transition between brain and the out-of-brain part of the probe.
%   It is only valid when the surface lies inside the recorded column: a boundary within
%   edgeGuardUm of the top site, or bands disagreeing by > bandSpreadUm, sets R.nearEdge /
%   R.bandsDisagree and warns.
%
%   Output R: metaFile, shank, xCoord, yUm, fs, freqs, psdDb, bandPowerDb, fits,
%   bandBoundaryUm, surfaceMeanUm, surfaceMedianUm, surfaceStdUm, cFactor, corr,
%   nearEdge, bandsDisagree, params.

arguments
    apMetaPath (1,1) string                         % full path of the SpikeGLX *.ap.meta file
    params.shank double = []                        % probe shank (0-based); [] = shank with most sites at xCoord
    params.xCoord (1,1) double = 27                 % column x position within the shank (um)
    params.minSites (1,1) double = 96               % fewest sites accepted in the column (48-site bank recordings are not whole-shank tests)
    params.profileStart_s (1,1) double = 99         % start of the analysed window (s from recording start)
    params.profileDur_s (1,1) double = 2            % duration of the analysed window (s)
    params.pad_s (1,1) double = 5                   % padding loaded on each side for the filters (s)
    params.lineNoiseHz (1,1) double = 50            % notch frequency (Hz)
    params.notchQ (1,1) double = 30                 % notch quality factor
    params.lowpassHz (1,1) double = 300             % low-pass cut-off (Hz)
    params.lowpassOrder (1,1) double = 4            % Butterworth order (collaborator's value unverified)
    params.fsNew (1,1) double = 1000                % target sampling rate (Hz)
    params.welchSeg (1,1) double = 256              % Welch segment length (samples, scipy default)
    params.nfft (1,1) double = 4096                 % FFT length
    params.bands double = [1.5 4; 4 10; 10 30; 30 80; 80 200]   % frequency bands (Hz), one per row
    params.nSegments (1,1) double = 6               % line segments of the piecewise fit
    params.cSigma (1,1) double = 3                  % Gaussian SD across channels for the C-factor
    params.nStarts (1,1) double = 10                % independent differential-evolution runs (best kept)
    params.seed (1,1) double = 1                    % random seed (reproducible fit)
    params.edgeGuardUm (1,1) double = 100           % boundary closer than this to the top site -> warn
    params.bandSpreadUm (1,1) double = 200          % band boundaries spread wider than this -> warn
    params.plot (1,1) logical = false               % draw the QC figures
    params.saveDir (1,1) string = ""                % folder for the QC figures ("" = do not save)
end

%% Metadata and probe column
[metaDir, metaName, metaExt] = fileparts(apMetaPath);                 % folder and file name of the meta file
assert(isfile(apMetaPath), 'estimateLfpSurface:noMeta', 'Meta file not found: %s', apMetaPath);   % fail loudly
binName = replace(metaName + metaExt, ".ap.meta", ".ap.bin");          % matching binary file name
assert(isfile(fullfile(metaDir, binName)), 'estimateLfpSurface:noBin', 'Binary not found next to %s', apMetaPath);
meta = ReadMeta(char(binName), char(metaDir));                         % SpikeGLX metadata (all values are strings)
fs = str2double(meta.imSampRate);                                      % AP sampling rate (Hz)
assert(isfield(meta, 'snsGeomMap'), 'estimateLfpSurface:noGeom', ...   % geometry is needed for site positions
    'No snsGeomMap in %s (older metadata are not supported).', apMetaPath);
geom = regexp(meta.snsGeomMap, '\((\d+):(\d+):(\d+):(\d+)\)', 'tokens');   % one (shank:x:y:used) entry per saved AP channel
geom = str2double(vertcat(geom{:}));                                   % [nAP x 4] numeric: shank, x, y, used
nAP = size(geom, 1);                                                   % number of AP channels (rows 1..nAP of the binary)
onCol = geom(:, 2) == params.xCoord;                                   % sites in the requested column (x)
if isempty(params.shank)                                               % shank not given: choose automatically
    shanks = unique(geom(onCol, 1));                                   % shanks that have sites at this x
    assert(~isempty(shanks), 'estimateLfpSurface:noColumn', 'No sites at x = %g um in %s.', params.xCoord, apMetaPath);
    nPer = arrayfun(@(s) sum(onCol & geom(:, 1) == s), shanks);        % sites per shank at this x
    [~, best] = max(nPer);                                             % shank with the longest column
    shank = shanks(best);                                              % chosen shank
    if sum(nPer == nPer(best)) > 1                                     % tie between shanks: report the choice
        warning('estimateLfpSurface:shankTie', 'Several shanks have %d sites at x = %g; using shank %d.', ...
            nPer(best), params.xCoord, shank);
    end
else
    shank = params.shank;                                              % shank given by the caller
end
sel = find(onCol & geom(:, 1) == shank);                               % channel rows (1-based) of the column
if numel(sel) < params.minSites                                        % too few sites for the fit: list what exists
    [u, ~, ic] = unique(geom(:, 1:2), 'rows');                         % available (shank, x) combinations
    avail = strjoin(compose('shank %d x %d: %d sites', [u accumarray(ic, 1)]), '; ');   % readable list (char)
    error('estimateLfpSurface:fewSites', 'Shank %d, x = %g has %d sites (< %d). Available: %s', ...
        shank, params.xCoord, numel(sel), params.minSites, avail);
end
[yUm, order] = sort(geom(sel, 3));                                     % site y positions (um from tip), ascending
chIdx = sel(order);                                                    % channel rows in depth order
assert(numel(unique(yUm)) == numel(yUm), 'estimateLfpSurface:dupY', 'Duplicate site positions in the column.');

%% Load the window with padding
fileDur = str2double(meta.fileTimeSecs);                               % recording length (s)
assert(params.profileStart_s + params.profileDur_s <= fileDur, 'estimateLfpSurface:window', ...
    'Window %g-%g s exceeds the recording (%g s).', params.profileStart_s, params.profileStart_s + params.profileDur_s, fileDur);
t0 = max(0, params.profileStart_s - params.pad_s);                     % load start (s), not before the recording
t1 = min(fileDur, params.profileStart_s + params.profileDur_s + params.pad_s);   % load end (s)
if t0 > params.profileStart_s - params.pad_s || t1 < params.profileStart_s + params.profileDur_s + params.pad_s
    warning('estimateLfpSurface:shortPad', 'Padding truncated at the recording edge; filter edge effects possible.');
end
samp0 = round(t0 * fs);                                                % first sample to read
nSamp = round((t1 - t0) * fs);                                         % samples to read
raw = ReadBin(samp0, nSamp, meta, char(binName), char(metaDir));       % [nSavedChans x nSamp] int16 values as double
raw = GainCorrectIM(raw * 1e6, chIdx, meta);                           % to microvolts on the column channels (as NPAPRecording.getData)
x = raw(chIdx, :);                                                     % [nSites x nSamp] column data, depth order
clear raw                                                              % free memory (all channels)
assert(all(chIdx <= nAP), 'estimateLfpSurface:sync', 'Column includes a non-AP channel.');   % never read the sync channel

%% Preprocessing: notch, low-pass, resample (all zero-phase or linear-phase)
[bN, aN] = notchCoeffs(params.lineNoiseHz, params.notchQ, fs);         % notch coefficients identical to scipy.signal.iirnotch
x = filtfilt(bN, aN, x')';                                             % zero-phase notch along time
[bL, aL] = butter(params.lowpassOrder, params.lowpassHz / (fs / 2), 'low');   % Butterworth low-pass
x = filtfilt(bL, aL, x')';                                             % zero-phase low-pass along time
q = round(fs / params.fsNew);                                          % integer decimation factor (30000.28 / 1000 -> 30)
fsOut = fs / q;                                                        % resulting rate (1000.009 Hz; within 1e-3 of target)
assert(abs(fsOut - params.fsNew) / params.fsNew < 1e-3, 'estimateLfpSurface:fs', ...
    'fs = %g Hz is not close to an integer multiple of %g Hz.', fs, params.fsNew);
x = resample(x', 1, q)';                                               % anti-aliased resampling along time
iStart = round((params.profileStart_s - t0) * fsOut) + 1;              % first sample of the window in the resampled data
nWin = floor(params.profileDur_s * fsOut);                             % window length in samples (as int() in the notebook)
assert(iStart + nWin - 1 <= size(x, 2), 'estimateLfpSurface:window', 'Window exceeds the loaded data.');
lfp = x(:, iStart:iStart + nWin - 1);                                  % [nSites x nWin] preprocessed window
clear x                                                                % free memory (padded data)

%% C-factor correction
[lfp, cFactor] = applyCFactor(lfp, params.cSigma);                     % median-centred and gain-equalised across sites
fprintf('estimateLfpSurface: shank %d, x %g, %d sites (%g-%g um), C-factor %.3f-%.3f\n', ...
    shank, params.xCoord, numel(yUm), yUm(1), yUm(end), min(cFactor), max(cFactor));   % progress / QC

%% Power per band
[freqs, psd] = welchScipy(lfp, fsOut, params.welchSeg, params.nfft);   % [nSites x nFreq] PSD, scipy-equivalent
psdDb = 10 * log10(max(psd, realmin));                                 % PSD in dB (floor avoids log of 0)
nBands = size(params.bands, 1);                                        % number of bands
bandPowerDb = nan(numel(yUm), nBands);                                 % [nSites x nBands] band power (dB)
for b = 1:nBands
    inBand = freqs >= params.bands(b, 1) & freqs <= params.bands(b, 2);   % frequency bins in the band (inclusive)
    assert(any(inBand), 'estimateLfpSurface:band', 'No frequencies in band %g-%g Hz.', params.bands(b, :));
    bandPowerDb(:, b) = 10 * log10(max(mean(psd(:, inBand), 2), realmin));   % mean linear power, then dB (as notebook)
end

%% Piecewise-linear fits and boundaries
stream = RandStream('mt19937ar', 'Seed', params.seed);                 % private stream: does not touch the global rng
fits = struct('breaks', cell(1, nBands), 'yHat', [], 'sse', []);       % one fit per band
bandBoundaryUm = nan(1, nBands);                                       % boundary per band (um)
for b = 1:nBands
    fits(b) = piecewiseFit(yUm, bandPowerDb(:, b), params.nSegments, params.nStarts, stream);   % pwlf-like fit
    bandBoundaryUm(b) = fits(b).breaks(end - 1);                       % topmost interior breakpoint (pwlf breaks[-2])
end
surfaceMeanUm = mean(bandBoundaryUm);                                  % surface estimate (notebook: mean of bands)
surfaceMedianUm = median(bandBoundaryUm);                              % robust alternative
surfaceStdUm = std(bandBoundaryUm, 1);                                 % population SD (np.nanstd, ddof = 0)

%% Validity guards
nearEdge = yUm(end) - surfaceMeanUm < params.edgeGuardUm;              % boundary at the top of the column: surface may lie above it
bandsDisagree = max(bandBoundaryUm) - min(bandBoundaryUm) > params.bandSpreadUm;       % bands point to different transitions
if nearEdge
    warning('estimateLfpSurface:boundaryAtEdge', ...
        '%s: boundary %.0f um is within %g um of the top site (%g um); the surface is probably not in the column.', ...
        metaName, surfaceMeanUm, params.edgeGuardUm, yUm(end));
end
if bandsDisagree
    warning('estimateLfpSurface:bandsDisagree', '%s: band boundaries span %.0f um (%s).', ...
        metaName, max(bandBoundaryUm) - min(bandBoundaryUm), mat2str(round(bandBoundaryUm)));
end

%% Output
corrMat = corrcoef(lfp');                                              % site x site LFP correlation (QC figure)
corrMat(~isfinite(corrMat)) = 0;                                       % flat channels -> 0 (as np.nan_to_num)
R = struct('metaFile', apMetaPath, 'shank', shank, 'xCoord', params.xCoord, 'yUm', yUm, 'fs', fsOut, ...
    'freqs', freqs, 'psdDb', psdDb, 'bandPowerDb', bandPowerDb, 'fits', fits, ...
    'bandBoundaryUm', bandBoundaryUm, 'surfaceMeanUm', surfaceMeanUm, 'surfaceMedianUm', surfaceMedianUm, ...
    'surfaceStdUm', surfaceStdUm, 'cFactor', cFactor, 'corr', corrMat, ...
    'nearEdge', nearEdge, 'bandsDisagree', bandsDisagree, 'params', params);   % everything needed to audit the estimate
fprintf('  band boundaries (um): %s -> mean %.1f, median %.1f, SD %.1f\n', ...
    mat2str(round(bandBoundaryUm, 1)), surfaceMeanUm, surfaceMedianUm, surfaceStdUm);   % result summary

if params.plot
    plotQC(R, metaName, params.saveDir);                               % figures A-C of the notebook
end
end

% =========================================================================
function [b, a] = notchCoeffs(f0, Q, fs)
% Second-order IIR notch, same formula as scipy.signal.iirnotch(f0, Q, fs).
w0 = 2 * f0 / fs;                                                      % notch frequency, normalised to Nyquist
bw = w0 / Q;                                                           % -3 dB bandwidth, normalised
beta = tan(bw * pi / 2);                                               % scipy: sqrt(1-gb^2)/gb * tan(bw/2) with gb = 1/sqrt(2)
gain = 1 / (1 + beta);                                                 % filter gain
b = gain * [1, -2 * cos(w0 * pi), 1];                                  % numerator
a = [1, -2 * gain * cos(w0 * pi), 2 * gain - 1];                       % denominator
end

% =========================================================================
function [y, c] = applyCFactor(x, sigma)
% C-factor correction as apply_c_factor in the notebook (no common-median reference).
x = x - median(x, 2);                                                  % remove each channel's median
r = floor(4 * sigma + 0.5);                                            % kernel radius (scipy truncate = 4)
k = exp(-0.5 * ((-r:r)' / sigma) .^ 2);                                % Gaussian kernel across channels
k = k / sum(k);                                                        % unit gain
xp = [repmat(x(1, :), r, 1); x; repmat(x(end, :), r, 1)];              % replicate edges (scipy mode = "nearest")
expected = conv2(xp, k, 'valid');                                      % depth-smoothed signal, same size as x
ratio = ones(size(x));                                                 % ratio defaults to 1 where the signal is 0
nz = x ~= 0;                                                           % samples with a non-zero signal
ratio(nz) = expected(nz) ./ x(nz);                                     % expected / observed per sample
c = median(ratio, 2);                                                  % one gain per channel
c(~isfinite(c)) = 1;                                                   % undefined gain -> no correction
y = x .* c;                                                            % corrected signal
end

% =========================================================================
function [f, P] = welchScipy(x, fs, nSeg, nfft)
% Welch PSD with scipy.signal.welch defaults: Hann (periodic) window, 50% overlap,
% constant detrend per segment, one-sided density, mean across segments.
nT = size(x, 2);                                                       % samples per channel
assert(nT >= nSeg, 'estimateLfpSurface:short', 'Window (%d samples) shorter than one Welch segment.', nT);
w = 0.5 - 0.5 * cos(2 * pi * (0:nSeg - 1) / nSeg);                     % periodic Hann (scipy get_window default)
step = nSeg - floor(nSeg / 2);                                         % hop = segment - overlap
starts = 1:step:nT - nSeg + 1;                                         % segment start samples
scale = 1 / (fs * sum(w .^ 2));                                        % density scaling
nF = floor(nfft / 2) + 1;                                              % one-sided frequency bins
P = zeros(size(x, 1), nF);                                             % accumulated PSD
for s = starts
    seg = x(:, s:s + nSeg - 1);                                        % [nSites x nSeg] segment
    seg = (seg - mean(seg, 2)) .* w;                                   % remove the mean, apply the window
    X = fft(seg, nfft, 2);                                             % zero-padded FFT along time
    P = P + abs(X(:, 1:nF)) .^ 2 * scale;                              % two-sided power of the kept bins
end
P = P / numel(starts);                                                 % mean across segments
if mod(nfft, 2) == 0
    P(:, 2:end - 1) = 2 * P(:, 2:end - 1);                             % one-sided: double all but DC and Nyquist
else
    P(:, 2:end) = 2 * P(:, 2:end);                                     % odd nfft: no Nyquist bin
end
f = (0:nF - 1) * fs / nfft;                                            % frequency of each bin (Hz)
end

% =========================================================================
function F = piecewiseFit(x, y, nSeg, nStarts, stream)
% Continuous piecewise-linear least-squares fit with nSeg segments (as pwlf.fit(nSeg)).
% Interior breakpoints are searched as pwlf does: differential evolution (scipy defaults:
% best1bin, population 15 x nBreaks, dithered mutation 0.5-1, crossover 0.7, tol 0.01),
% then a local polish (fminsearch). This is repeated nStarts times and the lowest error
% is kept. For fixed breakpoints the line coefficients are an ordinary least-squares solution.
x = x(:); y = y(:);                                                    % column vectors
lo = x(1); hi = x(end);                                                % breakpoints are bounded by the data range
nb = nSeg - 1;                                                         % number of interior breakpoints
ws = warning('off', 'MATLAB:rankDeficientMatrix');                     % coinciding breakpoints make the basis rank-deficient
wsS = warning('off', 'MATLAB:singularMatrix');                         % (the SSE is still correct; silence the solver)
restoreW = onCleanup(@() warning([ws wsS]));                           % restore warning states on exit
cost = @(bb) fitSSE(bb, x, y, lo, hi);                                 % error for a breakpoint vector
opts = optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000, 'TolX', 1e-3, 'TolFun', 1e-8);   % polish settings
best = inf; bestB = [];                                                % best fit so far
for s = 1:nStarts
    b0 = diffEvolution(cost, nb, lo, hi, stream);                      % global search (as scipy differential_evolution)
    [b, sse] = fminsearch(cost, b0, opts);                             % local polish
    if sse < best                                                      % keep the lowest error
        best = sse; bestB = b;
    end
end
breaks = [lo; sort(min(max(bestB, lo), hi)); hi];                      % full breakpoint list, as pwlf (ends included)
[~, yHat] = fitSSE(bestB, x, y, lo, hi);                               % fitted values at the best breakpoints
F = struct('breaks', breaks', 'yHat', yHat, 'sse', best);              % fit summary
end

function [sse, yHat] = fitSSE(b, x, y, lo, hi)
% Sum of squared errors of the continuous piecewise-linear fit for breakpoints b.
b = sort(min(max(b(:), lo), hi));                                      % keep breakpoints ordered and inside the range
A = [ones(size(x)), x - lo, max(0, x - b')];                           % pwlf basis: intercept, slope, hinge per breakpoint
coef = A \ y;                                                          % least squares (QR; basic solution if rank-deficient)
yHat = A * coef;                                                       % fitted values
sse = sum((y - yHat) .^ 2);                                            % error to minimise
end

function bBest = diffEvolution(cost, nb, lo, hi, stream)
% Differential evolution with scipy.optimize.differential_evolution defaults
% (strategy best1bin, popsize 15, mutation dithered in [0.5 1), recombination 0.7,
% maxiter 1000, tol 0.01), bounded to [lo hi] in every dimension.
nPop = 15 * nb;                                                        % population size (scipy: popsize x dimensions)
pop = lo + (hi - lo) * rand(stream, nPop, nb);                         % random initial population inside the bounds
E = zeros(nPop, 1);                                                    % cost of each member
for i = 1:nPop, E(i) = cost(pop(i, :)'); end                           % evaluate the initial population
for gen = 1:1000                                                       % scipy default maxiter
    [~, iBest] = min(E);                                               % current best member
    F = 0.5 + 0.5 * rand(stream);                                      % dithered mutation factor for this generation
    for i = 1:nPop
        r = randperm(stream, nPop, 3);                                 % candidates for the difference vector
        r(r == i) = [];                                                % never use the member itself
        if numel(r) < 2, r = setdiff(randperm(stream, nPop, 4), i, 'stable'); end   % refill if needed
        mutant = pop(iBest, :) + F * (pop(r(1), :) - pop(r(2), :));    % best1: best + F x (difference of two members)
        cross = rand(stream, 1, nb) < 0.7;                             % binomial crossover mask (recombination 0.7)
        cross(randi(stream, nb)) = true;                               % at least one parameter from the mutant
        trial = pop(i, :);                                             % start from the current member
        trial(cross) = mutant(cross);                                  % take the crossed parameters from the mutant
        out = trial < lo | trial > hi;                                 % parameters outside the bounds
        trial(out) = lo + (hi - lo) * rand(stream, 1, sum(out));       % scipy: re-draw them uniformly inside
        eT = cost(trial');                                             % cost of the trial
        if eT <= E(i)                                                  % greedy selection
            pop(i, :) = trial; E(i) = eT;                              % replace the member
            if eT < E(iBest), iBest = i; end                           % update the best immediately (scipy 'immediate')
        end
    end
    if std(E) <= 0.01 * abs(mean(E))                                   % scipy convergence: std <= atol + tol x |mean|
        break
    end
end
[~, iBest] = min(E);                                                   % final best member
bBest = sort(pop(iBest, :))';                                          % breakpoints, ascending column
end

% =========================================================================
function plotQC(R, name, saveDir)
% QC figures matching notebook panels A (band profiles + fits), B (PSD), C (correlation).
bandCol = [0 0 0; 0 0 0.55; 0 0 1; 0.39 0.58 0.93; 0.68 0.85 0.9];    % notebook band colours
nB = numel(R.bandBoundaryUm);                                          % number of bands
y = R.yUm;                                                             % site depths (um from tip)

fA = figure('Position', [100 100 900 600]); hold on;                   % A: band power vs depth
patch([-1e3 1e3 1e3 -1e3], [R.surfaceMeanUm R.surfaceMeanUm y(end) y(end)], [0.94 0.5 0.5], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % out-of-brain shading (above the surface)
patch([-1e3 1e3 1e3 -1e3], [y(1) y(1) R.surfaceMeanUm R.surfaceMeanUm], [0.68 0.85 0.9], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % in-brain shading
h = gobjects(1, nB);                                                   % legend handles
for b = 1:nB
    c = bandCol(min(b, size(bandCol, 1)), :);                          % band colour
    plot(R.bandPowerDb(:, b), y, '-', 'Color', [c 0.45], 'HandleVisibility', 'off');   % raw profile
    h(b) = plot(R.fits(b).yHat, y, '-', 'Color', c, 'LineWidth', 2, 'DisplayName', ...
        sprintf('%g-%g Hz (%.0f \\mum)', R.params.bands(b, :), R.bandBoundaryUm(b)));  % piecewise fit
    yline(R.bandBoundaryUm(b), '--', 'Color', c, 'HandleVisibility', 'off');            % band boundary
end
yline(R.surfaceMeanUm, 'r', 'LineWidth', 2.4, 'DisplayName', sprintf('Mean boundary: %.0f \\mum', R.surfaceMeanUm));
xlim([min(R.bandPowerDb(:)) max(R.bandPowerDb(:))] + [-1 1]);          % power range
ylim([y(1) y(end)]);                                                   % column range
xlabel('Power (dB)'); ylabel('Distance to probe tip (\mum)');          % axis labels with units
title('Power per frequency band with piecewise brain boundary', 'Interpreter', 'none');
legend('Location', 'eastoutside'); grid on; box off;                   % legend outside the data

fB = figure('Position', [100 100 900 600]);                            % B: PSD vs depth
keep = R.freqs <= 100;                                                 % frequencies shown (Hz)
imagesc(R.freqs(keep), y, R.psdDb(:, keep)); axis xy; colormap(fB, jet);   % depth x frequency map
cb = colorbar; cb.Label.String = 'Power (dB)';                         % colour scale with units
hold on; yline(R.surfaceMeanUm, 'r', 'LineWidth', 1.7);                % surface estimate
for b = 1:nB, yline(R.bandBoundaryUm(b), '--', 'Color', [0.5 0.5 0.5]); end   % band boundaries
xlabel('Frequency (Hz)'); ylabel('Distance to probe tip (\mum)'); title('Power spectrum');

fC = figure('Position', [100 100 700 600]);                            % C: LFP correlation
imagesc(y, y, R.corr, [-0.6 1]); axis xy; colormap(fC, jet);           % site x site correlation
cb = colorbar; cb.Label.String = 'Correlation (r)';                    % colour scale
hold on; yline(R.surfaceMeanUm, 'r', 'LineWidth', 1.5);                % surface estimate
for b = 1:nB, yline(R.bandBoundaryUm(b), '--', 'Color', [0.5 0.5 0.5]); end   % band boundaries
xlabel('Distance to probe tip (\mum)'); ylabel('Distance to probe tip (\mum)'); title('LFP correlation');

if strlength(saveDir) > 0                                              % save the three panels
    if ~isfolder(saveDir), mkdir(saveDir); end                         % create the output folder
    tags = ["A_bandPower", "B_psd", "C_corr"];                         % panel names
    figs = [fA, fB, fC];                                               % panel handles
    for i = 1:3
        exportgraphics(figs(i), fullfile(saveDir, name + "_lfpSurface_" + tags(i) + ".pdf"), 'ContentType', 'vector');
    end
end
end
