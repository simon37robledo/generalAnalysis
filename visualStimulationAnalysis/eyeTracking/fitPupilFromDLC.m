function T = fitPupilFromDLC(experimentNum, params)
% FITPUPILFROMDLC  Per-frame pupil ellipse from the DLC pupil points of one insertion.
%
%   T = fitPupilFromDLC(experimentNum)           % fit and write <outDir>\eye_<expNum>_<animal>_ins<N>.csv
%   T = fitPupilFromDLC(experimentNum, write=false)
%
%   Same rule as the original eye_tracking_analysis notebook (and Elena's threshold): the 12
%   Pupil_* points with DLC likelihood >= thr, at least minPoints of them, least-squares ellipse
%   (Halir & Flusser 1998, numerically stable form of Fitzgibbon et al. 1999; the method of the
%   Python LsqEllipse package). Coordinates are raw camera pixels of the video DLC ran on.
%   One row per DLC row (= saved video frame), frame_index 0-based, failed fits kept as NaN.
%
%   fit_valid: ellipse found, >= minPoints points, minor/major >= minAxisRatio, mean radial fit
%   error <= maxRelError x semi-major axis, and major diameter within [0.5, 2] x the session median.
%   Columns follow Elena's tables, so getEyePosition reads either (tableDir parameter).

arguments
    experimentNum (1,1) double                                                     % Excel row
    params.manifest string = "W:\Large_scale_mapping_NP\EyeTracking\manifest.csv"  % DLC file per insertion
    params.outDir string = "W:\Large_scale_mapping_NP\EyeTracking\tables"          % where our tables go
    params.thr double = 0.5                                                        % DLC likelihood threshold (Elena's)
    params.minPoints double = 6                                                    % points needed for a fit (notebook: > 5)
    params.minAxisRatio double = 0.5                                               % minor/major below this = implausible pupil
    params.maxRelError double = 0.1                                                % mean radial error / semi-major above this = bad fit
    params.write logical = true                                                    % write the CSV
end

M = readtable(params.manifest, 'TextType', 'string');                              % manifest
r = M(M.experimentNum == experimentNum, :);                                        % this insertion
assert(height(r) == 1, 'fitPupilFromDLC:row', 'experimentNum %d: %d manifest rows (need 1).', experimentNum, height(r));
assert(r.status_elena == "ready", 'fitPupilFromDLC:notReady', 'experimentNum %d: manifest status %s.', experimentNum, r.status_elena);   % has a DLC file

%% Read the pupil columns of the DLC CSV
fid = fopen(r.dlc_csv, 'r');                                                       % open
assert(fid > 0, 'fitPupilFromDLC:open', 'Cannot open %s', r.dlc_csv);
cleaner = onCleanup(@() fclose(fid));                                              % always close
fgetl(fid);                                                                        % line 1: scorer
parts = strtrim(split(string(fgetl(fid)), ","));                                   % line 2: bodyparts
coords = strtrim(split(string(fgetl(fid)), ","));                                  % line 3: x / y / likelihood
isPupil = startsWith(parts, "Pupil_");                                             % pupil point columns
fmtCells = repmat("%*f", 1, numel(parts));                                        % skip every column ...
fmtCells(isPupil) = "%f";                                                          % ... except the pupil columns
C = textscan(fid, char(strjoin(fmtCells, "")), 'Delimiter', ',', 'CollectOutput', true);   % read all rows
P = C{1};                                                                          % [frames x pupil columns]
pn = parts(isPupil); pc = coords(isPupil);                                         % names of the kept columns
names = unique(pn, 'stable');                                                      % the pupil points
nP = numel(names); nF = size(P, 1);                                                % points, frames
assert(nF == r.n_dlc_frames, 'fitPupilFromDLC:rows', 'experimentNum %d: read %d rows, manifest %d.', experimentNum, nF, r.n_dlc_frames);   % complete read
X = nan(nF, nP); Y = X; L = X;                                                     % x, y, likelihood per point
for j = 1:nP
    X(:, j) = P(:, pn == names(j) & pc == "x");                                    % x
    Y(:, j) = P(:, pn == names(j) & pc == "y");                                    % y
    L(:, j) = P(:, pn == names(j) & pc == "likelihood");                           % likelihood
end
clear P C                                                                          % free memory

%% Fit each frame
cx = nan(nF, 1); cy = cx; maj = cx; mnr = cx; phi = cx; err = cx;                  % ellipse parameters per frame
nUsed = sum(L >= params.thr, 2);                                                   % confident points per frame
for f = find(nUsed >= params.minPoints)'                                           % frames with enough points
    k = L(f, :) >= params.thr;                                                     % confident points
    [cx(f), cy(f), maj(f), mnr(f), phi(f), err(f)] = fitEllipse(X(f, k)', Y(f, k)');   % least-squares ellipse
end
medMaj = median(maj, 'omitnan');                                                   % session pupil size
valid = nUsed >= params.minPoints & isfinite(cx) & mnr ./ maj >= params.minAxisRatio ...
    & err ./ (maj / 2) <= params.maxRelError & maj >= 0.5 * medMaj & maj <= 2 * medMaj;   % plausible fits

T = table(repmat(experimentNum, nF, 1), repmat(r.animal, nF, 1), repmat(r.insertion, nF, 1), (0:nF-1)', ...
    cx, cy, maj, mnr, phi, err, nUsed, valid, repmat(params.thr, nF, 1), ...
    'VariableNames', {'experiment_num', 'animal', 'insertion', 'frame_index', 'center_x_px', 'center_y_px', ...
    'major_diameter_px', 'minor_diameter_px', 'major_axis_orientation_deg', 'mean_radial_fit_error_px', ...
    'n_pupil_points', 'fit_valid', 'pupil_likelihood_threshold'});                 % Elena's column names
if params.write
    if ~isfolder(params.outDir), mkdir(params.outDir); end                         % output folder
    writetable(T, fullfile(params.outDir, r.out_table));                           % one CSV per insertion
end
fprintf('experimentNum %d (%s ins %d): %d frames, %.1f%% valid, median pupil %.1f px.\n', ...
    experimentNum, r.animal, r.insertion, nF, 100 * mean(valid), medMaj);
end

function [x0, y0, major, minor, phiDeg, err] = fitEllipse(x, y)
% Halir & Flusser (1998) direct least-squares ellipse fit; returns centre, full axes, orientation,
% and the mean radial distance of the points from the fitted ellipse (px).
x0 = NaN; y0 = NaN; major = NaN; minor = NaN; phiDeg = NaN; err = NaN;             % default: no ellipse
mx = mean(x); my = mean(y); s = max(std(x), std(y));                               % normalise for numerical stability
if s == 0, return, end                                                             % degenerate points
u = (x - mx) / s; v = (y - my) / s;                                                % normalised coordinates
D1 = [u.^2, u.*v, v.^2]; D2 = [u, v, ones(size(u))];                               % quadratic / linear parts
S1 = D1' * D1; S2 = D1' * D2; S3 = D2' * D2;                                       % scatter blocks
if rcond(S3) < 1e-12, return, end                                                  % collinear points
Tm = -S3 \ S2';                                                                    % linear part as a function of the quadratic part
Mm = S1 + S2 * Tm;                                                                 % reduced scatter matrix
Mm = [Mm(3, :) / 2; -Mm(2, :); Mm(1, :) / 2];                                      % premultiply by inv(constraint matrix)
[V, ~] = eig(Mm);                                                                  % eigenvectors
cond4 = 4 * V(1, :) .* V(3, :) - V(2, :).^2;                                       % ellipse condition 4AC - B^2 > 0
a1 = V(:, cond4 > 0);                                                              % the ellipse solution
if size(a1, 2) ~= 1, return, end                                                   % none (or ambiguous)
a = real([a1; Tm * a1]);                                                           % conic A..F (normalised frame)
A = a(1); B = a(2); Cc = a(3); D = a(4); E = a(5); F = a(6);                       % names
den = B^2 - 4 * A * Cc;                                                            % < 0 for an ellipse
uc = (2 * Cc * D - B * E) / den; vc = (2 * A * E - B * D) / den;                   % centre (normalised)
F0 = A * uc^2 + B * uc * vc + Cc * vc^2 + D * uc + E * vc + F;                     % conic value at the centre
ev = eig([A B/2; B/2 Cc]);                                                         % quadratic-form eigenvalues
ax = sqrt(-F0 ./ ev);                                                              % semi-axes (normalised)
if any(~isfinite(ax)) || any(imag(ax) ~= 0) || any(ax <= 0), return, end           % not a real ellipse
th = 0.5 * atan2(B, A - Cc);                                                       % orientation of the first axis
x0 = mx + s * uc; y0 = my + s * vc;                                                % centre [px]
major = 2 * s * max(ax); minor = 2 * s * min(ax);                                  % full axes [px]
[~, im] = max(ax); phiDeg = rad2deg(th) + 90 * (im == 2);                          % orientation of the major axis [deg]
du = u - uc; dv = v - vc;                                                          % points relative to the centre
q = (A * du.^2 + B * du .* dv + Cc * dv.^2) / (-F0);                               % 1 on the ellipse
err = mean(abs(sqrt(max(q, 0)) - 1) .* hypot(du, dv) ./ max(sqrt(max(q, 0)), eps)) * s;   % mean radial distance [px]
end
