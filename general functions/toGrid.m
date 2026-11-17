function G = toGrid(A, nGrid)
% toGrid  Block-average the first two (square) dimensions of A to nGrid x nGrid;
% trailing dimensions are kept. Copy of the local toGrid in SpatialTuningIndex.m
% (the 9x9 map the stripe test uses), shared here for the RF display code.
%
%   A      [side, side, ...] map(s); side must be divisible by nGrid
%   nGrid  number of grid cells per side (9 for the MB and SB stimulus grids)
%   G      [nGrid, nGrid, ...] mean of each side/nGrid x side/nGrid block
sz = size(A);                                                        % [rfY, rfX, ...]
bs = sz(1) / nGrid;                                                  % pixels per grid cell
assert(sz(1) == sz(2) && bs == round(bs), 'toGrid: %dx%d map is not divisible into %dx%d cells.', sz(1), sz(2), nGrid, nGrid);  % square, whole blocks only
rest = sz(3:end);                                                    % trailing dimensions
G = reshape(A, [bs, nGrid, bs, nGrid, prod(rest)]);                 % split rows and columns into blocks
G = mean(mean(G, 1), 3);                                             % average within each block
G = reshape(G, [nGrid, nGrid, rest, 1]);                             % [nGrid, nGrid, ...]
end
