function T = getUnitDepths(exList, params)
% getUnitDepths  Depth (Z, um below the brain surface) of every good unit.
%
%   T = getUnitDepths(exList)                     % source "coorZ" (default)
%   T = getUnitDepths(exList, source="lfp")       % LFP-estimated brain surface
%
%   Unit position on the probe: the Kilosort channel of the unit (cluster_info.tsv
%   'ch', 0-based index into the kept channels) looked up in channel_positions.npy
%   of the same sorting folder. This avoids convertPhySorting2tIc's ic(1,:), which is
%   shifted to 1-based only when some unit sits on channel 0 (so most recordings are
%   off by one channel there).
%
%   source = "coorZ": depth = coor_Z (Experiment_Excel.xlsx, as getNeuronDepths)
%                     minus the unit's y position on the probe.
%   source = "lfp"  : depth = brain-surface y position estimated from the LFP
%                     (VStimAnalysis.getLfpSurface: coor_Z plus the animal's LFP offset
%                     from its depth tests) minus the unit's y position. Animals without
%                     a valid test keep coor_Z; their units get source "coorZ-noTest".
%
%   Output table: exp, animal, phyID, unitIdx (order of good units, as the stats),
%   chKS (Kilosort channel, 0-based), yUm (probe position), depthUm, source (per unit),
%   depthValid (false for every unit of an insertion with any unit above the surface:
%   its surface reference is wrong; rank tests within the insertion are still valid).

arguments
    exList double                                                  % recordings (master Excel rows)
    params.source (1,1) string {mustBeMember(params.source, ["coorZ", "lfp"])} = "coorZ"   % depth reference
    params.excelPath char = '\\sil3\data\Large_scale_mapping_NP\Experiment_Excel.xlsx'   % master table (coor_Z)
end

if params.source == "coorZ"
    X = readtable(params.excelPath);                               % master table with coor_Z per recording
end
T = table();
for ex = exList
    NP  = loadNPclassFromTable(ex);                                % recording
    obj = linearlyMovingBallAnalysis(NP, MultipleOffsets=false);   % any stimulus object: shared sorting folder (any MB session)
    p   = NP.convertPhySorting2tIc(obj.spikeSortingFolder);        % good-unit order used by every statistic
    good  = string(p.label') == 'good';                            % good somatic units
    phyID = p.phy_ID(good);                                        % their phy cluster IDs, in stats order

    sortDir = obj.spikeSortingFolder;                              % phy output folder
    if ~isfile(fullfile(sortDir, 'cluster_info.tsv'))              % older sessions keep phy files one level up
        sortDir = fileparts(sortDir);                              % (their kilosort/ subfolder is empty)
    end
    assert(isfile(fullfile(sortDir, 'cluster_info.tsv')) && isfile(fullfile(sortDir, 'channel_positions.npy')), ...
        'getUnitDepths: exp %d: no cluster_info.tsv / channel_positions.npy in %s or its parent.', ex, obj.spikeSortingFolder);
    ci   = readtable(fullfile(sortDir, 'cluster_info.tsv'), 'FileType', 'delimitedtext');   % phy cluster table
    pos  = readNPY(fullfile(sortDir, 'channel_positions.npy'));     % [nKeptChannels x 2] x, y (um)
    [ok, row] = ismember(phyID, ci.cluster_id);                    % each good unit's row in cluster_info
    assert(all(ok), 'getUnitDepths: exp %d: %d good units missing from cluster_info.tsv.', ex, sum(~ok));   % alignment
    ch = ci.ch(row);                                               % Kilosort channel (0-based)
    assert(all(ch >= 0 & ch < size(pos, 1)), 'getUnitDepths: exp %d: channel index outside channel_positions.', ex);
    y  = pos(ch + 1, 2);                                           % unit y position on the probe (um)

    switch params.source
        case "coorZ"
            surfaceY = X.coor_Z(ex);                               % probe y of the surface = manipulator depth (um)
            depth = surfaceY - y;                                  % as getNeuronDepths (coor_Z minus probe y)
            src = "coorZ";                                         % depth source of these units
        case "lfp"
            [surfaceY, src] = readLfpSurface(ex, obj);             % probe y of the brain surface (um) and its source
            depth = surfaceY - y;                                  % distance below the surface
    end
    animal = string(regexp(NP.recordingName, '(PV|SA)\d+', 'match', 'once'));   % animal ID
    valid = all(depth >= 0);                                       % a unit above the surface means the surface reference is wrong
    if ~valid                                                      % e.g. recs 1, 8: coor_Z ~ -21 (not a depth)
        warning('getUnitDepths:negativeDepth', ...
            'Exp %d (%s): surface %.1f um (%s) gives depths down to %.0f um (above the surface); depthValid = false.', ...
            ex, animal, surfaceY, src, min(depth));                 % explicit: absolute depths of this insertion are unusable
    end
    n = numel(phyID);
    T = [T; table(repmat(ex, n, 1), repmat(animal, n, 1), phyID(:), (1:n)', ch(:), y(:), depth(:), ...
        repmat(src, n, 1), repmat(valid, n, 1), ...
        'VariableNames', {'exp','animal','phyID','unitIdx','chKS','yUm','depthUm','source','depthValid'})];   %#ok<AGROW>
end
end

% =========================================================================
function [surfaceY, src] = readLfpSurface(ex, obj)
% readLfpSurface  Brain-surface position (probe y, um; frame of channel_positions.npy)
%   from the LFP depth estimate (VStimAnalysis.getLfpSurface: coor_Z calibrated by the
%   animal's whole-shank / depth tests). src is "lfp-offset", or "coorZ-noTest" when the
%   animal has no valid test (getLfpSurface warns).
S = obj.getLfpSurface(verbose=false);                              % surface for this insertion
assert(S.exp == ex, 'getUnitDepths:lfpRow', ...                    % the object must be this Excel row
    'Exp %d: getLfpSurface matched Excel row %d.', ex, S.exp);
surfaceY = S.surfaceY;                                             % probe y of the brain surface (um)
src = S.source;                                                    % depth source of these units
end
