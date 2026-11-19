function T = getUnitAtlasPositions(exList, params)
% getUnitAtlasPositions  Atlas position, brain region and depth of every good unit.
%
%   T = getUnitAtlasPositions(exList)
%
%   Source: the per-recording channel tables of the brain-clearing atlas alignment
%   (W:\Large_scale_mapping_NP\atlas_alignments\recording_tables\<animal>\*_channel_positions.csv).
%   One row per recorded channel (384). A unit takes the row of its Kilosort channel
%   (cluster_info.tsv 'ch', 0-based row of channel_positions.npy; same units and order as getUnitDepths).
%
%   Channels are paired by probe geometry (row order differs between the two files on 4-shank
%   probes): table contact (x, y) = Kilosort contact (x, y + one bank offset). Hard errors if the
%   row counts differ, a Kilosort contact has no table contact, or the pairing is not 1-to-1.
%
%   Output: exp, animal, phyID, penetration (table 'insertion' = physical penetration, unique
%   within animal), penID (globally unique code: animal number*100 + penetration), chKS, shank,
%   depthAtlasUm (um below the brain surface along the shank; NaN if undetermined),
%   regionId, region, regionKnown (false for ID 0 "Unlabeled atlas" and -1 "outside atlas"),
%   atlasX/Y/Z_um and scanX/Y/Z_um (TIFF-local axes, NOT stereotaxic AP/ML/DV).

arguments
    exList double                                                  % recordings (master Excel rows)
    params.tableDir char = fullfile('W:\Large_scale_mapping_NP', 'atlas_alignments', 'recording_tables')   % atlas channel tables
end

T = table();                                                       % output, one block of rows per recording
for ex = exList
    NP  = loadNPclassFromTable(ex);                                % recording
    obj = linearlyMovingBallAnalysis(NP, MultipleOffsets=false);   % any stimulus object: shared sorting folder
    p   = NP.convertPhySorting2tIc(obj.spikeSortingFolder);        % good-unit order used by every statistic
    good  = string(p.label') == 'good';                            % good somatic units
    phyID = p.phy_ID(good);                                        % their phy cluster IDs, in stats order

    sortDir = obj.spikeSortingFolder;                              % phy output folder
    if ~isfile(fullfile(sortDir, 'cluster_info.tsv'))              % older sessions keep phy files one level up
        sortDir = fileparts(sortDir);                              % (their kilosort/ subfolder is empty)
    end
    assert(isfile(fullfile(sortDir, 'cluster_info.tsv')) && isfile(fullfile(sortDir, 'channel_positions.npy')), ...
        'getUnitAtlasPositions: exp %d: no cluster_info.tsv / channel_positions.npy in %s or its parent.', ex, obj.spikeSortingFolder);
    ci  = readtable(fullfile(sortDir, 'cluster_info.tsv'), 'FileType', 'delimitedtext');   % phy cluster table
    pos = readNPY(fullfile(sortDir, 'channel_positions.npy'));      % [nKeptChannels x 2] Kilosort x, y (um)
    [ok, row] = ismember(phyID, ci.cluster_id);                    % each good unit's row in cluster_info
    assert(all(ok), 'getUnitAtlasPositions: exp %d: %d good units missing from cluster_info.tsv.', ex, sum(~ok));   % alignment
    ch = ci.ch(row);                                               % Kilosort channel (0-based)
    assert(all(ch >= 0 & ch < size(pos, 1)), 'getUnitAtlasPositions: exp %d: channel index outside channel_positions.', ex);

    animal = string(regexp(NP.recordingName, '(PV|SA)\d+', 'match', 'once'));   % animal ID
    f = dir(fullfile(params.tableDir, char(animal), [NP.recordingName '_g0_tcat.imec0.ap_channel_positions*.csv']));   % atlas table of this recording
    assert(numel(f) == 1, 'getUnitAtlasPositions: exp %d (%s): %d atlas tables found (need exactly 1; a name conflict gives two).', ...
        ex, NP.recordingName, numel(f));                           % none or several tables: fail loudly
    t = readtable(fullfile(f.folder, f.name), 'VariableNamingRule', 'preserve', 'TextType', 'string');   % channel table

    % ---- alignment invariants between the Kilosort channel map and the atlas table
    assert(height(t) == size(pos, 1), 'getUnitAtlasPositions: exp %d: table has %d channels, Kilosort %d.', ex, height(t), size(pos, 1));
    % Row order differs between the two files on 4-shank probes, so channels are paired by probe geometry:
    % the table's contact (x, y) equals the Kilosort contact (x, y + one bank offset), and is unique per contact.
    yOff = min(t.probe_y_um) - min(pos(:, 2));                     % bank offset between table y and Kilosort y (um)
    keyKS  = round([pos(:, 1), pos(:, 2) + yOff] * 10);            % Kilosort contacts in table coordinates (0.1 um grid)
    keyTbl = round([t.probe_x_um(:), t.probe_y_um(:)] * 10);       % table contacts (0.1 um grid)
    assert(size(unique(keyTbl, 'rows'), 1) == height(t), 'getUnitAtlasPositions: exp %d: table contacts are not unique.', ex);   % one row per contact
    [okG, rowOfKS] = ismember(keyKS, keyTbl, 'rows');              % table row of every Kilosort channel
    assert(all(okG), 'getUnitAtlasPositions: exp %d: %d Kilosort contacts have no atlas-table contact (offset %.1f um).', ex, sum(~okG), yOff);
    assert(numel(unique(rowOfKS)) == numel(rowOfKS), 'getUnitAtlasPositions: exp %d: two Kilosort channels map to one table contact.', ex);   % 1-to-1

    r = rowOfKS(ch + 1);                                           % table row of each unit (via its Kilosort channel)
    depth = t.depth_from_brain_surface_um(r);                      % um below the surface along the shank
    nNaN = sum(isnan(depth));                                      % units without a depth
    if nNaN > 0                                                    % explicit warning
        warning('getUnitAtlasPositions:noDepth', 'Exp %d (%s): %d of %d units have no atlas depth (outside the brain-scan mask).', ...
            ex, animal, nNaN, numel(r));
    end
    if any(depth < 0)                                              % units above the entry surface (axial depth, may be small)
        warning('getUnitAtlasPositions:aboveSurface', 'Exp %d (%s): %d units lie above the brain entry surface (min %.0f um).', ...
            ex, animal, sum(depth < 0), min(depth));
    end
    rid = double(t.brain_region_id(r));                            % atlas region ID (0 unlabelled, -1 outside atlas)
    n = numel(phyID);                                              % units of this recording
    pen = repmat(double(t.insertion(1)), n, 1);                    % physical penetration within the animal
    assert(all(t.insertion == t.insertion(1)), 'getUnitAtlasPositions: exp %d: several penetrations in one table.', ex);   % one table = one penetration
    animalNum = str2double(extract(animal, digitsPattern));        % number in the animal ID (PV140 -> 140)
    Tk = table(repmat(ex, n, 1), repmat(animal, n, 1), phyID(:), pen, animalNum * 100 + pen, ch(:), double(t.shank(r)), depth(:), ...
        rid(:), t.brain_region(r), rid(:) > 0, ...
        t.atlas_x_um(r), t.atlas_y_um(r), t.atlas_z_um(r), t.brain_scan_x_um(r), t.brain_scan_y_um(r), t.brain_scan_z_um(r), ...
        'VariableNames', {'exp','animal','phyID','penetration','penID','chKS','shank','depthAtlasUm', ...
        'regionId','region','regionKnown','atlasX_um','atlasY_um','atlasZ_um','scanX_um','scanY_um','scanZ_um'});   % per-unit rows
    T = [T; Tk];                                                   %#ok<AGROW>
end
end
