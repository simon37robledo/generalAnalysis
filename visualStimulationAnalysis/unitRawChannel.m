function ch = unitRawChannel(sortFolder, phyIDs)
% unitRawChannel  Recording channel (NP.getData numbering, 1-based) of phy units.
%
%   ch = unitRawChannel(obj.spikeSortingFolder, phyIDs)
%
%   Kilosort's cluster_info.tsv 'ch' is a 0-based index into the channels Kilosort
%   kept (channel_map.npy), not the recording channel. convertPhySorting2tIc's
%   ic(1,:) copies it, adding 1 only when some unit sits on channel 0, so using
%   ic(1,:) to read raw data picked a neighbouring channel (or, on 4-shank probes
%   with a non-identity channel map, a channel on another shank: up to 48 channels
%   off in recording 49). Here: ch = channel_map(ks_ch + 1) + 1.
%   Older sessions keep the phy files one folder above an empty kilosort/ folder.

sortDir = char(sortFolder);                                        % phy output folder
if ~isfile(fullfile(sortDir, 'cluster_info.tsv'))                  % older sessions: phy files one level up
    sortDir = fileparts(sortDir);
end
assert(isfile(fullfile(sortDir, 'cluster_info.tsv')) && isfile(fullfile(sortDir, 'channel_map.npy')), ...
    'unitRawChannel: no cluster_info.tsv / channel_map.npy in %s or its parent.', sortFolder);   % fail loudly
ci = readtable(fullfile(sortDir, 'cluster_info.tsv'), 'FileType', 'delimitedtext');   % phy cluster table
cm = double(readNPY(fullfile(sortDir, 'channel_map.npy')));        % kept channels -> recording channel (0-based)
[ok, row] = ismember(phyIDs(:), ci.cluster_id);                    % each unit's row
assert(all(ok), 'unitRawChannel: %d phy IDs not in cluster_info.tsv.', sum(~ok));   % unit alignment
ks = ci.ch(row);                                                   % Kilosort channel (0-based, kept channels)
assert(all(ks >= 0 & ks < numel(cm)), 'unitRawChannel: Kilosort channel outside channel_map.');
ch = cm(ks + 1) + 1;                                               % recording channel, 1-based (NP.getData)
end
