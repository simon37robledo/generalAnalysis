function T = windowMatchedDelta(vsObj, stim, params)
% windowMatchedDelta  Split-half delta firing rate on 300 ms segments, so that
%   stimuli tested with long windows can be compared with MB on the same time scale.
%
%   T = windowMatchedDelta(vsObj, "SDGm")    % StaticDriftingGratingAnalysis
%   T = windowMatchedDelta(vsObj, "MBR")     % linearlyMovingBarAnalysis
%
%   MB's DeltaRateCV is the split-half peak of a 300 ms window from the RF crossing.
%   The SDGm test window (1923 ms) and the MBR segments (stimDur/9) differ from that,
%   so their delta rates are not on the same scale. Here the response window is cut
%   into consecutive segments of segMs (300 ms) - SDGm: the drifting test window
%   from motion onset; MBR: the whole sweep - and each trial's segment rates minus
%   its blank pre-trial baseline rate enter splitHalfDelta: the best (direction x
%   segment) is chosen on odd repetitions and measured on even ones, and vice versa.
%   Units: spikes/s; ~0 for a unit with no response (cross-validated).
%   The responsiveness tests are unchanged; this is only the effect-size measure.
%
%   Output: table with phyID, unitIdx (stats order), dRate300 (spikes/s), nSeg.

arguments
    vsObj                                                          % stimulus analysis object
    stim (1,1) string {mustBeMember(stim, ["SDGm", "MBR"])}        % stimulus family
    params.segMs (1,1) double = 300                                % segment length (ms), as MB's window
end

S  = vsObj.StatisticsPerNeuron;                                    % figure-design statistics (window design)
rw = vsObj.ResponseWindow;                                         % condition matrix and timing
p  = vsObj.dataObj.convertPhySorting2tIc(vsObj.spikeSortingFolder);   % spike sorting
good  = string(p.label') == 'good';                                % good somatic units
goodU = p.ic(:, good);                                             % their spike indices
phyID = p.phy_ID(good);                                            % their phy IDs
cn = rw.colNames{1}(5:end);                                        % condition column names

if stim == "SDGm"
    C   = rw.C;                                                    % [onset angle tf sf] per trial
    on  = C(:, 1) + vsObj.VST.static_time * 1000;                  % motion onsets (ms)
    W   = S.params.BaseRespWindow;                                 % default window
    if isfield(S.params, 'SDGWindows') && ~isempty(S.params.SDGWindows)
        W = S.params.SDGWindows(2);                                % drifting test window (figure design)
    end
    W   = min(W, floor(rw.Moving.stimDur));                        % never past the drifting phase
    bl  = min(W, rw.stimInter - 100);                              % blank baseline length, as the stats
    bOn = C(:, 1) - bl;                                            % baseline before the static onset
    ci  = find(strcmpi(cn, 'angles'), 1) + 1;                      % angle column
    dirs = mod(pi/2 - deg2rad(C(:, ci)), 2*pi);                    % drift direction (360 = 0 merged)
else
    C   = rw.Speed1.C;                                             % [onset direction ...] per trial
    on  = C(:, 1);                                                 % sweep onsets (ms)
    W   = floor(rw.Speed1.stimDur);                                % whole sweep (ms)
    bl  = min(params.segMs, rw.stimInter - 100);                   % baseline as long as one segment
    bOn = on - bl;                                                 % baseline before onset
    ci  = find(strcmpi(cn, 'directions'), 1) + 1;                  % direction column
    dirs = mod(C(:, ci), 2*pi);                                    % bar direction
end
nSeg = floor(W / params.segMs);                                    % whole 300 ms segments in the window
assert(nSeg >= 1, 'windowMatchedDelta: %s window (%d ms) shorter than one segment.', stim, W);

Mr = BuildBurstMatrix(goodU, round(p.t), round(on(:)'), nSeg * params.segMs);   % [trials x units x ms]
Mb = BuildBurstMatrix(goodU, round(p.t), round(bOn(:)'), bl);                    % [trials x units x ms]
[nTr, nU, ~] = size(Mr);
seg = squeeze(mean(reshape(Mr, nTr, nU, params.segMs, nSeg), 3));   % [trials x units x segments] spikes/ms
if nU == 1, seg = reshape(seg, nTr, 1, nSeg); end                  % keep the unit dimension
if nSeg == 1, seg = reshape(seg, nTr, nU, 1); end                  % keep the segment dimension
X  = seg - mean(Mb, 3);                                            % minus each trial's baseline rate (spikes/ms)

[~, ~, grp] = uniquetol(dirs, 1e-4);                               % direction group per trial
half = zeros(nTr, 1);                                              % split: alternate repetitions per direction
for g = unique(grp)'
    ii = find(grp == g);                                           % this direction's trials, in time order
    half(ii) = 2 - mod(1:numel(ii), 2)';                           % 1, 2, 1, 2, ...
end
assert(all(ismember(half, [1 2])), 'windowMatchedDelta: split failed.');   % every trial assigned
dRate = splitHalfDelta(X, grp, half);                              % [1 x units] spikes/s

T = table(phyID(:), (1:nU)', dRate(:), repmat(nSeg, nU, 1), 'VariableNames', {'phyID','unitIdx','dRate300','nSeg'});
end
