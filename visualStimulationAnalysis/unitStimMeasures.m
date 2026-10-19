function M = unitStimMeasures(NP, stim)
% unitStimMeasures  Per-unit responsiveness measures of one stimulus, from the
%   figure-design StatisticsPerNeuron cache (bare call: errors if stale).
%
%   M = unitStimMeasures(NP, "MB")     % stim: MB, RG, SDGm, SDGs, MBR, FFF
%
%   One row per good unit, in the order every statistic uses:
%   phyID, unitIdx, p, responsive (p < 0.05), present, z (ZScoreU),
%   dRate (DeltaRateCV, spikes/s). MB takes, per unit, the speed with the lowest p
%   (Figure 1 rule) and requires presence during every speed's block.
%   Returns an empty table when the stimulus was not recorded.

ctor = struct('MB', @linearlyMovingBallAnalysis, 'RG', @rectGridAnalysis, 'SDGm', @StaticDriftingGratingAnalysis, ...
    'SDGs', @StaticDriftingGratingAnalysis, 'MBR', @linearlyMovingBarAnalysis, 'FFF', @fullFieldFlashAnalysis);   % analysis classes
M = table();
try
    o = ctor.(char(stim))(NP);                                     % stimulus object (errors if not recorded)
    o.ResponseWindow;                                              % cached response window must exist
catch
    return                                                         % stimulus absent from this recording
end
S = o.StatisticsPerNeuron;                                         % figure-design statistics
switch stim
    case "MB"
        f = fieldnames(S); f = f(startsWith(string(f), "Speed"));  % speed sub-fields
        P = cell2mat(cellfun(@(x) S.(x).pvalsResponse(:), f', 'UniformOutput', false));   % [units x speeds]
        Z = cell2mat(cellfun(@(x) S.(x).ZScoreU(:),       f', 'UniformOutput', false));
        D = cell2mat(cellfun(@(x) S.(x).DeltaRateCV(:),   f', 'UniformOutput', false));
        [p, b] = min(P, [], 2); li = sub2ind(size(P), (1:size(P,1))', b);   % lowest-p speed per unit
        z = Z(li); dr = D(li); pres = statPresence(S, "MB", "");   % values at that speed; presence in all speeds
    case {"SDGm", "SDGs"}
        fn = 'Moving'; if stim == "SDGs", fn = 'Static'; end       % grating phase
        p = S.(fn).pvalsResponse(:); z = S.(fn).ZScoreU(:); dr = S.(fn).DeltaRateCV(:);
        pres = logical(S.(fn).present(:));
    case "MBR"
        p = S.Speed1.pvalsResponse(:); z = S.Speed1.ZScoreU(:); dr = S.Speed1.DeltaRateCV(:);
        pres = statPresence(S, "MBR", "");
    otherwise                                                      % flat structs (RG, FFF)
        p = S.pvalsResponse(:); z = S.ZScoreU(:); dr = S.DeltaRateCV(:);
        pres = statPresence(S, string(stim), "");
end
sp = NP.convertPhySorting2tIc(o.spikeSortingFolder);               % unit order used by the statistics
phyID = sp.phy_ID(string(sp.label') == 'good');                    % good units
n = numel(phyID);
assert(numel(p) == n && numel(pres) == n, 'unitStimMeasures: %s statistics cover %d units, sorting has %d.', stim, numel(p), n);
M = table(phyID(:), (1:n)', p, p < 0.05, pres(:), z, dr, repmat(string(stim), n, 1), ...
    'VariableNames', {'phyID','unitIdx','p','responsive','present','z','dRate','stim'});
end
