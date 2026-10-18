function p = segSignFlipNull(DSeg, tc, nCats, nBoot)
% segSignFlipNull  Null-simulation mirror of the segmented max-statistic test.
%
%   p = segSignFlipNull(DSeg, tc, nCats, nBoot)
%
%   Used by checkMovingWindowNullFPR and checkLayoutNullFPR; one segment (nSegs = 1)
%   reduces it to the fixed-window max-category test.
% Segmented max-statistic sign-flip test of StatisticsPerNeuron (segMode branch,
% v6): observed = max over categories x segments of the category means; null =
% the same max after flipping each TRIAL's sign (shared across segments), with
% a running max over segments. DSeg: [nTrials x nNeurons x nSegs], rows in
% contiguous category blocks.
[nTr, nU, nSegs] = size(DSeg);                                     % dimensions
signs   = 2 * randi(2, nTr, nBoot) - 3;                            % one sign per trial and permutation
signsR  = reshape(signs, tc, nCats, nBoot);                        % [trialsCat x nCats x nBoot]
R       = reshape(DSeg, tc, nCats, nU, nSegs);                     % [trialsCat x nCats x nNeurons x nSegs]
cm      = reshape(mean(R, 1), nCats, nU, nSegs);                   % observed category means per segment
obs     = max(reshape(permute(cm, [1 3 2]), nCats*nSegs, nU), [], 1);  % max over cats x segs, one neuron per column
nullMax = -Inf(nBoot, nU);                                         % running max over segments
for seg = 1:nSegs
    cmP = pagemtimes(permute(R(:,:,:,seg), [3 1 2]), permute(signsR, [1 3 2])) / tc;  % [nNeurons x nBoot x nCats]
    mx  = max(permute(cmP, [3 1 2]), [], 1);                       % max over categories [1 x nNeurons x nBoot]
    nullMax = max(nullMax, reshape(mx, nU, nBoot).');              % fixed (per-neuron) reshape, then running max
end
p = mean(nullMax >= obs, 1);                                       % right-tailed p per neuron
end
