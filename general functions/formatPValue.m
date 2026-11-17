function [txt, op] = formatPValue(p, pRes)
% FORMATPVALUE  Render a p-value as plain decimal text for figures and reports.
%
%   [txt, op] = formatPValue(p)
%   [txt, op] = formatPValue(p, pRes)
%
%   txt  - the number alone, as plain decimal text (no 'p', no operator).
%   op   - the relational operator to pair with it, '=' or '<'.
%   pRes - smallest p-value the estimator can resolve, i.e. 2/nBoot for a
%          two-tailed bootstrap p (the minimum is 2/nBoot, not 1/nBoot,
%          because the two-tailed value doubles the smaller tail). Sets the
%          floor used when p == 0. Defaults to 1e-3, a conservative value for
%          callers that do not know their nBoot.
%
%   Split this way so callers control spacing while never building a malformed
%   label: a cramped figure bracket writes ['p' op txt] -> 'p<0.001', console
%   output writes sprintf('p %s %s', op, txt) -> 'p < 0.001'.
%
%   Formatting rules:
%     * Always plain decimal - never scientific notation. 0.000300, not 3.00e-04.
%     * Three significant digits, trailing zeros kept (since 2026-11-13, so labels
%       look alike: 0.0160 next to 0.00942); powers of ten (the p == 0 floor) stay
%       short. Three rather than
%       two matters at the conventional thresholds: at two digits p=0.0499
%       would render as '0.05' and read as non-significant, and p=0.00099 as
%       '0.001'. Examples:
%         0.012   ->  '0.0120'
%         0.0499  ->  '0.0499'
%         3.0e-4  ->  '0.000300'
%         1.2e-5  ->  '0.0000120'
%     * A bootstrap p of exactly 0 means no resample crossed the null, so the
%       value is below the resolution of the run, not truly zero. It returns
%       op='<' with the floor set to the next power of ten at or above pRes:
%         nBoot = 1e5 -> pRes = 0.00002 -> '<0.0001'
%         nBoot = 1e4 -> pRes = 0.0002  -> '<0.001'
%       Deriving the floor this way keeps it consistent with the exact values
%       printed alongside it. A hard-coded '<0.001' would be incoherent next to
%       a resolved 'p=0.0003' on the same figure: it would assert the estimator
%       cannot see below 0.001 while the neighbouring label proves it can.
%       Rounding UP to a power of ten (rather than to pRes itself) keeps the
%       statement true: with p == 0 the run establishes only that the estimate
%       fell below its resolution. The raw 0 stays in the stats tables.
%     * NaN or empty returns txt='' and op='' so callers can draw nothing.
%
%   Deliberately no significance stars: a thresholded symbol flips category
%   (e.g. *** to **) on Monte Carlo noise alone when p sits near 0.001 or 0.05,
%   while the underlying effect is unchanged.

% Default resolution for callers that cannot supply their nBoot
if nargin < 2 || isempty(pRes) || ~isfinite(pRes) || pRes <= 0
    pRes = 1e-3;                            % conservative fallback
end

% Nothing to report
if isempty(p) || isnan(p)
    txt = '';
    op  = '';
    return
end

% Below the estimator's resolution: report the floor implied by pRes, rounded
% up to the next power of ten (e.g. pRes=0.00002 -> 0.0001), not a bare zero.
if p == 0
    txt = renderDecimal(10^ceil(log10(pRes)));
    txt = regexprep(txt, '0+$', '');        % a power-of-ten floor has one digit: '0.000100' -> '0.0001'
    txt = regexprep(txt, '\.$', '');       % drop a bare trailing dot
    op  = '<';
    return
end

op = '=';                                   % an actual value was resolved

txt = renderDecimal(p);

end


function s = renderDecimal(v)
% renderDecimal  Plain fixed-point text for a positive value < 1, at three
%   significant digits, trailing zeros kept. Never uses an exponent.

% Decimal places for 3 significant digits, never fewer than 3 places.
% e.g. v=0.0499 -> ceil(-log10(v)) = 2 -> 4 decimals -> '0.0499'
nDec = max(3, ceil(-log10(v)) + 2);

% Plain fixed-point rendering at that precision. Trailing zeros are kept
% (2026-11-13) so every value shows exactly three significant digits on a
% figure ('0.0160', not '0.016', next to '0.00942'). The p == 0 floor is
% shortened by its caller branch above.
s = sprintf('%.*f', nDec, v);

end
