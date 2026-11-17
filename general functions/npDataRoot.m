function root = npDataRoot()
% npDataRoot  Root of every analysis output: W:\Large_scale_mapping_NP (W: = \\sil3\data).
%
%   root = npDataRoot()
%
%   Paper figures (Paper_figs), pooled caches (lizards\Combined_lizard_analysis) and plots
%   are written under this root for every recording, wherever its raw data lives.
%   Before 2026-11-18 the root was derived from a recording's own path, so recordings on
%   Y: (\\sil1, e.g. SA5) wrote figures to Y:\...\Paper_figs and, lacking a 'lizards'
%   folder, pooled caches to a relative path inside the recording folder.
root = 'W:\Large_scale_mapping_NP';                                   % sil3
assert(isfolder(root), 'npDataRoot:unreachable', ...
    '%s (sil3) is not reachable; outputs would be lost or misplaced.', root);   % fail loudly
end
