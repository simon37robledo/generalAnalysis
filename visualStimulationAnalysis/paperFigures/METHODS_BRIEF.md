# Methods brief — analysis and statistics (state of 2026-10-20)

Source material for rewriting the Methods. Facts only; prose and style come from the
existing Methods text. Every number below is produced by the code in
`visualStimulationAnalysis/` (figure scripts in `paperFigures/`, shared design in
`figureDesign.m`) and the tables in `W:\Large_scale_mapping_NP\Paper_figs\` and
`...\Paper_figs\Stats summary tables\`. Items marked **[fill]** need a value from the lab
records; items marked **[verify]** were read from one or two sessions only.

## How to use this brief (new session)
Paste into the new session: (1) your current Methods text, (2) this file, (3) CLAUDE.md.
Ask it to rewrite the Methods in your style, adding the stimuli below (NI pending, NV out).

---

## 1. Data set and unit inclusion
- Species/recordings: *Pogona vitticeps* (PV), 8 animals, 44 recordings (one insertion
  each), Neuropixels 2.0 (some 4-shank) in aDVR. *Laudakia vulgaris* (SA): Figure 7, pending.
- Spike sorting: Kilosort4; unit classification BombCell; only "good" somatic units.
  877 somatic units recorded in the analysed recordings.
- **Per-block unit presence (new).** Units drift in and out within a session. A unit counts
  for a stimulus only if present throughout that stimulus block: the block (first onset to
  last onset + stimulus duration) is split into 10 equal bins; presence ratio = fraction of
  bins with >= 1 spike; present if >= 0.7 (BombCell's presence ratio, computed per block
  rather than over the whole recording, because the issue is drift between blocks). In
  every comparison a unit must be present in ALL compared blocks. Recordings 87, 92, 94
  have no unit present in any compared block and drop out of the per-insertion statistics.
  Units recorded vs analysable are reported per figure (e.g. Figure 1: 877 recorded,
  730 present in both MB and SB blocks).

## 2. Stimuli and response windows
All windows are fixed before analysis and identical across figures (one statistics cache
per stimulus per recording). Baselines are blank pre-trial periods, capped by the
inter-trial interval (ITI - 100 ms). Luminosity locked to 255 (white) for MB and SB.

| Stimulus | Protocol (per trial) | Response window(s) | Baseline |
|---|---|---|---|
| Moving ball (MB) | ball 120 px, 500 px/s [verify, px -> deg **[fill]**], 4 directions [verify older sessions 1, 8], parallel offsets, sweep ~2.3 s, ITI 1 s | 300 ms from the ball crossing each of 9 x 9 screen cells (spatial-grid test) | 300 ms pre-trial |
| Static ball / rect grid (SB, RG) | 81 positions (9 x 9), 500 or 1000 ms (two protocols), ITI 1 s | 300 ms ON from onset + 300 ms OFF from each recording's measured offset | 300 ms pre-trial |
| Moving bar (MBR) | bar sweeping the screen, 4 directions, ~2.3 s | 9 equal consecutive segments of the sweep (9 = columns/rows of the grid) | same length as one segment |
| Static grating (SDGs) | first phase of the SDG trial: static grating 1000 ms (1500 ms in recs 40-43) | first 1000 ms of the static phase (longest window inside every static phase) | same length, capped by ITI (900 ms; 400 ms in recs 40-43) |
| Drifting grating (SDGm) | 4 angles (0/90/180/270; 45-degree steps in rec 1), TF 0.52 and 2.08 Hz, SF 0.0042/0.0083 cyc/px [verify, **[fill]** cyc/deg], drift 2-5 s | first 1923 ms after motion onset (= one period of the slowest TF, 1/0.52 Hz; fixed a priori) | blank pre-trial (before the static phase), same length capped by ITI |
| Full-field flash (FFF) | 500/1000 ms, ITI 2 s | 300 ms from onset | 300 ms pre-trial |
| Natural images (NI) | **pending** (lizard posture images; trials [1:15 61:75]; stimuli in a subfolder) | 300 ms from onset | 300 ms pre-trial |
| Natural videos (NV) | **excluded** from the paper (too few present units, late in sessions) | — | — |

Why the grating windows differ from 300 ms: with 300 ms, drifting-grating responses sat
at chance (10 of 258 units); with one temporal period they did not (56/258 at ~2 s), while
static gratings stayed at chance at every window (diagnostic run 2026-10-18; null FPR at
the long windows 4.3-5.9%).

Display conventions (figures only): PSTHs drawn as steps between bin edges; 50 ms bins
for windows <= 1 s, 100 ms for 2 s windows; SDGs drawn only over its 1000 ms test window;
MB population PSTH aligned to the RF crossing (t = 0 = ball enters the preferred cell).

## 3. Stimulus timing quality control
- Onsets from the photodiode (sync square). For every recording and stimulus the median
  lag between the diode onset and the trial TTL was checked: RG ~46 ms, FFF ~48, MBR ~54,
  NI ~50, MB 60-83 ms (consistent within each stimulus).
- Two recordings (68 RG, 84 FFF) had onsets 100.2 ms late in every trial: the image-trial
  sync requires an exact flip count, a few spurious analog crossings triggered a fallback
  extraction that returned a constant 100 ms offset. Fixed: after any fallback, every
  onset/offset is compared with the raw analog diode edges and snapped to them if the
  offset is consistent (error otherwise). Both recordings resynced; figure numbers
  unchanged (recording 68 has no unit present in both MB and SB blocks).
- SDG: the static-phase onset is derived as the first drifting flip minus the static
  duration (the static frame and its diode flip are verified to coincide in recordings
  52, 84, 97; recordings 40-43 about 35 ms late, minor).

## 4. Responsiveness test (all stimuli)
- Per unit and trial: mean rate in the response window minus that trial's baseline mean.
- Statistic: maximum over conditions (and over cells x directions for MB, positions x
  {ON, OFF} for SB, conditions x segments for MBR) of the condition means, with the
  null-mean subtraction for the max (winner's curse). Null: 10,000 sign flips of the
  per-trial differences, the same maximum taken AFTER each flip (Nichols & Holmes 2002).
  Right-tailed p = fraction of null maxima >= observed. Responsive: p < 0.05, no FDR
  (the max-statistic controls the family-wise error within a stimulus).
- A design that took each trial's best window before the sign flip (MB sliding window)
  gave ~100% false positives and was discarded.

## 5. Null calibration with artificial neurons (supplementary figure)
Simulated units with no response (homogeneous Bernoulli spiking, 1 ms bins, at the same
rate before and during the stimulus) on the trial layout of a real recording, through each
test. Per test design: 240 simulated units at each of 6 firing rates (0.5, 1, 2, 5, 10,
20 Hz) = 1,440 units, 2,000 sign flips each; 9 designs (7 used, 2 discarded) = 12,960
simulated units in total. Layouts: MB and SB from recording 97 (MB 36 conditions x 10
trials, 2367 ms sweep; SB 81 positions x 10 trials); FFF, SDGs and SDGm from recording 52;
MBR from recording 70; NI from recording 92. Pooled false-positive rates at alpha = 0.05
(false positives / 1,440): MB grid 4.9%, SB ON+OFF 3.5%, MBR 9
segments 5.0%, SDGs 5.3%, SDGm 4.1%, FFF 3.1%, NI 3.7%; discarded MB sliding window 99.9%.
Clopper-Pearson CIs per rate in `Null_calibration_*.csv`; figure `figS_NullCalibration.jpg`.
(The FP/FN raster screen was a display-concordance QC and is not reported.)

## 6. Response measures
- z-score (ZScoreU): response minus baseline at the preferred condition, divided by the SD
  of the per-trial baseline means (across trials).
- Delta firing rate (spikes/s), split-half cross-validated: the preferred condition (and
  window/cell) is selected on odd repetitions and measured on even ones, and vice versa;
  averaged. Unbiased under H0 (`Null_calibration_delta_rate_bias.csv`). MB with two
  speeds (recordings 49-54): the speed with the lowest p is used, a mild selection effect.
- Fraction responsive: per insertion, over the units present in the compared blocks.

## 7. Cross-animal inference
- Hierarchical bootstrap (Saravanan et al. 2020, `hierBootMatchFreq`): animals, then
  insertions within animals, then units within insertions, with replacement, matching the
  observed counts; 10,000 replicates; percentile 95% CIs; two-tailed p (or one-tailed where
  stated). Paired comparisons use per-unit differences matched by unit identity.
- Leave-one-animal-out for every comparison: `Leave_one_animal_out_summary.csv` (p range,
  least robust animal, runs still at p < 0.05).
- Benjamini-Hochberg across the depth tests of Figure 6.
- Insertions with no analysable unit are excluded from fraction comparisons (a 0/0
  fraction had deflated earlier p-values; fixed 2026-10-19).

## 8. Figure-specific analyses
**Figure 2 — spatial tuning and stripe RFs.** RFs from spike-triggered stimulus maps
(MB: convolution of the spike train with the ball trajectories; SB: rate per position),
9 x 9 grid. Spatial tuning index L_amplitude_diff (top 5% of cells vs the rest, shuffle
subtracted). Split half: MB preferred direction chosen on odd repetitions, index measured
on even (and vice versa, averaged); SB averaged over the same halves. Stripe RFs: Radon
anisotropy on the 9 x 9 map vs a pixel-permutation null; MB takes the maximum over
directions in the observed score and in every surrogate; accepted orientations 20-55 deg.
Pixel permutation is calibrated for SB (non-responsive units 1.4% stripes) but not fully
for MB (convolution smoothing; non-responsive units 5.4% stripes), so stripe fractions are
reported against these empirical baselines, not as p-values. Comparisons that select and
measure on the same RFs are descriptive.

**Figure 5 — direction tuning.** Baseline-subtracted rate per direction in the same windows
(MB at the RF-cell crossing, cell chosen on the direction-averaged response; SDGm 1923 ms;
MBR whole sweep). DSI/OSI by vector sum (Ringach et al. 2002; Mazurek et al. 2014);
per-unit permutation test shuffling direction labels; DSI shown minus its shuffled mean
(raw DSI of an untuned unit ~0.65 with these trial counts). Grating angle a -> drift
direction pi/2 - a in the ball convention (from the stimulus code: procedural sine
grating, phase increasing, rotated clockwise). Same tuning across stimuli: correlation of
tuning curves (units responsive to both) and preferred-direction agreement (units
selective to both), against within-insertion re-pairing of units.

**Figure 6 — depth.** Unit depth below the brain surface from Kilosort's channel positions
(cluster_info.tsv + channel_positions.npy) and coor_Z (to be replaced by the LFP-based
surface estimate of [collaborator, **[fill]**]; validation against coor_Z and DiI tracks).
Test: within-insertion rank correlation (depth and measure ranked inside each insertion);
p from permuting depths within insertions (calibrated: 6.3% FPR with between-insertion
confounds); CI from the hierarchical bootstrap of the same statistic; BH across 22 tests.

## 9. Key results (current numbers)
- **Dataset:** 354 of 774 testable units (45.7%) respond to >= 1 stimulus
  (MB 253, SB 132, SDGm 59, MBR 31, SDGs 7, FFF 3).
- **Fig 1 (MB vs SB, 297 units, 8 animals):** z 1.40 [0.20, 2.15] p = 0.015;
  delta rate 1.43 [0.08, 2.89] spikes/s p = 0.037; fraction responsive 0.141
  [0.015, 0.283] p = 0.033. LOAO: z holds in 6/8 runs, delta rate and fraction in 4/8
  (largest p without PV140, 0.12).
- **Fig 2:** spatial tuning index MB vs SB not different (unpaired one-tailed P = 0.27;
  LOAO 0.18-0.52); paired MB - SB in 76 shared units 0.89 [0.03, 3.53] p = 0.041, not
  robust (4/5 LOAO runs n.s.). Stripe RFs: SB 9/126 (7.1%; baseline 1.4%), MB 26/232
  (11.2%; baseline 5.4%; 15 of the 26 from recording 81).
- **Fig 3 (SDGm vs SDGs, 64 units, 6 animals):** z 0.84 [0.22, 1.22] p < 0.001; delta rate
  1.36 [0.17, 2.31] p = 0.007; fraction 0.23 [0.07, 0.38] p = 0.0002; robust in all LOAO runs.
- **Fig 4:** MB vs MBR (121 units, 3 animals): z 2.26, delta rate 2.96, fraction 0.26, all
  p < 0.001. MB vs SDGm (104 units, 6 animals): z 1.00 [0.55, 1.41] p < 0.001, delta rate
  1.74 [0.36, 3.47] p = 0.008, fraction n.s. (p = 0.63). SB vs SDGs (26 units, 6 animals):
  z 1.70 p = 0.037, delta rate 1.91 p = 0.036, not robust (3/6 LOAO runs). MBR vs SDGm:
  descriptive (PV97 only).
- **Fig 5:** MB direction-selective 111/253, per-insertion fraction 0.35 [0.13, 0.52],
  p = 0.003 vs 5% (robust); SDGm 7/59, MBR 3/31 not above chance. Same tuning:
  MB-SDGm r = -0.13 (15 units) p = 0.34; MB-MBR r = 0.23 (20 units, 3 animals) p = 0.044.
- **Fig 6 (coor_Z depths):** no robust depth relation (22 tests, BH); MB stripe yes/no
  q = 0.050 driven by PV140 / recording 81.

## 10. References
Nichols & Holmes 2002 (max-statistic permutation); Saravanan, Berman & Sober 2020
(hierarchical bootstrap); Benjamini & Hochberg 1995; Ringach, Shapley & Hawken 2002 and
Mazurek, Kager & Van Hooser 2014 (circular variance OSI/DSI); Clopper & Pearson 1934;
Pachitariu et al. 2024 (Kilosort4); Fabre et al. 2023 (BombCell); Siegle et al. 2021
(presence ratio).
