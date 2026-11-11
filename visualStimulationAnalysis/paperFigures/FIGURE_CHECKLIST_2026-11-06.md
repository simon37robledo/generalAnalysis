# Which panels to check (Paper_figs root, state 6 Nov 2026 23:00)

Rule of thumb: grating panels are current only if their name has **Ex_40-97** (recordings 1 and 8 out).
The same panels with **Ex_1-97** for SDGm/SDGs pairs are from the night run (00:35-01:10) and are stale.
Times are file times on 6 Nov.

## Current — check these
| Figure | Panels | Time |
|---|---|---|
| 1 MB vs SB | `fig1_Zscore-Swarm/-Scatter`, `fig1_DeltaRate-Swarm/-Scatter`, `fig1_ResponsiveUnits` (all `Ex_1-97_Combined_MB-RG`; SpkRate panels removed 2026-11-10); `fig1_PSTH-MB-RG_allStims`, `fig1_PSTH-RG_offset`; rasters `fig1_Raster-MB-RG-all_union`, `fig1_Raster-RG-all_offset`; examples `fig1_PV140_..._2-*eNeuron-246*` (SB, phy 768) and `fig1_PV233_..._6-*eNeuron-1*` (MB, phy 357) | 10:27-10:43 |
| 2 spatial tuning | `fig2_SpatialTuningIndex-...allResp1-tp-5` (unpaired), `...allResp0-tp-5` (paired); `fig2_StripeFraction_vs_baseline`; `fig2_StripeNeurons*` (8); SB example `fig2_PV140_..._9-*eNeuron-15*` (phy 288) | 09:33-10:37 |
| 3 SDGm vs SDGs | `fig3_*-Ex_40-97_Combined_SDGm-SDGs` (6), `fig3_PSTH-SDGm-SDGs` | 18:10-18:11 |
| 4 locality | SB vs SDGs: `fig4_*-Ex_40-97_Combined_RG-SDGs` (6), `fig4_PSTH-RG-SDGs_allStims`. MB vs MBR: `fig4_*-Ex_1-97_Combined_MB-MBR` (6), `fig4_PSTH-MB-MBR_allStims`, `fig4_Raster-MB-MBR-all_union` (10:45, unaffected). MB vs SDGm: `fig4_*-Ex_40-97_Combined_MB-SDGm` (6), `fig4_PSTH-MB-SDGm_allStims`. MBR vs SDGm (PV97, descriptive): `fig4_*-Ex_49-54_Combined_MBR-SDGm` | 18:12-18:25 |
| 5 direction | `fig5_DirectionSelective_fraction`, `fig5_DSI_swarm`, `fig5_PreferredDirection`, `fig5_ExampleUnit_rec74_phy808_tuning` (18:35) and `_raster` (10:00) | 18:35 |
| 6 depth | `fig6_z_vs_depth`, `fig6_z_depth_rho_per_insertion`, `fig6_MBminusSB_z_vs_depth`, `fig6_responsive_vs_depth` | 18:31 |
| 7 Laudakia | `fig7_*-Ex_28-98_Combined_MB-RG`; `fig7_*-Ex_28-56_Combined_{SDGm-SDGs, MB-SDGm, RG-SDGs}` | 01:37-01:43 |
| Dataset | `summary_ResponsiveUnion-Ex_1-97_...` | 18:26 |
| Supplementary | `figS_NullCalibration` (redrawn 22:53 without NI); `figS_SpatialTuningIndex-...tp-5` to `tp-30` (threshold sweep; tp-5 duplicates Fig 2); `figS_LfpSurface_vs_coorZ` (3 Nov, unaffected) | |
| Numbers | `Stats summary tables\*.csv` (18:32-19:13) | |

## 2026-11-10: population PSTHs/rasters redrawn (split-half trials, onset-aligned)
- Main: `fig1_PSTH-MB-RG_allStims_cv_mbOnset_onOff` (two panels), `fig1_Raster-MB-RG-all_union_cv_mbOnset_onOff`
  (11 Nov: onset | offset axis, MB 0-2.3 s and SB 0-500 ms from onset, then 0-500 ms from each recording's offset);
  `fig3_PSTH-SDGm-SDGs_cv`, `fig3_Raster-SDGs-SDGm-all_union_cv` (two tiles, replaces the combined SDG raster);
  `fig4_PSTH-RG-SDGs_allStims_cv`, `fig4_Raster-RG-SDGs-all_union_cv`, `fig4_*-MB-MBR*_cv_mbOnset`, `fig4_*-MB-SDGm*_cv_mbOnset`.
- Supplementary (named by main figure, `sfigN_`): `sfig1_` crossing-aligned MB vs SB PSTH/raster (`_cv`),
  `sfig1_RG-OnOff-PSTH`, `sfig1_RG-OnOff-Raster`, `sfig1_PSTH-RG_offset_cv`, `sfig1_Raster-RG-all_offset_cv`; `sfig1_NullCalibration`;
  `sfig2_SpatialTuningIndex-*` (threshold sweep); `sfig4_` crossing-aligned MB-MBR / MB-SDGm; `sfig6_LfpSurface_vs_coorZ`.
- Scripts renamed: `SFig1_NullCalibration.m`, `SFig2_TuningIndexThreshold.m`; `FigS_LfpDepthCheck.m` keeps its name (prefix `sfig6`).
- Fraction-responsive panels (Figs 1, 3, 4): y-label "Fraction responsive", two-line footer with paper names (SB, not RG).
- MOVED to `Paper_figs\_superseded_2026-11-10` (32 files, list in its README): the 6 Nov `fig1_/fig3_/fig4_`
  PSTH/raster files without `_cv`, `fig1_*RG*_offset` (now `sfig1_`), and every `*_SpkRate-*` panel
  (SpkR = plain max inflated by the number of categories; replaced by `DeltaRate-Scatter`).
- Fig 7 (SA) rerun 2026-11-10 21:06-21:08: `fig7_DeltaRate-Scatter-*` added (4 comparisons), new fraction-panel
  footer; SA leave-one-animal-out tables byte-identical to before.

## Stale — MOVED to `Paper_figs\_superseded_2026-11-06` (65 files + `Sensitivity_noPR`; list in its README)
- `fig3_*`, `fig4_*` and structs with `Ex_1-97_Combined_SDGm-SDGs`, `..._RG-SDGs`, `..._MB-SDGm` (00:35-01:10): include recordings 1 and 8.
- Grating rasters `fig3_Raster-SDG-all_union`, `fig4_Raster-RG-SDGs-all_union`, `fig4_Raster-MB-SDGm-all_union` (10:44-10:46): include recordings 1 and 8.
- Old Fig 2 MB example `fig2_PV132_..._7-*eNeuron-23*` (phy 180): replaced by rec 74 phy 727 (not drawn yet).
- Everything dated before 5 Nov in the root (40 files: old Fig 7 runs `Ex_44-98`/`Ex_56-56`, old example RFs eNeuron-111/-19/-14, `fig5_ExampleUnit_rec75_phy32_*`, `fig6_z_vs_depth_MB_SB`, Aug `SupplementaryStatsTable`), except `figS_LfpSurface_vs_coorZ.jpg`.
- Folder `Sensitivity_noPR`: noPR is now the main analysis, so this duplicates the root (5 Nov run, before the grating fixes).
- Folder `Archive_standardBC_20261105` = the standard-criterion sensitivity. Only its Fig 1 (MB vs SB) is a clean comparison; its grating panels also have the old grating onsets and recordings 1 and 8.

## Drawn 6 Nov 23:52-23:55 (MATLAB desktop)
- Fig 2 MB example (7 Nov, rule: tuning index >= median, no stripe, interior RF cell, highest MB z) = rec 74 phy 292:
  `fig2_PV140_Experiment_4_6_25_2-Speed1-MovBall-SelectedTrials-eNeuron-8`, `...-RF-sep-Dir-3.1416-...-eNeuron-8`.
  (Earlier picks rec 74 phy 727 and rec 75 phy 66 moved to `_superseded_2026-11-06`.)
- `fig3_Raster-SDG-all_union` (98 units), `fig4_Raster-RG-SDGs-all_union`, `fig4_Raster-MB-SDGm-all_union` (recordings 1 and 8 out).
  Note: the dashed motion-onset line is at 1 s; recordings 40-43 start motion at 1.5 s.
