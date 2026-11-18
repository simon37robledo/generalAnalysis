# Eye tracking: what we need back

Thanks for running your DLC-analysis script on these recordings. This page only describes the
**output** we need, so the MATLAB side can read your files without changes to your code.

## Input
`W:\Large_scale_mapping_NP\EyeTracking\manifest.csv` (same storage as `\\sil3\data\...`).
One row per insertion. The columns you need:

| column | meaning |
|---|---|
| `experimentNum`, `animal`, `insertion` | key of the insertion (do not renumber) |
| `dlc_csv` | the DLC output to analyse (raw DLC CSV, not filtered/inverted) |
| `dlc_video` | the video DLC ran on; coordinates are in its pixels (`frame_w` x `frame_h`) |
| `n_dlc_frames` | number of frames in the DLC CSV |
| `had_inverted_file` | an older flipped copy exists for this insertion (mirror / upside-down camera) |
| `status_elena` | only rows marked `ready` can be analysed; the others show why they are missing (e.g. `needs_dlc`) |
| `note` | e.g. the DLC was run on a cropped/enhanced copy of the video |

(`status_sync` is for our side only: you can ignore it.)
| `out_table` | file name we expect back for this row |

## Output: one CSV per insertion
Write to `W:\Large_scale_mapping_NP\EyeTracking\tables\<out_table>`
(e.g. `eye_49_PV97_ins1.csv`).

**Required**
1. **One row per video frame, for every frame of the DLC CSV**, in frame order. Frames where the
   pupil fit fails stay in the table with NaN values. Please do not drop them: we map row *i* to the
   *i*-th camera trigger, so a dropped row shifts every later frame in time.
2. A frame-index column, and say whether it starts at 0 or 1.
3. Pupil centre x and y **in the pixels of `dlc_video`**, with no flips. If your script flips or
   crops, tell us which flip or crop it applied.
4. Pupil size: ellipse axes or diameter, in pixels.
5. Any quality measure your script gives (DLC likelihood, number of points used, valid flag).

**Please also tell us** (in a short text file next to the tables or in an email):
- the likelihood threshold and any filtering/smoothing applied,
- the column names and units, if they differ from the above.

**Optional, if easy**: for one frame per animal, the eyelid-retractor length (3 mm) and the visible
eye-globe width, both in pixels. We use them to convert pupil shifts to degrees. If not, we do it.

## Start with one
Please send one finished insertion first (e.g. `experimentNum` 49, PV97 ins 1). We will check it
against our frame count and timing before you run the rest.
