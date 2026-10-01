# Paper 2 — Late S-wave energy and decorrelation in Utah FORGE downhole microearthquakes: elastic tests of smooth, layered, fracture-zone and stochastic media

Follow-up to Paper 1 (`../abg_framework_2026`). Tests the GFI concept (Leary, *Matlab Tools for ID Geothermal
System Flow Analysis*, Web Session 4): crustal properties are spatially correlated **pink-noise** fields
(well-log S(k) ~ 1/k; lognormal k ~ exp(αφ)). The pink-noise crust is compared with homogeneous, layered and
anisotropic (VTI) crusts, using 3D elastic finite-difference synthetics and the Utah FORGE April 2022 (well 16A
stages 1–3) microearthquake waveforms recorded in monitoring well 56-32.

## What is (and is not) in this repository
- Included: all MATLAB code (`code/`), figures (`figs/`), manuscript builder and outputs (`manuscript/`), and every
  number quoted in the manuscript as CSV (`results/*.csv`).
- Not included (size or data-release reasons, see `.gitignore`): strain Green's tensors `sgt/`, `sgt_hf/`
  (~75 MB per model; regenerate with `run_sgt(<model>)` / `run_sgt_hf(<model>)`), the extracted FORGE 2022 56-32
  waveforms `data/` (regenerate from the GES SEG-2 files with `s01_extract_forge_events.m`), intermediate
  `results/*.mat` and run logs.
- Requirements: MATLAB R2025b/R2026b (Signal Processing, Statistics toolboxes); Python 3.12 with `python-docx`, `pandas`,
  `lxml`, `latex2mathml`; Microsoft Word (for `MML2OMML.XSL`, used by `manuscript/springer_doc.py`, and for PDF export).
- The FD solver builds on Igel, Mora & Riollet (1995) and H. Igel's acoustic FD MATLAB codes (Igel 2016).

## Revision after internal review (2026-10-01) — CURRENT STATE
The sections further below document the original analysis and are partly superseded (window overlap, 150-event
subset, σ interpretation, "pink crust" framing).

**Documents** (`manuscript/`, Python 3.12): `build_paper2_rev.py` → `Paper2_revised.docx` (EMBED=1 →
`Paper2_revised_review.docx/.pdf`) and `build_additional_file1.py` → `Paper2_Additional_file_1.docx/.pdf`. All numbers are read from `results/*.csv`.

**New simulations** (`queue_revision.sh`, horizontal forces only): smooth 3-D gradients `GRAD3`, 6-m sonic-log
layering `L6`, 12 planar fracture zones `FZ`, variable Vp/Vs `VSD`/`VSO`, and a common family grid
`X_<pink|exp|gau>_<σ>_<p|a>_<seed>` (σ 0.09/0.13/0.18 × 2 shapes × seeds 11/33); convergence cases
`run_convergence('C5'|'C25'|'CBIG'|'CSP')` (5 m, 2.5 m, domain +150 m, sponge 40). Post-processing: `post_revision.m`.

| script | purpose | output |
|---|---|---|
| `s19_coda_prediction.m`, `mt_fit_variants.m`, `s27_summarise_s19.m` | per-event late-energy test, non-overlapping windows (fit P−10…S+15 ms, coda S+40…100 ms), all 416 two-sensor events, full/dev/DC MT, window grid | `table_coda_prediction.csv`, `table_coda_window_sensitivity.csv`, `table_mt_variants.csv` |
| `s20_family_inference.m` | Φ (explicit weights) with event bootstrap, family best member, holdout | `table_phi_bootstrap.csv`, `table_family_bootstrap.csv`, `table_holdout*.csv` |
| `s21_ringing_test.m` | in-band spectral smoothness, null-ringing test | `table_ringing_prominence.csv`, `table_null_ringing.csv` |
| `s22_calib_bootstrap.m` | sensor calibration bootstrap + leave-one-stage-out | `table_calib_*.csv` |
| `s23_source_spectra.m` | apparent corners, Brune source-duration sensitivity | `table_source_*.csv` |
| `s24_stage_and_fields.m` | stage comparison, field statistics after truncation | `table_stage_coda.csv`, `table_field_stats.csv` |
| `s25_d2_uncertainty.m` | D2 cluster bootstrap + location perturbation | `table_d2_uncertainty.csv` |
| `s26_convergence.m` | coda convergence (grid, domain, sponge) | `table_convergence.csv` |
| `s28_s_timing_check.m` | S-timing check, well-timed subsets | `table_s_timing_subsets.csv` |
| `s29_sigma_eff.m` | σ_eff per family/sensor | `table_sigma_eff.csv`, `table_sigma_curves.csv` |
| `s30_selection_table.m` | event-selection flow (800 → 555 → 416 → 150) | `table_event_selection.csv` |
| `s18_log_spectrum_uncertainty.py` | sonic-log slope variants, spectral-model AIC | `table_log_*.csv` |

**Key revised results** (median log10 predicted/observed late energy, sensor A/B, 416 events):
smooth H -1.80/-1.09, L -1.64/-1.17, GRAD3 -1.67/-1.11, VTI -1.91/-1.19; deterministic small-scale L6 -1.28/-0.53, FZ -1.04/-0.39;
pink σ .045 -0.74/+0.02, σ .09 -0.14/+0.53, σ .13 +0.31/+0.78. Smooth media under-predict by ×12–82; L6/FZ supply much of
the energy but keep the coda of neighbouring events coherent (D2), unlike the data; stochastic media reproduce both.
σ_eff (reference domain) 0.05–0.13 by family/sensor; the domain test (+0.40 log units for a +150 m domain) makes these
upper bounds (÷ ~1.7, i.e. ~0.03–0.08 vs sonic log 0.045). Grid (2.5 m) and sponge changes are small (+0.01, +0.13).
The 20–80 Hz data do not select a spectral family (seed scatter in Φ up to 0.33 > family differences).
Corrected errors: window overlap; MT residual/overlay shift sign (VR and coda ratios unaffected); FD interior size.

## Final revision (second internal review, 2026-10-01) — release `v1.0-submission`
- D2 (inter-event coda coherence) recomputed in the pure late-coda windows S+40–100 and S+50–110 ms
  (`s11_shape_diagnostics(name,'forge',win)`, `run_d2_windows.m`, `s31_d2_windows.m`), with a noise-only reference
  window (S+550–610 ms), all-pairs medians and a location stress test of 5–50 m applied to data and media alike
  (`fig_d2_perturb.m`). Result: in the late coda all media approach the noise-limited coherence; smooth media and planar
  fracture zones stay more coherent than the data (19 of 20 cases), the 6-m layering and stochastic media come closest.
  The claim is therefore "structure at 5–10 m is required; decorrelation favours distributed 3-D heterogeneity without
  excluding fine layering".
- Main figures redrawn at print size (17.4 cm, 7–9 pt, 300 dpi PNG + vector PDF) via `pub_export.m`;
  new geometry figure `fig_geometry_rev.m`.
- Authors: Wiwit Suryanto (UGM) and Peter Leary (Geoflow Imaging).

## Pipeline (MATLAB R2025b/R2026b, run from `code/`)
| step | script | output |
|---|---|---|
| FD solver | `fd3d_elastic_vti.m` — velocity–stress staggered grid O(2,4), VTI (Thomsen), moment-tensor or point-force source, strain-rate recording | — |
| validation | `test_fd_homog.m` (P moveout 5833 vs 5800 m/s), `test_reciprocity.m` (SGT vs direct, corr 0.995) | `figs/val_reciprocity.png` |
| pink noise | `pinknoise3d.m` (GFI recipe, slide 7), `build_model.m`, `model_catalog.m` | — |
| S01 | `s01_extract_forge_events.m` — match GES 2022 catalogue to 56-32 SEG-2 (GPS time), 800 events × 6 ch × 1 s | `data/forge2022_events.mat` |
| S02/S03 | timing QC; AIC picks; **sensor self-calibration** (position, clock, Vp) | `data/forge_picks_sensors.mat` |
| S04 | 56-32 ThruBit sonic porosity → σ(lnV), S(k) slope, 1D layering | `data/log5632_stats.mat`, `figs/fig_log5632.png` |
| S05 | travel-time inversion H / gradient / VTI with AIC/BIC; residual spatial correlation | `data/tt_inversion.mat` |
| S06 | data waveform features (horizontal energy, Ricker-convolved), noise library | `data/forge_features.mat` |
| SGT runs | `run_sgt.m <model>` — 2 sensors × 3 forces, strain rates at 555 event nodes | `sgt/<model>.mat` |
| S07 | synthetic features (FORGE-like DC mechanisms, real noise at matched SNR) | `results/feat_<model>_<mech>.mat` |
| S08 | deterministic waveform inversion: per-event moment tensor, VR, unexplained coda | `results/mtinv_<obs>.mat` |
| S09 | statistical waveform inversion: KS + envelope + TT-statistics misfit; model ranking | `results/compare_<obs>_<mech>.mat` |

## Models (`model_catalog.m`)
H homogeneous (Vp 5808, Vp/Vs 1.73) · L layered (56-32 log, 50-m blocks) · A VTI (TT-inverted ε, δ) ·
P1 pink σ=0.045 p=1.5 · P2 σ=0.045 p=1.8 · P3 σ=0.09 p=1.5 · P4 σ=0.02 p=1.5 · P5 σ=0.045 p=1.2 ·
T1 synthetic truth = pink σ=0.045 p=1.5 (other realisation) · T2 synthetic truth = layered + VTI (ε .06, δ .03, γ .06).
Grid 113×137×149 at 5 m (sponge 20 cells), dt 0.25 ms, 1300 steps, Ricker 45 Hz; analysis band 20–80 Hz.

## Key results (2026-09-27)
| test | best | deterministic H / L / A | notes |
|---|---|---|---|
| synthetic truth T1 (pink σ .045) | R1 Φ=0.09 (same stats, other seed) | 0.29 / 0.29 / 0.29 | σ ≥ 0.09 rejected (Φ ≥ 0.51); p weakly resolved |
| synthetic truth T2 (layered+VTI) | L 0.13, A 0.16 | H 0.22 | all pink rejected (Φ ≥ 0.27) → no pink bias |
| FORGE data | P9 σ .18 0.59, P6 σ .13 0.62 | 1.42 / 1.32 / 1.30 | plateau σ≈0.13–0.18; σ=.045 realisations 1.09–1.27; unchanged with random mechanisms |

Deterministic MT waveform inversion (VR) always prefers the smoothest model (T1: H 0.93 vs P1 0.88), i.e. it cannot identify
a pink-noise crust; the statistical (coda/envelope) inversion can. Tables: `results/table_*.csv`; figures: `figs/fig_*.png`.

## Manuscript
`manuscript/build_paper2.py` (Python 3.12: `py -3.12 build_paper2.py`, `EMBED=1` for review copy) reads `results/*.csv`
→ `Paper2_manuscript.docx`, `Paper2_manuscript_review.docx/.pdf`. Placeholders: co-authors, funding, data URL.
Rebuild order: s01 → s03 → s04 → s05 → s06 → s06b → run_sgt(each model) → s07 → s09 → s08 → s10 → fig_* → build_paper2.py.

## Update 2026-09-28 — larger σ and competing random media
- New models: P12 (σ .22), P10 (σ .25; dt sub-stepped ×2 by `run_sgt`), P11 (σ .18, p 1.8), R5 (σ .18 realisation);
  competitors E15/E50 (exponential ACF) and G15/G50 (Gaussian ACF), σ = 0.13, `randmedium3d.m`.
- FORGE Φ: pink optimum σ≈0.18 (0.59; bracketed by 0.13→0.62 and 0.22→0.72, 0.25→0.79). Competitors at σ .13: G50 0.47,
  E15/E50 0.55, G15 0.61 → **coda constrains heterogeneity strength, not spectral shape**.
- Per-event observed-vs-predicted coda (`fig_seismograms.m`, `fig_coda_ratio.m`, `results/table_coda_ratio.csv`):
  smooth crusts under-predict coda ×12–30 (1–11 % events within ×3); pink σ .09 unbiased; G50 72/64 % within ×3.
- Shape from the 56-32 log (`fig_log_vs_media.m`): rms vs log spectrum 0.29 pink, 0.50–0.56 exponential, ≥2.5 Gaussian.
- Manuscript retitled and rewritten accordingly (`build_paper2_v1_backup.py` = previous version).

## Update 2026-09-28 (evening) — equal-σ shape-resolution test and shape diagnostics
- Truths at σ = 0.13, seed 33: R6 (pink), E15b (exponential a=15), G50b (Gaussian a=50); candidates seed 11.
- Φ recovers pink (P6 0.25 vs 0.33) and Gaussian (G50 0.36 vs 0.38) truths but mistakes exponential for pink
  (P6 0.21 vs E15 0.24). D1/D2 (`s11_shape_diagnostics.m`) carry no shape information (combined score always → G50).
- D2 inter-event coda coherence is new independent evidence against smooth crusts: data 0.69→0.49 (0–5→80–160 m),
  H 0.86→0.84, VTI 0.88→0.79; all scattering media decorrelate like the data.
- FORGE Φ leans to G50 (0.47 vs P6 0.62), opposite to the log. Manuscript states: seismic data reject smooth crust and
  require strong scattering; pink spectral form rests on the borehole data. Tables: `results/table_shape_test.csv`,
  `table_coherence.csv`; figure `figs/fig_shape_data.png`.

## Update 2026-09-29 — repeated truth realisations (`s13_confusion.m`)
- Truths: pink / exponential a=15 / Gaussian a=50, seeds 33, 44, 55 (R6, R6c, R6d, E15b–d, G50b–d); candidates seeds 11 & 33.
- Confusion: 5/9 correct; small-scale-rich vs Gaussian 8/9; pink vs exponential 3/6 (chance). Cause: pink and
  exponential a=15 fields from the same noise are 96 % correlated over 10–700 m (`results/table_field_corr.csv`).
- FORGE family means (2 realisations): pink 0.61, exponential 0.54, Gaussian 0.51 (G50 seed11 0.47 vs seed33 0.56) →
  spread 0.09 < median truth-test spread 0.14 → data prefer no family. Earlier "leans Gaussian" was a realisation effect.
- Tables: `table_confusion_runs.csv`, `table_confusion.csv`, `table_data_families.csv`. Manuscript Table 3 replaced.

## Update 2026-09-29 — bandwidth, multi-octave coda, high-frequency runs
- `fig_event_spectrograms.m` → figs/fig_event_spectrograms.png, fig_bandwidth.png, results/table_bandwidth.csv:
  SNR>3 to ~1.8 kHz for all events; spectral peaks ~330/480/920 Hz = suspected tool/coupling resonances.
- `s14_coda_multioctave.m`, `s14b_clean_band_test.m`: 7 octave bands 25–1300 Hz; clean band 25–200 Hz:
  Qc^-1 ∝ f^-1.48 (A) / f^-1.39 (B), single power law (break disfavoured). Tables: table_coda_multioctave.csv,
  table_coda_cleanband.csv; fig_coda_multioctave.png.
- `s15_model_coda_multioctave.m` + `coda_band_measures.m`, `s16_qc_slope_models.m`: Qc^-1 slope 25→50 Hz, data −1.46
  [−1.58, −1.26]; pink −1.29…−1.68, exp15 −1.35…−1.81, Gauss50 −0.59…−1.44 (1/4 inside), smooth −0.83…−0.87.
  table_qc_slope.csv, fig_qc_slope.png.
- High-frequency FD (`forge_setup_hf.m`, `run_sgt_hf.m`, `synth_records_hf.m`): dx 2.5 m, Ricker 90 Hz, sensor B,
  models H P6 E15 G50 R6 E15b G50b → sgt_hf/ (running). Manuscript: new Results subsections + Figs 13–16, Table 4.

## Update 2026-09-30 — high-frequency runs and correlation-length break test (`s17_break_test.m`)
- 2.5-m grid (sensor B, 50/100 Hz) for H, P6, R6, E15, E15b, G50, G50b and exponential a = 5, 9, 17, 34, 67, 134 m
  (`queue_aseries.sh`); a = 9–134 m also on the 5-m grid. Table `results/table_break_test.csv`, `figs/fig_break_test.png`.
- Qc^-1 curvature (slope 50–100 minus 25–50 Hz) is NOT diagnostic: models −0.08…+0.49, no trend with a; data −0.04.
- 100-Hz coda level (log10 coda/S, sensor B): data −1.58; pink −1.20/−1.50; exponential −1.03 (a 9 m) … −1.82 (a 134 m);
  Gaussian a 50 m −2.20/−3.01; homogeneous −3.17 → large-scale Gaussian and smooth crusts rejected at 100 Hz.
- 50–100 Hz Qc^-1 slope: data −1.50 [−1.70, −1.18]; only exponential a ≤ 9 m and pink (edge) inside; pink is the only
  medium close to both criteria (single realisations, σ fixed → suggestive).
