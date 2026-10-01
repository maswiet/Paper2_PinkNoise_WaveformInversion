"""Supporting Information of the GJI version of Paper 2 (Figures S1-S11, Tables S1-S18).
Run:  py -3.12 build_supporting_gji.py   -> Paper2_GJI_Supporting_Information.docx (figures embedded)
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from springer_doc import SDoc  # noqa: E402

R = P / "results"; F = P / "figs"
PENDING = []


def rd(f):
    fn = R / f
    if not fn.exists():
        PENDING.append(f); return None
    return pd.read_csv(fn)


d = SDoc(line_numbers=False, spacing=1.15, size=11)
FIGN = [0]; TABN = [0]


def fig(path, caption):
    FIGN[0] += 1
    if Path(path).exists():
        d.picture(path, 16.5)
    else:
        PENDING.append(Path(path).name)
    p = d.d.add_paragraph(); r = p.add_run(f"Figure S{FIGN[0]}. "); r.bold = True
    d.runs(p, caption)


def tab(title, header, rows, widths, legend=None, size=8):
    TABN[0] += 1
    d.table(title, header, rows, widths, legend=legend, size=size, label=f"Table S{TABN[0]}.")


def pending_tab(title, what):
    TABN[0] += 1
    d.p(f"**Table S{TABN[0]} {title}** [PENDING: {what}]")


f2 = lambda x: "—" if pd.isna(x) else f"{x:+.2f}"
TITLE = ("Supporting Information for \"Late S-wave coda of Utah FORGE downhole microearthquakes: elastic tests of "
         "smooth, deterministic and stochastic media\"")
d.p(f"**{TITLE}**", align="center", size=13)
d.p("Wiwit Suryanto (Universitas Gadjah Mada) and Peter Leary (Geoflow Imaging)", align="center")
d.p("This file contains the diagnostics referred to in the main text. Figures S3, S4 and S11 and Table S14 "
    "come from a preliminary analysis (fixed σ = 0.13 for the exponential and Gaussian media, 150 events); they are "
    "retained for completeness and are superseded in the main text by the per-event test on all events, the "
    "common-grid comparison and the cluster-bootstrap coherence analysis.")

# ================================================================= FIGURES
d.h("Supplementary figures", 1)
fig(F / "val_reciprocity.png",
    "Reciprocity validation in a heterogeneous VTI pink-noise medium (σ = 0.045, ε = 0.05, δ = 0.02, γ = 0.05). "
    "a–c) East, north and vertical velocity at a sensor from a direct double-couple simulation (black) and from the "
    "strain-Green's-tensor reciprocity synthetic (red dashed). Correlation and amplitude ratio are 0.995.")
fig(F / "fig_log5632.png",
    "The 56-32 sonic log. a) Sonic-derived *V*_{P} (grey) and its 50-m block average (black); dashed lines mark the "
    "depths of sensors A and B. b) Distribution of the detrended, 6-m-averaged δln *V*. c) Periodogram of δln *V* "
    "with a power-law fit over 1–200 m. Slope uncertainty is given in Table S15 and main-text Fig. 9.")
fig(F / "fig_misfit_ranking.png",
    "Earlier statistical-misfit ranking (Φ) for a) the synthetic pink truth T1 (σ = 0.045), b) the layered-VTI "
    "truth T2 and c) FORGE (preliminary model set: exponential and Gaussian media only at σ = 0.13). Grey: smooth "
    "media; orange: pink; blue: exponential and Gaussian. Lower is better. Superseded for FORGE by main-text Fig. 8.")
fig(F / "fig_envelopes.png",
    "Median S-coda envelopes aligned on the S peak and normalized by it: a) sensor A, 20–40 Hz; b) sensor A, "
    "40–80 Hz; c) sensor B, 40–80 Hz. Black: FORGE; H homogeneous; A VTI; E50/G50 exponential/Gaussian "
    "(*a* = 50 m, σ = 0.13); P3/P9 pink (σ = 0.09/0.18).")
fig(F / "fig_event_spectrograms.png",
    "Four example events (magnitude −0.61 to +0.24) at sensor B, 4 kHz. a, d, g, j) Horizontal velocity (> 5 Hz) "
    "with P (blue) and S (red). b, e, h, k) Spectrograms to 2 kHz on a common dB scale relative to each event's "
    "maximum; white dashed and dotted lines mark 100 and 50 Hz, the Nyquist frequencies of 200 and 100 sample/s recording. c, f, i, l) S-window "
    "spectrum (black) and pre-event noise (grey); shaded: FD analysis band.")
fig(F / "fig_bandwidth.png",
    "a, b) Median (solid) and interquartile range (dotted) of the S-window amplitude SNR spectrum for sensors A and "
    "B; horizontal line at SNR = 3. c) Peak frequency of the S-window velocity spectrum for all events. The "
    "clustering near 330, 480 and 920 Hz indicates tool or coupling resonances; the share of recorded velocity "
    "energy above 100 Hz is therefore not a source or path property.")
fig(F / "fig_coda_multioctave.png",
    "FORGE coda measures in seven octave bands (25 Hz to 1.3 kHz): a) coda level, b) single-scattering Q_{c}^{-1} "
    "(finite-window diagnostic; at 25 Hz the window spans one to three periods), c) rms S-pulse width. Median and "
    "interquartile range for sensors A and B. Shaded (> 200 Hz): affected by the resonances.")
fig(F / "fig_qc_slope.png",
    "Q_{c}^{-1} frequency slope log_{2}[Q_{c}^{-1}(50 Hz)/Q_{c}^{-1}(25 Hz)] at sensor A, FORGE (black line, grey "
    "95 % bootstrap band) and simulated media (dots). A finite-window diagnostic: data and synthetics are measured "
    "identically, but the absolute values are not medium properties.")
fig(F / "fig_break_test.png",
    "Exploratory high-frequency series (2.5-m grid, sensor B, reduced volume, noise-free synthetics, σ = 0.13). "
    "a) Q_{c}^{-1} normalized at 50 Hz from the 25–50 Hz (sensor A, 5-m grid) and 50–100 Hz (sensor B, 2.5-m grid) "
    "slopes. b) Curvature of Q_{c}^{-1} versus correlation length. c) Coda level at 100 Hz versus correlation "
    "length; black line: FORGE. Grid, domain, sensor and event subset change together between the two octaves, "
    "so the comparison is not a pure frequency effect.")
fig(F / "fig_d2_perturb.png",
    "Location stress test of the inter-event coda coherence (window S + 40 to S + 100 ms): binned medians after perturbing the hypocentres used for binning by isotropic Gaussian errors of 5 m (solid), 20 m (dashed) and 50 m (dotted), for FORGE (black, thick) and the media (colours as in main-text Fig. 7). a) Sensor A. b) Sensor B. Large errors flatten the curves of data and media alike; the levels at ≥ 20 m separation are preserved.")
fig(F / "fig_log_vs_media.png",
    "Vertical fluctuation spectra of the 56-32 sonic log (black) and of the model media over 12–200 m, normalized at "
    "50 m. The rms log difference in the legend is a descriptive distance, not a likelihood.")

# ================================================================= TABLES
d.h("Supplementary tables", 1)
WS = rd("table_coda_window_sensitivity.csv")
mods = [m for m in ("H", "L", "L6", "GRAD3", "FZ", "A", "P1", "P3", "P6", "G50") if m in set(WS.model)]
rows = []
for fe in (15, 10, 20, 60):
    for cw in ("40-100", "50-110", "30-100"):
        r = [f"S + {fe}", f"S + {cw.replace('-', '–')}"]
        for m in mods:
            w = WS[(WS.model == m) & (WS.fit_end_ms == fe) & (WS.coda_ms == cw)].iloc[0]
            r.append(f"{w.med_A:+.2f}/{w.med_B:+.2f}")
        rows.append(r)
tab("Window sensitivity of the per-event coda test", ["Fit end (ms)", "Coda (ms)"] + mods, rows,
    [1.4, 1.6] + [1.35] * len(mods),
    legend="Median log_{10}(predicted/observed late energy), sensor A/B. The fit window starts at P − 10 ms. Rows with "
           "fit end S + 60 ms overlap the coda window (for comparison only).")
MV = rd("table_mt_variants.csv")
rows = []
for m in [m for m in ("H", "L", "L6", "GRAD3", "FZ", "A", "P1", "P3", "P6", "E15", "E50", "G15", "G50", "VSD", "VSO")
          if m in set(MV.model)]:
    g = MV[MV.model == m].set_index("mt")
    rows.append([m] + [f"{g.loc[v,'med_A']:+.2f}/{g.loc[v,'med_B']:+.2f}" for v in ("full", "dev", "dc")]
                + [f"{g.loc[v,'VR']:.2f}" for v in ("full", "dev", "dc")]
                + [f"{g.loc['full','cond_median']:.1f} ({g.loc['full','cond_p90']:.1f})",
                   f"{g.loc['full','pct_ISO_abs']:.0f}/{g.loc['full','pct_CLVD_abs']:.0f}/{g.loc['full','pct_DC']:.0f}"])
tab("Moment-tensor variants: coda prediction, fit, conditioning and decomposition",
    ["Medium", "Full A/B", "Deviatoric A/B", "DC A/B", "VR full", "VR dev", "VR DC", "Cond. median (p90)", "|ISO|/|CLVD|/DC (%)"],
    rows, [1.4, 1.8, 1.8, 1.8, 1.1, 1.1, 1.1, 2.0, 2.2],
    legend="Median log_{10}(pred/obs) at the primary windows. Condition number of the column-normalized design matrix. "
           "Decomposition of the full solutions (median absolute percentages); with two sensors in one well the "
           "non-DC parts are poorly resolved and are not interpreted.")
ST = rd("table_s_timing_subsets.csv")
rows = []
for m in [m for m in ("H", "L", "L6", "GRAD3", "FZ", "A", "P1", "P3", "P6", "E15", "E50", "G15", "G50") if m in set(ST.model)]:
    g = ST[ST.model == m].set_index("subset")
    rows.append([m] + [f"{g.loc[s,'med_log10_pred_obs_A']:+.2f}/{g.loc[s,'med_log10_pred_obs_B']:+.2f}"
                       for s in ("all", "within10ms", "within5ms")])
n = ST[ST.model == "H"].set_index("subset").n
tab("Per-event coda test for well-timed event subsets",
    ["Medium", f"All ({n['all']})", f"S peak ±10 ms ({n['within10ms']})", f"S peak ±5 ms ({n['within5ms']})"],
    rows, [2.0, 3.0, 3.5, 3.5],
    legend="Median log_{10}(pred/obs), sensor A/B. Subsets: events whose observed 40–80 Hz S-envelope peak lies within "
           "the stated interval of the predicted S time on both sensors.")
CB = rd("table_calib_bootstrap.csv").set_index("param"); LO = rd("table_calib_loso.csv").set_index("case")
lab = {"E_m": "East (m)", "N_m": "North (m)", "zA_m": "Depth A (m)", "zB_m": "Depth B (m)", "offA_ms": "Clock A (ms)",
       "offB_ms": "Clock B (ms)", "Vp_m_s": "*V*_{P} (m/s)"}
rows = [[lab[p], f"{CB.loc[p,'estimate']:.1f}", f"{CB.loc[p,'bootstrap_std']:.1f}",
         f"[{CB.loc[p,'lo95']:.1f}, {CB.loc[p,'hi95']:.1f}]",
         f"{LO.loc['without stage 2', p]:.1f}", f"{LO.loc['without stage 3', p]:.1f}"] for p in lab]
tab("Sensor self-calibration: event bootstrap and leave-one-stage-out",
    ["Parameter", "Estimate", "Bootstrap s.d.", "95 % interval", "Without stage 2", "Without stage 3"], rows,
    [2.8, 2.0, 2.2, 3.4, 2.6, 2.6],
    legend="Event bootstrap (500 resamples). Leaving out stage 3 (most events) shifts both sensors ~50 m deeper and "
           "lowers *V*_{P} by 8 %: depth, clock and velocity trade off with the catalogue locations, so the bootstrap "
           "intervals understate the true uncertainty.")
SS = rd("table_source_spectra.csv")
rows = []
for s in ("A", "B"):
    x = SS[SS.sensor == s]
    for lo, hi, lab_ in ((-9, 9, "all"), (-1, -0.3, "−1 to −0.3"), (-0.3, 0.3, "−0.3 to 0.3"), (0.3, 9, "≥ 0.3")):
        y = x[(x.Mw >= lo) & (x.Mw < hi)]
        rows.append([s, lab_, len(y), f"{100*np.isfinite(y.fc_apparent_Hz).mean():.0f}",
                     f"{y.fc_apparent_Hz.median():.0f}", f"{100*(y.fc_apparent_Hz < 80).mean():.0f}",
                     f"{y.disp_slope_30_150.median():.2f}"])
tab("Apparent corner frequencies of S-window displacement spectra",
    ["Sensor", "Magnitude", "Events", "% corner ≤ 300 Hz", "Median corner (Hz)", "% corner < 80 Hz", "Slope 30–150 Hz"],
    rows, [1.4, 2.4, 1.4, 2.4, 2.4, 2.4, 2.2],
    legend="Apparent corner: first frequency above 50 Hz at which the smoothed displacement spectrum falls 6 dB below its "
           "30–50 Hz level (searched to 300 Hz, below the resonances). It includes path attenuation and is therefore "
           "a lower bound on the source corner. Slope: log–log displacement-spectrum slope where amplitude SNR > 3.")
SD = rd("table_source_duration_sensitivity.csv")
rows = [[r.model, "impulsive" if np.isinf(r.brune_fc_Hz) else f"{r.brune_fc_Hz:.0f}", r.sensor,
         f"{r.log10_coda_over_S:.2f}", f"{r.S_width_ms:.1f}"] for r in SD.itertuples()]
tab("Sensitivity of synthetic late energy and S width to source duration",
    ["Medium", "Brune corner (Hz)", "Sensor", "log_{10} coda/S", "S width (ms)"], rows, [2.0, 3.0, 1.6, 3.0, 2.6],
    legend="Homogeneous (H) and pink σ = 0.13 (P6) synthetics with double-couple sources, convolved with Brune "
           "moment-rate functions; median over events; windows as in the main text.")
RP = rd("table_ringing_prominence.csv"); NR = rd("table_null_ringing.csv")
rows = [[r.sensor, r.band_Hz, f"{r.S_spectrum_peak_prominence_dB:.2f}", f"{r.coda_over_S_peak_prominence_dB:.2f}"]
        for r in RP.itertuples()]
rows += [["null test", f"*f*_{{r}} {r.f_r_Hz:.0f} Hz, *Q* {r.Q:.0f}", f"deficit A {r.deficit_A:+.2f}",
          f"deficit B {r.deficit_B:+.2f}"] for r in NR.itertuples()]
tab("Ringing tests: spectral-peak prominence and null-ringing deficits",
    ["Sensor / test", "Band / resonance", "S-spectrum peak prominence (dB)", "Coda/S ratio peak prominence (dB)"],
    rows, [2.6, 4.0, 4.4, 4.4], size=7,
    legend="Upper rows: maximum prominence of peaks in the smoothed median S spectrum and coda/S spectral ratio. Lower "
           "rows: homogeneous synthetics convolved with a single resonance scaled so that its spectral peak equals "
           "the observed 30–120 Hz prominence; deficit = log_{10} of predicted/observed median coda/S.")
STG = rd("table_stage_coda.csv")
rows = [[r.sensor, int(r.stage), int(r.n), f"{r.median_log10_coda_over_S:.2f}", f"{r.median_R_m:.0f}"] for r in STG.itertuples()]
tab("Observed late energy by stimulation stage", ["Sensor", "Stage", "Events", "Median log_{10} coda/S", "Median distance (m)"],
    rows, [2.0, 2.0, 2.0, 4.0, 3.5],
    legend="Stage 2 has only 14 events and differs from stage 3 also in distance; a stage-by-stage test of "
           "stimulation-dependent scattering is not possible with these data.")
FS = rd("table_field_stats.csv")
rows = [[r.medium, f"{r.sigma_lnV_after_truncation:.4f}", f"{r.pct_cells_at_truncation:.2f}",
         f"{r.vp_min:.0f}–{r.vp_max:.0f}", f"{r.vs_min:.0f}–{r.vs_max:.0f}"] for r in FS.itertuples()]
tab("Random-field statistics after truncation", ["Medium", "σ(ln *V*) after truncation", "% cells truncated",
                                                  "*V*_{P} range (m/s)", "*V*_{S} range (m/s)"],
    rows, [3.4, 3.2, 2.6, 3.2, 3.2],
    legend="Truncation at ±3σ (±2.5σ for σ > 0.15) without renormalisation; density is constant (2650 kg/m³).")
CV = rd("table_convergence.csv")
if CV is not None:
    rows = [[r.case, r.measure, r.band, f"{r.median_C5:.3g}", f"{r.median_case:.3g}", f"{r.median_diff:+.3f}",
             f"{r.robust_std_diff:.3f}", int(r.n)] for r in CV.itertuples()]
    tab("Numerical convergence of coda measures", ["Case", "Measure", "Band", "Median 5-m reference", "Median case",
                                                  "Median difference", "Robust s.d. of difference", "Events"],
        rows, [1.4, 3.6, 1.6, 1.9, 1.7, 1.9, 2.0, 1.2], size=7,
        legend="C25: 2.5-m grid; CBIG: domain enlarged by 150 m on every side; CSP: sponge of 40 instead of 20 nodes. "
               "Same band-limited medium (σ = 0.13), sensor B, events and double-couple mechanisms as the 5-m "
               "reference C5. Differences are per event.")
else:
    pending_tab("Numerical convergence of coda measures", "S26")
SC = rd("table_sigma_curves.csv"); SE = rd("table_sigma_eff.csv")
if SC is not None and SE is not None:
    rows = [[r.family, f"{r.sigma:.3f}", int(r.n_realisations), f2(r.med_log10_pred_obs_A), f2(r.med_log10_pred_obs_B)]
            for r in SC.itertuples()]
    rows += [[r.family, "σ_{eff}", "", f"{r.sigma_eff_A:.3f}" if r.sigma_eff_A > 0 else "outside grid",
              f"{r.sigma_eff_B:.3f}" if r.sigma_eff_B > 0 else "outside grid"] for r in SE.itertuples()]
    tab("Median predicted/observed late energy versus σ, and σ_{eff}", ["Family", "σ", "Realisations",
                                                                        "Median log_{10}(pred/obs), A", "Median, B"],
        rows, [3.0, 1.8, 2.2, 3.6, 2.6], size=7,
        legend="Averaged over realizations; σ_{eff} interpolated linearly in ln σ at the zero crossing.")
else:
    pending_tab("Median predicted/observed late energy versus σ, and σ_{eff}", "S29")
PB = rd("table_phi_bootstrap.csv")
if PB is not None:
    rows = [[r.model, r.family, "" if pd.isna(r.sigma) else f"{r.sigma:.3f}", "" if pd.isna(r.shape) else f"{r.shape:g}",
             f"{r.PHI:.3f}", f"[{r.PHI_lo95:.3f}, {r.PHI_hi95:.3f}]"] for r in PB.itertuples()]
    tab("Statistical misfit Φ with event-bootstrap 95 % intervals", ["Medium", "Family", "σ", "*p* or *a* (m)", "Φ",
                                                                   "95 % interval"],
        rows, [4.2, 2.2, 1.4, 2.0, 1.6, 3.2], size=7,
        legend="200 event-bootstrap resamples (identical for every medium). Φ as defined in the main text.")
else:
    pending_tab("Statistical misfit Φ with event-bootstrap 95 % intervals", "S20")
FB = rd("table_family_bootstrap.csv"); HS = rd("table_holdout_summary.csv")
if FB is not None and HS is not None:
    hs = HS.set_index("family")
    rows = [[r.family, f"{r.best_PHI:.3f}", f"{r.boot_median:.3f} [{r.lo95:.3f}, {r.hi95:.3f}]", f"{r.pct_boot_wins:.0f}",
             f"{hs.loc[r.family,'median_PHI_test']:.3f} [{hs.loc[r.family,'min_PHI_test']:.3f}, {hs.loc[r.family,'max_PHI_test']:.3f}]"
             if r.family in hs.index else "—"] for r in FB.sort_values("best_PHI").itertuples()]
    tab("Family comparison: best member, bootstrap and holdout", ["Family", "Best Φ", "Bootstrap median [95 %]",
                                                                "% bootstrap wins", "Held-out Φ median [range]"],
        rows, [2.4, 1.6, 4.2, 2.4, 4.4],
        legend="Best member per bootstrap resample; holdout: best (σ, shape, seed) chosen on a training half and scored "
               "on the other half (20 random splits + near/far split, both directions).")
else:
    pending_tab("Family comparison: best member, bootstrap and holdout", "S20")
CF = rd("table_confusion_runs.csv"); DF = rd("table_data_families.csv").set_index("family")
rows = [[f"{r.truth_family} ({r.truth})", f"{r.score_pink:.2f}", f"{r.score_exponential:.2f}", f"{r.score_Gaussian:.2f}",
         r.selected + (" ✓" if r.selected == r.truth_family else " ✗")] for r in CF.itertuples()]
rows.append(["FORGE data (preliminary set)", f"{DF.loc['pink','mean']:.2f}", f"{DF.loc['exponential','mean']:.2f}",
             f"{DF.loc['Gaussian','mean']:.2f}", "—"])
tab("Equal-σ family test (σ = 0.13): family scores for nine synthetic truths",
    ["Truth", "Pink", "Exponential (a = 15 m)", "Gaussian (a = 50 m)", "Selected"], rows, [4.4, 2.0, 3.2, 3.2, 3.0],
    legend="Family score = mean Φ over candidate realizations whose seed differs from the truth. 5 of 9 correct; "
           "P(≥ 5 of 9 | p = 1/3) = 0.145. Pink and exponential fields built from the same seed are strongly correlated.")
LSV = rd("table_log_slope_variants.csv"); LSM = rd("table_log_spectral_models.csv")
rows = [[r.detrend, r.psd, r.segment, f"{r.beta_12_200:.2f}", f"{r.beta_1_200:.2f}"] for r in LSV.itertuples()]
rows += [[r.band, r.model, r.note if isinstance(r.note, str) else "", f"rms {r.rms_log:.3f}", f"ΔAIC {r.dAIC:.1f}"]
         for r in LSM.itertuples()]
tab("Sonic-log spectral slope and spectral-model comparison", ["Detrend / band", "PSD / model", "Segment / parameter",
                                                               "β (12–200 m) / rms", "β (1–200 m) / ΔAIC"],
    rows, [2.6, 4.2, 4.0, 2.4, 2.4], size=7,
    legend="Upper rows: least-squares slope of ln PSD versus ln *k* (unbinned). Lower rows: fits to log-binned "
           "multitaper spectra (equal weight per log bin), which give a flatter power-law slope over 12–200 m "
           "(0.54) than the unbinned fit; the slope over this 1.2-decade band is poorly determined.")
CP = rd("table_coda_prediction.csv")
rows = []
for m in CP.model.unique():
    a = CP[(CP.model == m) & (CP.events == "all")].iloc[0]; t = CP[(CP.model == m) & (CP.events == "top150")].iloc[0]
    rows.append([m, f"{a.med_log10_pred_obs_A:+.2f}/{a.med_log10_pred_obs_B:+.2f}", f"{a.pct_within3_A:.0f}/{a.pct_within3_B:.0f}",
                 f"{a.VR:.2f}", f"{t.med_log10_pred_obs_A:+.2f}/{t.med_log10_pred_obs_B:+.2f}"])
tab("Per-event coda test for every medium and seed", ["Medium", f"Median A/B (all {int(CP.n.max())})", "% within ×3 A/B",
                                                      "VR", "Median A/B (150 highest SNR)"],
    rows, [4.4, 3.0, 2.6, 1.4, 3.6], size=7,
    legend="Primary windows (fit P − 10 ms to S + 15 ms; coda S + 40 to S + 100 ms). Common-grid media are named "
           "X_<family>_<σ>_<p or a>_<seed>.")
DA = rd("table_d2_allpairs.csv")
if DA is not None:
    rows = [[r.model, r.window.replace("-", "–"), r.sensor, f"{r.median_coh_all:.3f}", f"[{r.lo95:.3f}, {r.hi95:.3f}]",
             int(r.n_pairs)] for r in DA.itertuples()]
    tab("Inter-event coda coherence over all pairs within 160 m, by window",
        ["Medium", "Window", "Sensor", "Median coherence", "95 % interval", "Pairs"], rows, [4.2, 3.0, 1.4, 2.6, 3.4, 1.6],
        size=7, legend="Cluster-bootstrap intervals over events (500 resamples). Pairs subsampled to 20 000 where more "
                       "exist. The noise window (S + 550 to S + 610 ms) contains no coda and gives the coherence of "
                       "noise alone; every synthetic record contains its event's own recorded noise.")
DP = rd("table_d2_perturb.csv")
if DP is not None:
    rows = []
    for (m, w, g), x in DP.groupby(["model", "window", "sensor"], sort=False):
        if m not in ("data", "H", "L6", "FZ", "P6", "E15", "G50"):
            continue
        r = [m, w.replace("-", "–"), g]
        for sp in (5, 20, 50):
            y = x[(x.perturb_m == sp) & (x.d_lo_m >= 20)].median_coh.values
            r.append(" / ".join(f"{v:.2f}" for v in y))
        rows.append(r)
    tab("Location stress test: binned coherence (20–40, 40–80, 80–160 m) after hypocentre perturbation",
        ["Medium", "Window", "Sensor", "σ_{loc} = 5 m", "σ_{loc} = 20 m", "σ_{loc} = 50 m"], rows,
        [1.8, 2.6, 1.3, 3.4, 3.4, 3.4], size=7,
        legend="Mean over 50 perturbation draws of the median coherence in the bins 20–40, 40–80 and 80–160 m; "
               "perturbations applied identically to data and media.")
print(f"figures S1-S{FIGN[0]}, tables S1-S{TABN[0]}")
print("PENDING:", sorted(set(PENDING)))
d.save(P / "manuscript" / "Paper2_GJI_Supporting_Information.docx")
