"""Paper 2, Geophysical Journal International (GJI) version (2026-10-01).
All numbers are read from results/*.csv.
Run:  py -3.12 build_paper2_gji.py            -> Paper2_GJI.docx (text, tables and captions)
      EMBED=1 py -3.12 build_paper2_gji.py    -> Paper2_GJI_review.docx (figures embedded)
Supporting Information: build_supporting_gji.py.
"""
import os, sys
from pathlib import Path
from math import comb
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gji_doc import GDoc  # noqa: E402
from refs_gji import REFS  # noqa: E402

R = P / "results"; F = P / "figs"
PENDING = []


def rd(f, **kw):
    fn = R / f
    if not fn.exists():
        PENDING.append(f); return None
    return pd.read_csv(fn, **kw)


def pend(label):
    PENDING.append(label); return f"[PENDING: {label}]"


K = rd("key_numbers.csv").set_index("key").value
f2 = lambda x: f"{x:.2f}"
x3 = np.log10(3)

# ---------------------------------------------------------------- per-event coda prediction (S19/S27)
CP = rd("table_coda_prediction.csv")
CPa = CP[CP.events == "all"].set_index("model")
CPt = CP[CP.events == "top150"].set_index("model")
has = lambda m: m in CPa.index
cpa = lambda m, c: CPa.loc[m, c]
SMOOTH = [m for m in ("H", "L", "GRAD3", "A") if has(m)]          # no structure below ~50 m
DET = [m for m in ("L6", "FZ") if has(m)]                          # deterministic small-scale structure
NONSTOCH = SMOOTH + DET
smA = CPa.loc[SMOOTH, "med_log10_pred_obs_A"]; smB = CPa.loc[SMOOTH, "med_log10_pred_obs_B"]
smwA = CPa.loc[SMOOTH, "pct_within3_A"]; smwB = CPa.loc[SMOOTH, "pct_within3_B"]
dtA = CPa.loc[DET, "med_log10_pred_obs_A"]; dtB = CPa.loc[DET, "med_log10_pred_obs_B"]
STOCH = [m for m in CPa.index if m not in NONSTOCH]
fac = lambda v: 10 ** (-v)
n_ev = int(CPa.loc["H", "n"])
REAL = [["P1", "R1", "R2", "T1"], ["P3", "R3"], ["P6", "R6"], ["P9", "R5"]]
REAL += [[f"X_{fm}_{s}_{a}_11", f"X_{fm}_{s}_{a}_33"] for fm in ("exp", "gau") for s in ("0.09", "0.18") for a in (15, 50)]
def spread(groups):
    v = [CPa.loc[[m for m in g if has(m)], c].max() - CPa.loc[[m for m in g if has(m)], c].min()
         for g in groups if sum(has(m) for m in g) > 1 for c in ("med_log10_pred_obs_A", "med_log10_pred_obs_B")]
    return max(v) if v else np.nan


REAL += [["P8", "X_pink_0.09_1.8_33"], ["X_pink_0.13_1.8_11", "X_pink_0.13_1.8_33"], ["P11", "X_pink_0.18_1.8_33"],
         ["E15", "E15b"], ["E50", "X_exp_0.13_50_33"], ["G15", "X_gau_0.13_15_33"], ["G50", "G50b"]]
fam_of = lambda g: ("pink" if g[0].startswith(("P", "R", "T", "X_pink")) else
                    "exp15" if ("exp" in g[0] and "_15_" in g[0]) or g[0].startswith("E15") else
                    "exp50" if ("exp" in g[0] and "_50_" in g[0]) or g[0].startswith("E50") else
                    "gau15" if ("gau" in g[0] and "_15_" in g[0]) or g[0].startswith("G15") else "gau50")
SPR = {f: spread([g for g in REAL if fam_of(g) == f]) for f in ("pink", "exp15", "exp50", "gau15", "gau50")}
rspread = max(SPR["pink"], SPR["exp15"])          # scatter for the small-scale-rich families
# sensitivity of the median ratio to ln(sigma), pink p = 1.5, sigma 0.045 -> 0.09 (mean of sensors)
rslope = np.mean([(cpa("P3", c) - cpa("P1", c)) / np.log(2) for c in ("med_log10_pred_obs_A", "med_log10_pred_obs_B")])
WS = rd("table_coda_window_sensitivity.csv")
MV = rd("table_mt_variants.csv")
ST = rd("table_s_timing_subsets.csv")
SEL = rd("table_event_selection.csv")
SELC = rd("table_selection_coda.csv")
NR = rd("table_null_ringing.csv")
RP = rd("table_ringing_prominence.csv")
SD = rd("table_source_duration_sensitivity.csv")
SS = rd("table_source_spectra.csv")
CB = rd("table_calib_bootstrap.csv").set_index("param")
LO = rd("table_calib_loso.csv").set_index("case")
STG = rd("table_stage_coda.csv")
FS = rd("table_field_stats.csv").set_index("medium")
D2 = rd("table_d2_uncertainty.csv")
D2A = rd("table_d2_allpairs.csv"); D2W = rd("table_d2_windows.csv"); D2P = rd("table_d2_perturb.csv")
W15, W40, W50, WNZ = "S+15-100 ms", "S+40-100 ms", "S+50-110 ms", "S+550-610 ms"
da = lambda m, w, g, c="median_coh_all": D2A[(D2A.model == m) & (D2A.window == w) & (D2A.sensor == g)][c].iloc[0]
def dar(ms, w, g):
    v = [da(m, w, g) for m in ms if ((D2A.model == m) & (D2A.window == w)).any()]
    return f"{min(v):.2f}–{max(v):.2f}" if len(v) > 1 else f"{v[0]:.2f}"
D2SM, D2DT, D2ST = ["H", "L", "GRAD3", "A"], ["L6", "FZ"], ["P3", "P6", "E15", "G50"]
dnoise = lambda g: da("data (noise window)", WNZ, g)
_cases = [(m, w, g) for m in ["H", "L", "GRAD3", "A", "FZ"] for w in (W40, W50) for g in "AB"]
nsep = sum(da(m, w, g, "lo95") > da("data", w, g, "hi95") for m, w, g in _cases); ntot = len(_cases)
LS = rd("table_log_slope_summary.csv").iloc[0]
LSV = rd("table_log_slope_variants.csv")
LSM = rd("table_log_spectral_models.csv")
CF = rd("table_confusion_runs.csv")
MT1 = rd("table_mt_T1.csv").set_index("model")
M1 = rd("table_misfit_T1_forge.csv").set_index("model")
M2 = rd("table_misfit_T2_forge.csv").set_index("model")
FC = rd("table_field_corr.csv").set_index("medium").corr_with_pink_seed33
BW = rd("table_bandwidth.csv").set_index("sensor")
BT = rd("table_break_test.csv").set_index("model")
SE = rd("table_sigma_eff.csv"); SC = rd("table_sigma_curves.csv")
PB = rd("table_phi_bootstrap.csv"); FB = rd("table_family_bootstrap.csv"); HS = rd("table_holdout_summary.csv")
HO = rd("table_holdout.csv"); CV = rd("table_convergence.csv")
if CV is not None:
    cvs = lambda case, b: CV[(CV.measure.str.startswith("log10 coda/S")) & (CV.case == case) & (CV.band == b)].median_diff.iloc[0]
    cvd = lambda case, b="40-80 Hz": CV[(CV.measure.str.startswith("decay")) & (CV.case == "CBIG") & (CV.band == b)]["median_C5" if case == "C5" else "median_case"].iloc[0]
    fdom = np.exp(cvs("CBIG", "40-80 Hz") / rslope)          # sigma_eff over-estimate factor from the domain test

ws = lambda m, fe, cw, c: WS[(WS.model == m) & (WS.fit_end_ms == fe) & (WS.coda_ms == cw)][c].iloc[0]
mv = lambda m, v, c: MV[(MV.model == m) & (MV.mt == v)][c].iloc[0]
stt = lambda m, s, c: ST[(ST.model == m) & (ST.subset == s)][c].iloc[0]
d2 = lambda m, s, lo, c="median_coh": D2[(D2.model == m) & (D2.sensor == s) & (D2.d_lo_m == lo)][c].iloc[0]
sd = lambda m, fc, s: SD[(SD.model == m) & (SD.brune_fc_Hz == fc) & (SD.sensor == s)].log10_coda_over_S.iloc[0]
n_ok = int((CF.truth_family == CF.selected).sum())
p_binom = sum(comb(9, k) * (1 / 3) ** k * (2 / 3) ** (9 - k) for k in range(n_ok, 10))
pl = CF.truth_family.isin(["pink", "exponential"])
n_pe = int(((CF.truth_family == CF.selected) & pl).sum())
n_pl = int((pl & CF.selected.isin(["pink", "exponential"])).sum() + (~pl & (CF.selected == "Gaussian")).sum())
rexp = FC["exponential a=15 m (seed 33)"]
sA = SS[SS.sensor == "A"]; sB = SS[SS.sensor == "B"]
ms = lambda T, m, b: LSM[(LSM.band == b) & (LSM.model == m)][T].iloc[0]


def sigeff(fam, g):
    """sigma_eff text for family fam, sensor g ('A'/'B')."""
    if SE is None:
        return pend("sigma_eff")
    v = SE.set_index("family").loc[fam, f"sigma_eff_{g}"]
    if pd.isna(v):
        return "n/a"
    if v < -0.5:
        return f"> {(-v - 1):.2f}"
    if v < 0:
        return f"< {-v:.2f}"
    return f"{v:.3f}"


RELEASE = os.environ.get("RELEASE", "v1.0-submission")
TITLE = ("Late S-wave coda of Utah FORGE downhole microearthquakes: elastic tests of smooth, deterministic and "
         "stochastic media")
RUNNING = "Late S coda at Utah FORGE"

d = GDoc(REFS)
d.embed = bool(os.environ.get("EMBED"))
d.p(f"**{TITLE}**", align="center", size=14)
d.p("Wiwit Suryanto^{1} and Peter Leary^{2}", align="center")
for a in ["^{1} Department of Physics, Faculty of Mathematics and Natural Sciences, Universitas Gadjah Mada, "
          "Sekip Utara BLS 21, Yogyakarta 55281, Indonesia. E-mail: ws@ugm.ac.id",
          "^{2} Geoflow Imaging, Auckland 1010, New Zealand"]:
    d.p(a, align="left", size=11)
d.p("Corresponding author: Wiwit Suryanto (ws@ugm.ac.id)", align="left", size=11)
d.p(f"Abbreviated title: {RUNNING}", align="left", size=11)

# ================================================================= SUMMARY
d.p("**SUMMARY**", align="left")
D2M = set(D2.model)
d2det = [m for m in NONSTOCH if m in D2M]
sev = pd.concat([SE.sigma_eff_A, SE.sigma_eff_B]); sev = sev[sev > 0]
se_lo, se_hi = sev.min(), sev.max()
se_txt = f"{se_lo:.2f}–{se_hi:.2f}"
summary = (
    f"Downhole records of microearthquakes induced by the April 2022 stimulation at Utah FORGE contain long S-wave "
    f"codas. We compare {n_ev} events recorded at 4 kHz by two three-component geophones in well 56-32 with 3-D "
    f"elastic finite-difference wavefields computed by reciprocity for smooth media (homogeneous, sonic-log layering "
    f"in 50-m blocks, 3-D gradients, transverse isotropy), simple deterministic small-scale media (sonic-log layering "
    f"at 6-m resolution, thin fracture zones) and lognormal random media with power-law, exponential and Gaussian "
    f"spectra. Moment tensors fitted to the direct waves predict the 40–80 Hz energy 40–100 ms after the S arrival. "
    f"The smooth media under-predict this late energy by median factors of {fac(smA.max()):.0f}–{fac(smA.min()):.0f} "
    f"and {fac(smB.max()):.0f}–{fac(smB.min()):.0f} at the two sensors, a deficit that is robust to the analysis "
    f"windows, the source parametrization, event timing and a worst-case instrument resonance. The deterministic "
    f"small-scale media reduce the deficit to factors of {fac(dtB.max()):.0f}–{fac(dtA.min()):.0f}, and the random "
    f"media span the observations. In the early coda, records of neighbouring events decorrelate with separation as "
    f"in the random media, whereas the smooth and deterministic media stay coherent; in the late coda the coherence "
    f"is limited by noise and separates the data only from the smooth and fracture-zone media. Within isotropic "
    f"lognormal media the late energy implies a fluctuation strength (standard deviation of ln V) of {se_txt}, "
    f"probably an upper bound and of the same order as the sonic-log value of {K.log_sig6:.3f}. Neither the "
    f"waveforms nor the log resolve the spectral family at these frequencies. The late S energy therefore requires "
    f"structure at the 5–10 m scale, unless an unmodelled site response of comparable size is present. Synthetic "
    f"tests show that comparing waveforms with an arbitrary random realization favours smooth media, so "
    f"heterogeneous media should be compared through wavefield statistics.")
d.p(summary)
NW_ABS = len(summary.split())
d.p("**Key words:** Coda waves; Wave scattering and diffraction; Computational seismology; Induced seismicity; "
    "Body waves", align="left")

# ================================================================= INTRODUCTION
d.h("Introduction", 1)
d.p("Microseismic monitoring of engineered geothermal system (EGS) stimulation relies on a velocity model "
    "[@majer2007]. In routine practice that model is homogeneous, layered or weakly anisotropic, smooth enough for "
    "ray theory to locate events and for a few source parameters to fit the direct waves. Everything after the "
    "direct S wave, the coda, is usually left unmodelled. Coda waves are generated by scattering on small-scale "
    "heterogeneity [@aki1975; @sato2012]. Their envelopes and attenuation have been modelled with random media "
    "whose fluctuations have exponential, Gaussian or power-law autocorrelation functions, by 2-D finite "
    "differences [@frankel1986; @roth1993; @hestholm1994], by 3-D simulations that include intrinsic attenuation "
    "[@takemura2017] and by envelope theory [@shapiro1993; @saito2003; @sato2012]. Which heterogeneity produces "
    "the coda of a particular reservoir, and how strong it must be, remains much less clear, because near-source "
    "records with enough bandwidth to test it are rare.")
d.p("Borehole logs provide an independent constraint. Sonic and resistivity logs in crystalline and sedimentary "
    "crust show power-law fluctuation spectra, *S*(*k*) ∝ *k*^{-β} with β ≈ 1–1.5, over metre-to-kilometre scales "
    "[@leary1991; @wu1994; @holliger1996], and the heterogeneity inferred from seismic scattering has likewise been "
    "described as fractal [@wu1985]. Core porosity and permeability co-vary lognormally [@leary1997; @leary2002], so "
    "a lognormal, spatially correlated field with a power-law spectrum (a \"pink-noise\" field) is a natural "
    "candidate for both flow and wave propagation; such a field has been used to describe the flow data of the "
    "FORGE reservoir [@leary2026]. Whether the seismic wavefield of the FORGE microearthquakes (MEQs) requires strong "
    "small-scale heterogeneity, and whether it can distinguish a power-law spectrum from spectra with a "
    "correlation length, has not been tested.")
d.p("A second observable is the similarity of the coda of neighbouring events. Coda-wave interferometry exploits "
    "the sensitivity of multiply scattered waves to small changes of the medium or of the source position "
    "[@snieder2002; @planes2014; @margerin2016], and coda decorrelation has recently been used to follow fracture "
    "growth in EGS experiments [@pradhan2025]. Here we use the decorrelation of the coda between neighbouring "
    "events as a descriptor of how the scattering structure varies in space.")
d.p("We use the 4-kHz records of two downhole sensors in well 56-32 during the April 2022 stimulation and 3-D "
    "elastic finite-difference (FD) wavefields for several hundred event–sensor paths per model. We ask three "
    "questions, in decreasing order of what the data can answer. (1) Do smooth media, or media with deterministic "
    "small-scale structure (thin layers, planar fracture zones), reproduce the late S-wave energy and its "
    "inter-event coherence? (2) If not, what fluctuation strength do random media require? (3) Do the waveforms "
    "prefer a power-law spectrum over exponential or Gaussian spectra? The first question is the main test, and we "
    "subject it to controls for window overlap, source, instrument ringing, event selection and numerical "
    "convergence. The second and third are answered with their limitations, and the third is answered negatively; "
    "Section 5.4 discusses what would be needed to resolve it.")

# ================================================================= DATA
d.h("Data", 1)
d.h("Records and event selection", 2)
sel = SEL.set_index("step")
rows = [["GES catalogue, April 2022 stimulation of 16A(78)-32", "36 641", "—"],
        ["60-s SEG-2 files read (18–24 April 2022; 4 kHz, 6 channels)", "3556 files", "—"]]
lab = {"extracted (M >= -1, largest 800)": "Largest events with *M* ≥ −1 whose 1-s window lies inside a file (extraction cap, ordered by magnitude)",
       "inside FD box (either sensor)": "Clean P onset on at least one sensor, hypocentre inside the FD interior",
       "both sensors (coda test)": "Clean P onsets on both sensors: per-event coda test, D2, Φ",
       "150 highest P-SNR": "150 highest P-wave SNR (sensitivity check only)"}
for k, v in lab.items():
    r = sel.loc[k]
    rows.append([v, f"{int(r.n)}", f"{r.M_med:.2f} [{r.M_p10:.2f}, {r.M_p90:.2f}]; {r.R_med_m:.0f} [{r.R_p10_m:.0f}, "
                 f"{r.R_p90_m:.0f}] m; {r.pct_stage3:.0f} %"])
d.table("Event selection", ["Step", "Events", "Magnitude; distance to sensor mid-point; share of stage 3"], rows,
        [8.0, 1.8, 6.2],
        legend="Median [10th, 90th percentile]. Stage 3 dominates at every step; stage 2 contributes 14 of the "
               f"{int(sel.loc['both sensors (coda test)'].n)} events used, and stage 1 none (Section 2.1).")
sc = SELC.set_index("subset")
d.p(f"Geo-Energie Suisse recorded the April 2022 stimulation of well 16A(78)-32 in 60-s SEG-2 files from "
    f"monitoring well 56-32, which carried a two-level geophone string from stage 2 onwards [@whidden2023], and "
    f"published a catalogue of 36 641 located events [@dyer2022]. Fig. 1 shows the "
    f"geometry, Fig. 2 example records, and Table 1 the selection steps. The only magnitude criterion is the extraction cap of 800 events; no step uses the coda. "
    f"The per-event test uses all {int(sel.loc['both sensors (coda test)'].n)} events with clean P onsets on both "
    f"sensors; the 150 events with the highest P-wave SNR serve as a sensitivity check. They have a slightly weaker, "
    f"not stronger, observed coda than the remaining events (median log_{{10}} coda/S "
    f"{sc.loc['150 highest P-SNR','obs_log10_coda_S_A']:.2f} against {sc.loc['remaining 266','obs_log10_coda_S_A']:.2f} "
    f"at sensor A), so the selection did not inflate the coda.")
d.h("Sensor calibration", 2)
d.p(f"The file headers carry no receiver coordinates. Automatic P onsets [@maeda1985] on the two three-component "
    f"sensors (channels 1–3, sensor A; 4–6, sensor B) were inverted with the catalogue hypocentres for sensor "
    f"position (both sensors in one vertical well), clock offsets and *V*_{{P}} (L1 misfit). The solution places "
    f"the sensors at depths of {CB.loc['zA_m','estimate']:.0f} m (A) and {CB.loc['zB_m','estimate']:.0f} m (B) with "
    f"*V*_{{P}} = {CB.loc['Vp_m_s','estimate']:.0f} m/s. The deeper sensor lies within 12 m of the reported "
    f"reservoir depth of the string (8315 ft, 2534 m; <@whidden2023>), and *V*_{{P}} agrees with the 5.84 km/s "
    f"used for the FORGE granitoid in tomographic models [@nakata2023]. Event-bootstrap 95 % intervals are narrow "
    f"(±{(CB.loc['zA_m','hi95']-CB.loc['zA_m','lo95'])/2:.0f} m for A, ±{(CB.loc['Vp_m_s','hi95']-CB.loc['Vp_m_s','lo95'])/2:.0f} m/s "
    f"for *V*_{{P}}), but they understate the true uncertainty, because the catalogue locations were themselves "
    f"computed with a smooth velocity model. Leaving out stimulation stage 3 moves both sensors "
    f"{LO.loc['without stage 3','zA_m']-CB.loc['zA_m','estimate']:.0f} m deeper and lowers *V*_{{P}} by "
    f"{100*(1-LO.loc['without stage 3','Vp_m_s']/CB.loc['Vp_m_s','estimate']):.0f} %, which shows that depth, "
    f"clock and velocity trade off structurally. We therefore regard sensor depths as uncertain by ~50 m and "
    f"*V*_{{P}} by ~8 %, and test the effect of timing errors on the coda measurements directly (see below). The "
    f"P travel-time residuals of neighbouring events are correlated up to ~40 m separation (Fig. 2c), which may "
    f"reflect structure or correlated catalogue location error; we do not interpret them further. "
    f"Channels 1 and 4 carry only ~3 % of the P energy and are not used; channels 2–3 and 5–6 behave as horizontal "
    f"pairs, and all measurements use rotation-invariant horizontal energy. Sensor A is the shallower and sensor B "
    f"the deeper of the two.")
d.h("Sonic log", 2)
bmin, bmax = LSV.beta_12_200.min(), LSV.beta_12_200.max()
b1min, b1max = LSV.beta_1_200.min(), LSV.beta_1_200.max()
d.p(f"The 56-32 sonic log (1950–2780 m, calliper-screened) gives a detrended ln *V* with σ = {K.log_sig6:.3f} at "
    f"6-m and {K.log_sig50:.3f} at 50-m averaging. Its spectral slope depends on band and estimator (Fig. 9). Over "
    f"12–200 m the multitaper estimate is β = {LS.beta_ref:.2f} ± {LS.beta_se_jackknife:.2f} (taper jackknife), and "
    f"across detrending orders, periodogram, Welch and multitaper estimators and three log segments it ranges from "
    f"{bmin:.2f} to {bmax:.2f}. Over 1–200 m it is {b1min:.1f}–{b1max:.1f}. On log-binned spectra a power law is the "
    f"best of four models over 12–200 m, but exponential and Gaussian autocorrelation spectra with correlation "
    f"lengths of 2–4 m lie within ΔAIC = {ms('dAIC','exponential ACF (Lorentzian)','12-200 m'):.1f} and "
    f"{ms('dAIC','Gaussian ACF','12-200 m'):.1f}. Over 1–200 m the Gaussian spectrum is excluded "
    f"(ΔAIC = {ms('dAIC','Gaussian ACF','1-200 m'):.0f}), but the exponential one is not "
    f"(ΔAIC = {ms('dAIC','exponential ACF (Lorentzian)','1-200 m'):.1f}). The log is thus compatible with power-law "
    f"fluctuations, but a single 830-m vertical log cannot establish scale invariance or exclude a short "
    f"correlation length. It is also a sample of intact rock along one line, not of the stimulated volume that the "
    f"MEQ waves traverse.")
d.figure("Utah FORGE 2022 geometry: stimulation MEQs, 56-32 sensors and FD volume",
         "a) Map view and b) east–depth section of the 416 analysed MEQs (orange: stage 3; blue: stage 2) and of "
         "the other catalogued events inside the FD interior (grey), with the calibrated sensors A and B in well 56-32 "
         "(red triangles) and the FD interior (dotted). Coordinates relative to the 16A(78)-32 wellhead. "
         "c) Hypocentral distances of the 416 events to sensors A and B.",
         F / "fig_geometry_rev.png")
d.figure("FORGE 56-32 records and S envelopes",
         "a, d) Horizontal records (one component, 40–80 Hz) of 40 events ordered by hypocentral distance for sensors "
         "A and B; the dashed line is the predicted S time. b, e) Median (solid) and interquartile range (dotted) of "
         "the horizontal S envelope, 40–80 Hz. c) Correlation of P travel-time residuals versus inter-event "
         "separation.", F / "fig_data_overview.png")

# ================================================================= METHODS
d.h("Methods", 1)
d.h("3-D elastic finite differences and reciprocity", 2)
d.p(f"We solve the velocity–stress equations on a staggered grid [@virieux1986] with fourth-order spatial and "
    f"second-order temporal differences [@levander1988] for a transversely isotropic medium with a vertical symmetry "
    f"axis (VTI; stiffnesses from *V*_{{P}}, *V*_{{S}}, ρ and the Thomsen parameters ε, δ, γ; <@thomsen1986>). The treatment of "
    f"anisotropy on the staggered grid follows {{@igel1995}}, and the MATLAB implementation "
    f"grew out of the acoustic codes of {{@igel2016}}, extended to 3-D elastic VTI media, moment-tensor and "
    f"point-force sources and strain-rate recording. A Cerjan sponge [@cerjan1985] absorbs outgoing energy. The "
    f"grid is {int(K.nx)}×{int(K.ny)}×{int(K.nz)} nodes at {K.dx_m:.0f} m (a 0.36 × 0.48 × 0.54 km interior plus "
    f"20-node sponges), with Δ*t* = {K.dt_s*1e3:.2f} ms, {int(K.nt)} steps and a {K.f0_Hz:.0f}-Hz Ricker source. The "
    f"analysis band is 20–80 Hz, with ≥ 5.4 grid points per S wavelength at 80 Hz even for the slowest truncated "
    f"S velocity of the σ = 0.18 media (media with σ = 0.22 and 0.25, which fall below 5 points, are shown only in "
    f"Supporting Information). Three point-force runs per sensor record the strain-rate "
    f"tensor at every event node, so that by reciprocity [@zhao2006] the velocity at sensor component *n* from a "
    f"moment tensor **M** at **x**_{{e}} is")
d.eq(r"v_{n}(\mathbf{x}_{s},t) = M_{pq}\int_{0}^{t}\dot{\varepsilon}^{(n)}_{pq}(\mathbf{x}_{e},\tau)\,d\tau .")
d.p("Because only horizontal components are used, the two horizontal-force runs per sensor suffice. Against "
    "direct moment-tensor simulations in a heterogeneous VTI medium, the reciprocity synthetics agree with "
    "correlation 0.995 and amplitude ratio 0.995 (Supporting Information, Fig. S1).")
d.p("Coda amplitudes in a small domain could be affected by the grid, the domain edges and the sponge. We therefore "
    "simulated one band-limited stochastic medium (σ = 0.13, power-law spectrum with a 780-m outer scale, low-passed "
    "at 0.06–0.08 cycles/m so that it is identical on every grid) four times with sensor B, the same "
    "358 events and the same double-couple mechanisms: on the reference 5-m grid, on a 2.5-m grid, in a domain "
    "extended by 150 m on every side, and with a sponge of 40 instead of 20 nodes. The coda measures of the three "
    "variants are compared event by event with the reference.")
d.h("Media", 2)
d.p("All media have constant density ρ = 2650 kg/m³ and background *V*_{P} = 5808 m/s, *V*_{P}/*V*_{S} = 1.73, "
    "except where stated. Four smooth media contain no structure below ~50 m: a homogeneous medium (H), a layered "
    "medium from the 56-32 sonic log averaged in 50-m blocks (L), a medium with smooth 3-D gradients of ±3 % "
    "vertically and ±2 % horizontally across the box (GRAD3), and a homogeneous VTI medium with Thomsen parameters "
    "inverted from the P travel times (A). Two deterministic media contain small-scale structure: the sonic log "
    "sampled at its 6-m resolution on the 5-m grid (L6), and 12 planar fracture zones, one cell (5 m) thick, about "
    "240 m tall, striking N25°E and dipping 75°, with *V*_{P} reduced by 10 % and *V*_{S} by 20 % (FZ, which also "
    "varies *V*_{P}/*V*_{S}).")
d.p("The stochastic media are lognormal, *V*(**x**) = *V̄* exp[σ *f*(**x**)], where *f* is a zero-mean, "
    "unit-variance Gaussian random field obtained by filtering white noise *ŵ* in the wavenumber domain. For the "
    "power-law (pink) field the filter is")
d.eq(r"\hat{f}(\mathbf{n}) \propto \frac{\hat{w}(\mathbf{n})}{\left(1+|\mathbf{n}|\right)^{p}}, \qquad |\mathbf{n}| = L\,|\mathbf{k}|/2\pi ,")
d.p("where **n** is the dimensionless FFT index, *L* the box length along each axis and **k** the wavenumber. The "
    "roll-off wavenumber is therefore *k*_{0} = 2π/*L* (outer scale equal to the box, 0.6–0.75 km), and the power "
    "spectrum decays as |**k**|^{-2p} above it; *p* = 1.5 gives one-dimensional spectra ∝ *k*^{-1.05} and *p* = 1.8 "
    "gives ∝ *k*^{-1.48} along lines through the field. Because the 1-D spectrum of an isotropic 3-D field with "
    "spectrum ∝ |**k**|^{-2p} decays approximately as *k*^{-(2p-2)}, the sonic-log slopes (β ≈ 0.85 over 12–200 m "
    "and 1.2–1.7 over 1–200 m; Section 2.3) correspond to *p* ≈ 1.4–1.85, the range of the simulated media. "
    "Exponential and Gaussian fields use the filters "
    "(1 + *k*²*a*²)^{-1} and exp(−*k*²*a*²/8) with the correlation length *a* in metres. All fields are periodic "
    "over the FD box, which is larger than any source–receiver path. After scaling, ln *V* is truncated at ±3σ "
    "(±2.5σ for σ > 0.15) without renormalization. Truncation affects "
    f"{FS['pct_cells_at_truncation'].min():.2f}–{FS['pct_cells_at_truncation'].max():.1f} % of cells and reduces σ by "
    f"≤ {100*max(1-FS.loc['pink s.13 p1.5','sigma_lnV_after_truncation']/0.13, 1-FS.loc['pink s.18 p1.5','sigma_lnV_after_truncation']/0.18):.1f} % "
    f"(Supporting Information, Table S9). The truncated velocity range is wide: at σ = 0.13, *V*_{{P}} spans "
    f"{FS.loc['pink s.13 p1.5','vp_min']/1e3:.1f}–{FS.loc['pink s.13 p1.5','vp_max']/1e3:.1f} km/s, and at σ = 0.18 "
    f"{FS.loc['pink s.18 p1.5','vp_min']/1e3:.1f}–{FS.loc['pink s.18 p1.5','vp_max']/1e3:.1f} km/s. The extremes "
    "occupy a small fraction of the volume, but they show that σ ≥ 0.13 is not a small perturbation.")
d.p("In the default media *V*_{P} and *V*_{S} share one field, so *V*_{P}/*V*_{S} is constant. Two variants "
    "relax this at σ = 0.13: independent fields for *V*_{P} and *V*_{S} (VSD), and a *V*_{S}-dominated medium with "
    "σ(ln *V*_{S}) = 0.13 and σ(ln *V*_{P}) = 0.065 (VSO). Because σ absorbs any unmodelled contribution to the "
    "late energy, including density contrasts, anisotropic fracture scattering, intrinsic attenuation and site "
    "effects, we call the fitted value an *effective scattering-strength parameter* σ_{eff} within this isotropic "
    "model class, not a measured property of the rock.")
d.p("Media were simulated in two sets (Table 2). A first set explores σ and *p* for the pink family. A second set "
    "forms a common grid for all three families: σ = 0.09, 0.13 and 0.18; *p* = 1.5 and 1.8 (pink), "
    "*a* = 15 and 50 m (exponential, Gaussian); two white-noise seeds each (11 and 33). Every family thus has the "
    "same number of candidates on the common grid. Fields of different families generated from the same seed "
    "share their phases and are strongly correlated (the pink and *a* = 15 m exponential fields with seed 33 by "
    f"{rexp*100:.0f} %); this is a controlled comparison of spectral shape, not a test with independent fields.")
rows = [["H", "smooth: homogeneous", "*V*_{P} 5808 m/s, *V*_{P}/*V*_{S} 1.73"],
        ["L", "smooth: 1-D layered", "56-32 sonic log, 50-m blocks"],
        ["GRAD3", "smooth: 3-D gradients", "±3 % vertical, ±2 % E and N"],
        ["A", "smooth: homogeneous VTI", "travel-time inverted ε, δ"],
        ["L6", "deterministic: 1-D layered", "56-32 sonic log at 6-m resolution"],
        ["FZ", "deterministic: fracture zones", "12 planar zones, 5 m thick, N25°E/75°, *V*_{P} −10 %, *V*_{S} −20 %"],
        ["P4, P1, P3, P6, P9", "pink, *p* = 1.5 (seed 11)", "σ = 0.02, 0.045, 0.09, 0.13, 0.18"],
        ["T1, R1, R2; R3; R6; R5", "pink, *p* = 1.5", "further seeds at σ = 0.045; 0.09; 0.13; 0.18"],
        ["P5, P2; P7, P8; P11", "pink, other *p*", "*p* = 1.2, 1.8 at σ = 0.045, 0.09; *p* = 1.8 at 0.18"],
        ["P12, P10", "pink, *p* = 1.5", "σ = 0.22, 0.25 (Supporting Information only)"],
        ["Common grid", "pink, exponential, Gaussian", "σ = 0.09, 0.13, 0.18 × (*p* = 1.5, 1.8 or *a* = 15, 50 m) × seeds 11, 33"],
        ["VSD, VSO", "pink, variable *V*_{P}/*V*_{S}", "σ = 0.13; independent *V*_{S} field; σ_{VS} = 0.13, σ_{VP} = 0.065"],
        ["R6c, R6d, E15c/d, G50c/d", "synthetic truths", "σ = 0.13, seeds 44 and 55 (shape-resolution test)"],
        ["T2", "synthetic truth", "50-m layering + VTI (ε = 0.06, δ = 0.03, γ = 0.06)"]]
d.table("Media simulated", ["Model", "Class", "Parameters"], rows, [4.4, 4.2, 7.4],
        legend="All media: ρ = 2650 kg/m³. Stochastic media: lognormal, truncated at ±3σ (±2.5σ for σ > 0.15). "
               "The common grid comprises 36 media (12 per family), 14 of which also belong to the first set.")
d.h("Source and preprocessing", 2)
fcA = sA.fc_apparent_Hz; fcB = sB.fc_apparent_Hz
d.p(f"FD synthetics carry a 45-Hz Ricker source. To give data and synthetics the same source spectrum, each observed "
    f"trace is convolved with the same Ricker wavelet. This is valid if the true source spectrum is flat in the "
    f"analysis band, i.e. if corner frequencies exceed ~80 Hz. The S-window displacement spectra (raw data, "
    f"5–300 Hz) support this for most but not all events. The apparent corner, the first frequency above 50 Hz "
    f"at which the spectrum falls 6 dB below its 30–50 Hz level, has a median of {fcA.median():.0f} Hz at "
    f"sensor A and {fcB.median():.0f} Hz at B. Because it includes path attenuation, it is a lower bound on the "
    f"source corner; {100*(fcA < 80).mean():.0f} % (A) and {100*(fcB < 80).mean():.0f} % (B) of events have "
    f"apparent corners below 80 Hz. To bound the effect of finite source duration, we convolved the homogeneous "
    f"and σ = 0.13 pink synthetics with {{@brune1970}} source-time functions with corners of 150, 100 and 60 Hz. Even the "
    f"60-Hz source raises the homogeneous coda-to-S ratio by only {sd('H',60,'A')-sd('H',np.inf,'A'):.2f} log units "
    f"(A) and changes the pink one by {sd('P6',60,'A')-sd('P6',np.inf,'A'):+.2f}. Source duration therefore cannot "
    f"close a deficit of 1–2 log units.")
d.h("Per-event coda prediction test", 2)
d.p("For each event and medium, the moment tensor is fitted by linear least squares to the four horizontal traces "
    "(two sensors, 40–80 Hz) in the window P − 10 ms to S + 15 ms, with a per-sensor time shift of up to ±6 ms and "
    "the tool azimuth and handedness of each sensor fixed from a previous grid search. The fitted moment tensor "
    "then predicts the whole record. After a guard band of 25 ms, the test quantity is the ratio of mean "
    "horizontal energy in S + 40 to S + 100 ms to the peak S energy (S − 8 to S + 15 ms), compared per event "
    "between prediction and observation as log_{10}(predicted/observed). The fit and coda windows do not overlap. "
    "The sensitivity analysis varies the end of the fit window (S + 10, 15, "
    "20 and, for comparison, an overlapping 60 ms) and the coda window (S + 30–100, 40–100 and 50–110 ms). It also "
    "compares a full moment tensor, a deviatoric one and a double couple (grid search over strike, dip and rake, "
    "with a least-squares scalar moment). Condition numbers of the column-normalized design matrix and the ISO, "
    "CLVD and DC fractions of the full solutions are reported. A timing check compares the observed S-envelope "
    "peak with the predicted S time. The test is repeated on the events whose S peak lies within ±10 ms (and "
    "±5 ms) of the prediction on both sensors.")
d.h("Tests for instrument and borehole ringing", 2)
d.p("A resonant sensor or borehole response *h*(*t*) moves energy from the S peak into later windows, so a coda-to-S "
    "ratio does not cancel it. We applied three tests. First, the S-window spectra and the coda/S spectral ratio "
    "were examined for resonance peaks. In the analysis band (30–120 Hz), any peak in the smoothed S spectrum is at "
    f"most {RP[RP.band_Hz=='30-120'].S_spectrum_peak_prominence_dB.max():.2f} dB above its surroundings. The same "
    f"measure detects the resonances near 330, 480 and 920 Hz (prominence up to {RP.S_spectrum_peak_prominence_dB.max():.0f} dB), "
    "which serve as a positive control. Second, in a null-ringing test the homogeneous synthetics were convolved "
    "with the strongest single resonance (40–80 Hz, *Q* = 5–40) whose spectral peak does not exceed the observed "
    "in-band prominence, and the coda-to-S ratio was recomputed. Third, the two sensors, which have different "
    "coupling and different in-band peaks, are analysed separately throughout. No calibration shots, perforation "
    "shots or third sensor were available for 2022, so an empirical transfer function could not be measured.")
d.h("Statistical comparison of media", 2)
d.p("With two sensors the realization of a random medium cannot be recovered, so stochastic media are compared with "
    "the data through wavefield statistics. On the horizontal energy envelope we measure three S-normalized "
    "features per event: the coda level log_{10}(*E*_{coda}/*E*_{S}) 30–100 ms after the S peak, the coda decay "
    "rate and the rms S-pulse width. We use three sensor–band combinations: A 20–40 Hz, A 40–80 Hz and B 40–80 Hz. "
    "B 20–40 Hz is excluded because it carries a persistent non-event plateau of unknown, probably mechanical or "
    "borehole, origin; the 40–80 Hz band used elsewhere is free of it. The misfit is")
d.eq(r"\Phi = \Phi_{\mathrm{KS}} + \Phi_{E} ,")
d.p("with", align="left")
d.eq(r"\Phi_{\mathrm{KS}} = \frac{1}{9}\sum_{k=1}^{9} D_{\mathrm{KS},k} ,")
d.eq(r"\Phi_{E} = \frac{1}{3}\sum_{c=1}^{3}\left[\frac{1}{N_{\tau}}\sum_{\tau}\left(\log_{10}\tilde{E}^{\mathrm{obs}}_{c}(\tau)-\log_{10}\tilde{E}^{\mathrm{mod}}_{c}(\tau)\right)^{2}\right]^{1/2} ,")
d.p("where *D*_{KS,k} is the two-sample Kolmogorov–Smirnov distance between the observed and synthetic "
    "distributions for the *k*-th of the nine pairs of feature (three) and sensor–band combination *c* (three), *Ẽ*_{c} is the median post-peak S envelope, and τ runs "
    "over the *N*_{τ} samples 0–100 ms after the S peak. All nine KS "
    "terms and all three envelope terms are weighted equally, events are weighted equally, and events lacking a "
    "measurement in one combination are omitted from that combination only. The features are correlated, so Φ is a "
    "ranking statistic, not a likelihood. Its uncertainty is estimated by an event bootstrap (200 resamples, the "
    "same resample for every medium). Family comparisons use the best member of each family on the common grid per "
    "bootstrap sample, and a holdout test in which each family's best (σ, shape, seed) is chosen on half of the "
    "events and scored on the other half (20 random splits plus a near/far split by distance).")
d.h("Inter-event coda coherence", 2)
d.p("For pairs of events recorded on the same sensor, D2 is the maximum normalized two-component cross-correlation "
    "(lag ±8 ms) of a 40–80 Hz coda window, summarized as the median in bins of inter-event separation (0–5 to "
    "80–160 m) and as the median over all pairs within 160 m. Three windows are used: S + 15 to S + 100 ms (early "
    "coda, including the tail of the S pulse), S + 40 to S + 100 ms (the window of the late-energy test) and "
    "S + 50 to S + 110 ms. The coherence of noise alone is measured on the same records in a window S + 550 to "
    "S + 610 ms. Every synthetic record contains that event's own recorded noise at its observed SNR, so the noise "
    "contribution is common to data and media. Both events of a pair are recorded by the same sensor in the same "
    "(unknown) orientation, so the two-component correlation is a scalar product of horizontal vectors and does not "
    "depend on the sensor azimuth. Pairs that share an event are not independent, so 95 % intervals "
    "come from a cluster bootstrap over events (each resampled event carries all its pairs). To test sensitivity "
    "to location error, the hypocentres used for binning were perturbed by isotropic Gaussian errors of 5, 10, "
    "20, 30 and 50 m, for data and media alike. The all-pairs median does not depend on the binning and is "
    "therefore insensitive to location error.")
d.h("Synthetic truth tests", 2)
d.p("The ability of the statistical comparison to select a medium class was tested with synthetic truths. Two "
    "truths were used in the main test: a pink medium (T1, σ = 0.045) and 50-m layering with VTI (T2). A further "
    "nine truths at equal σ = 0.13 (pink, exponential *a* = 15 m and Gaussian *a* = 50 m, seeds 33, 44 and 55) "
    "tested resolution of the spectral family. Synthetic records use double-couple mechanisms consistent with "
    "the FORGE stress field, and real 56-32 noise at each event's observed SNR.")

# ================================================================= RESULTS
d.h("Results", 1)
d.h("Smooth media under-predict the late S energy", 2)
c = lambda m: f"{cpa(m,'med_log10_pred_obs_A'):+.2f}/{cpa(m,'med_log10_pred_obs_B'):+.2f}"
d.p(f"Fig. 3 compares observed and predicted sensor-A seismograms for three events. The moment tensor fitted to "
    f"the direct waves reproduces the P and S pulses in all media. In the homogeneous medium the prediction then "
    f"decays to near zero within the guard band, whereas the observed traces keep oscillating through the coda "
    f"window. Table 3 and Fig. 4 give the statistics over {n_ev} events. The four smooth media under-predict the "
    f"coda-to-S ratio by a median factor of {fac(smA.max()):.0f}–{fac(smA.min()):.0f} at sensor A and "
    f"{fac(smB.max()):.0f}–{fac(smB.min()):.0f} at sensor B, and only {smwA.min():.0f}–{smwA.max():.0f} % (A) and "
    f"{smwB.min():.0f}–{smwB.max():.0f} % (B) of events fall within a factor of three. Smooth 3-D gradients "
    f"({c('GRAD3')}, A/B) and VTI anisotropy ({c('A')}) alter the direct waves but add no late energy.")
d.p(f"The two simple deterministic small-scale media tested change this substantially. The sonic log at its 6-m resolution (L6) "
    f"gives {c('L6')}, and the thin fracture zones (FZ) give {c('FZ')}: under-prediction by factors of "
    f"{fac(cpa('L6','med_log10_pred_obs_A')):.0f} and {fac(cpa('FZ','med_log10_pred_obs_A')):.0f} at sensor A but "
    f"only {fac(cpa('L6','med_log10_pred_obs_B')):.1f} and {fac(cpa('FZ','med_log10_pred_obs_B')):.1f} at sensor B, "
    f"where {cpa('L6','pct_within3_B'):.0f} % and {cpa('FZ','pct_within3_B'):.0f} % of events fall within a factor "
    f"of three. Averaging the log in 50-m blocks therefore removes most of the late "
    f"energy that the logged layering can generate. For the level of late energy, structure at the 5–10 m scale is "
    f"the essential ingredient, whether it is layered, planar or random. Stochastic media span the observations: "
    f"pink media with σ = 0.02 under-predict ({c('P4')}), σ = 0.09 is close at sensor A "
    f"({cpa('P3','med_log10_pred_obs_A'):+.2f}) and σ = 0.045 at sensor B ({cpa('P1','med_log10_pred_obs_B'):+.2f}), "
    f"and the seed-11 stochastic media with σ = 0.13 over-predict "
    f"({min(cpa(m,'med_log10_pred_obs_A') for m in ('P6','E15','E50','G15','G50')):+.2f} to "
    f"{max(cpa(m,'med_log10_pred_obs_B') for m in ('P6','E15','E50','G15','G50')):+.2f}). The best stochastic "
    f"media bring {max(cpa(m,'pct_within3_A') for m in STOCH):.0f} % (A) and "
    f"{max(cpa(m,'pct_within3_B') for m in STOCH):.0f} % (B) of events within a factor of three. The residual "
    f"scatter reflects the unknown realization. No single medium matches the median at both sensors; for the pink "
    f"family, sensor B requires about half the σ of sensor A.")
d.figure("Observed and predicted sensor-A seismograms with non-overlapping fit and coda windows",
         "Horizontal component 1, 40–80 Hz, for three events at increasing distance (rows; distance on the left). "
         "Columns: homogeneous medium, deterministic fracture zones (FZ), pink media with σ = 0.045 and 0.09. Black: "
         "observed; red: predicted with the moment tensor fitted in the light-blue window (P − 10 ms to S + 15 ms); "
         "grey: coda window (S + 40 to S + 100 ms); dotted: predicted S time. Variance reduction of the fit in each "
         "panel.", F / "fig_seis_rev.png")
rows = []
order = [("H", "smooth: homogeneous"), ("L", "smooth: layered, 50 m"), ("GRAD3", "smooth: 3-D gradients"),
         ("A", "smooth: VTI"), ("L6", "deterministic: layered, 6 m"), ("FZ", "deterministic: fracture zones"),
         ("P4", "pink σ = 0.02"), ("P1", "pink σ = 0.045"),
         ("P3", "pink σ = 0.09"), ("P6", "pink σ = 0.13"), ("P9", "pink σ = 0.18"), ("VSD", "pink σ = 0.13, independent *V*_{S}"),
         ("VSO", "pink, *V*_{S}-dominated"), ("X_exp_0.09_15_11", "exponential *a* = 15 m, σ = 0.09"),
         ("E15", "exponential *a* = 15 m, σ = 0.13"), ("X_exp_0.09_50_11", "exponential *a* = 50 m, σ = 0.09"),
         ("E50", "exponential *a* = 50 m, σ = 0.13"), ("X_gau_0.09_15_11", "Gaussian *a* = 15 m, σ = 0.09"),
         ("G15", "Gaussian *a* = 15 m, σ = 0.13"), ("X_gau_0.09_50_11", "Gaussian *a* = 50 m, σ = 0.09"),
         ("G50", "Gaussian *a* = 50 m, σ = 0.13")]
for m, name in order:
    if has(m):
        rows.append([name, f"{cpa(m,'med_log10_pred_obs_A'):+.2f}", f"{cpa(m,'med_log10_pred_obs_B'):+.2f}",
                     f"{cpa(m,'pct_within3_A'):.0f}", f"{cpa(m,'pct_within3_B'):.0f}", f2(cpa(m, 'VR'))])
d.table(f"Predicted versus observed late S energy for {n_ev} FORGE events",
        ["Medium", "Median log_{10}(pred/obs), A", "Median, B", "% within ×3, A", "% within ×3, B", "VR"], rows,
        [4.8, 3.4, 1.8, 2.0, 2.0, 1.2],
        legend="Late energy = mean horizontal energy S + 40 to S + 100 ms divided by the S-peak energy, 40–80 Hz. Moment "
               "tensors fitted in P − 10 ms to S + 15 ms (no overlap with the coda window). VR = median variance "
               "reduction of the direct-wave fit. Stochastic media: seed 11; other seeds in the Supporting Information.")
d.figure("Per-event predicted/observed late S energy for all media",
         f"a) Sensor A. b) Sensor B. {n_ev} events, 40–80 Hz, windows as in Table 3. Boxes: interquartile range; "
         "whiskers: 10–90 %; diamonds: median; dotted lines: factor of three. Grey: smooth media; purple: "
         "deterministic small-scale media (6-m layering, fracture zones); orange: pink media (light orange: "
         "variable *V*_{P}/*V*_{S}); blue: exponential; green: Gaussian.",
         F / "fig_coda_ratio_rev.png")
d.h("Robustness of the deficit", 2)
wsA = WS[(WS.model.isin(SMOOTH)) & (WS.overlap == 0)]
rows = []
def rng(vals):
    return f"{min(vals):+.2f} to {max(vals):+.2f}"
def four(fA, fB):
    return [rng([fA(m) for m in SMOOTH]), rng([fB(m) for m in SMOOTH]), rng([fA(m) for m in DET]), rng([fB(m) for m in DET])]


rows.append(["Primary (fit to S + 15 ms; coda S + 40–100 ms)"] + four(lambda m: cpa(m, 'med_log10_pred_obs_A'),
                                                                     lambda m: cpa(m, 'med_log10_pred_obs_B')))
for fe in (10, 20):
    rows.append([f"Fit to S + {fe} ms"] + four(lambda m: ws(m, fe, '40-100', 'med_A'), lambda m: ws(m, fe, '40-100', 'med_B')))
for cw in ("50-110", "30-100"):
    rows.append([f"Coda S + {cw.replace('-', '–')} ms"] + four(lambda m: ws(m, 15, cw, 'med_A'), lambda m: ws(m, 15, cw, 'med_B')))
rows.append(["Overlapping fit (S + 60 ms), coda S + 30–100 ms"]
            + four(lambda m: ws(m, 60, '30-100', 'med_A'), lambda m: ws(m, 60, '30-100', 'med_B')))
for v, lab_ in (("dev", "Deviatoric moment tensor"), ("dc", "Double couple")):
    rows.append([lab_] + four(lambda m: mv(m, v, 'med_A'), lambda m: mv(m, v, 'med_B')))
rows.append(["150 highest-SNR events"] + four(lambda m: CPt.loc[m, 'med_log10_pred_obs_A'],
                                              lambda m: CPt.loc[m, 'med_log10_pred_obs_B']))
n10 = int(stt('H', 'within10ms', 'n')); n5 = int(stt('H', 'within5ms', 'n'))
for s_, lab_ in (("within10ms", f"S peak within ±10 ms on both sensors ({n10} events)"),
                 ("within5ms", f"S peak within ±5 ms on both sensors ({n5} events)")):
    rows.append([lab_] + four(lambda m: stt(m, s_, 'med_log10_pred_obs_A'), lambda m: stt(m, s_, 'med_log10_pred_obs_B')))
rows.append(["Homogeneous + worst in-band resonance (null-ringing test)", rng(NR.deficit_A), rng(NR.deficit_B), "", ""])
rows.append(["Homogeneous, Brune source *f*_{c} = 60 Hz (change)", f"{sd('H',60,'A')-sd('H',np.inf,'A'):+.2f}",
             f"{sd('H',60,'B')-sd('H',np.inf,'B'):+.2f}", "", ""])
d.table("Robustness of the late-energy deficit of smooth and deterministic media",
        ["Variant", "Smooth, A", "Smooth, B", "L6 and FZ, A", "L6 and FZ, B"], rows, [7.0, 2.3, 2.3, 2.3, 2.3],
        legend=f"Median log_{{10}}(predicted/observed late energy); range over the smooth media ({', '.join(SMOOTH)}) and "
               f"over the deterministic small-scale media ({', '.join(DET)}). Null-ringing test: range over "
               "resonance frequencies 40–80 Hz and *Q* = 5–40, for the median over events of the homogeneous "
               "double-couple synthetics relative to the observed median. Brune row: change of the homogeneous "
               "coda-to-S ratio relative to an impulsive source.")
ov = max(abs(ws(m, 60, '40-100', c) - ws(m, 15, '40-100', c)) for m in SMOOTH for c in ('med_A', 'med_B'))
d.p(f"Table 4 shows that the deficit does not depend on the analysis choices. Ending the fit window at S + 10 or "
    f"S + 20 ms, or moving the coda window, changes the smooth-model medians by a few tenths of a log unit. Even a "
    f"fit window that overlaps the coda window (to S + 60 ms) changes them by at most {ov:.2f}, and a coda window "
    f"beginning at S + 30 ms gives a smaller deficit only because it includes more of the early coda. "
    f"Deviatoric and double-couple solutions give "
    f"the same result, although they fit the direct waves less well (median VR in the homogeneous medium "
    f"{f2(mv('H','full','VR'))} for the full, {f2(mv('H','dev','VR'))} for the deviatoric and "
    f"{f2(mv('H','dc','VR'))} for the double-couple solution). The moment-tensor problem "
    f"is moderately conditioned (median condition number {min(mv(m,'full','cond_median') for m in SMOOTH):.0f}–"
    f"{max(mv(m,'full','cond_median') for m in SMOOTH):.0f} for smooth and "
    f"{min(mv(m,'full','cond_median') for m in STOCH if m in set(MV.model)):.0f}–"
    f"{max(mv(m,'full','cond_median') for m in STOCH if m in set(MV.model)):.0f} for stochastic media, decreasing with σ). The full "
    f"solutions carry large non-double-couple parts (median |ISO| {mv('H','full','pct_ISO_abs'):.0f} %, "
    f"|CLVD| {mv('H','full','pct_CLVD_abs'):.0f} % in the homogeneous medium). With two sensors in one well these "
    f"parts are poorly resolved and we do not interpret them. They could absorb part of the path effects on the "
    f"direct waves, but the double-couple solutions cannot, and they give the same deficit. The lower condition "
    f"numbers in random media arise because scattered energy inside the fit window makes the six Green's-function "
    f"columns less similar. The point is that no source constraint changes the "
    f"coda deficit. The observed S peak lies within ±10 ms of the predicted time for "
    f"{n10} events on both sensors, and the deficit remains in that subset and in the ±5 ms subset. Timing and "
    f"location errors therefore do not explain it.")
d.p(f"The ringing tests constrain the instrument explanation (Fig. 5). Between 30 and 120 Hz the S spectra of both "
    f"sensors are smooth, and the coda/S spectral ratio has no peak. The worst resonance compatible with this "
    f"smoothness leaves the homogeneous medium {-NR.deficit_A.max():.1f}–{-NR.deficit_A.min():.1f} log units short "
    f"at sensor A and {-NR.deficit_B.max():.1f}–{-NR.deficit_B.min():.1f} at B. The resonances above 300 Hz are "
    f"strong and clearly visible, so an in-band resonance large enough to create the observed coda would not have "
    f"escaped detection. This test excludes a single resonant mode, but not a broadband, non-resonant coupling "
    f"response with a long impulse response. That possibility can be excluded only with an empirical transfer "
    f"function, which is not available.")
d.figure("Tests for instrument or borehole ringing",
         "a, b) Sensors A and B: median S-window spectrum normalized over 30–150 Hz (black) and median coda/S "
         "spectral ratio (red, normalized). Dotted lines mark the analysis band 40–80 Hz and grey lines the "
         "resonances near 330, 480 and 920 Hz. c) Null-ringing test: log_{10} coda/S of the homogeneous synthetics "
         "(sensor A) after convolution with the strongest resonance compatible with the observed in-band spectral "
         "smoothness, versus resonance frequency, for *Q* = 5–40. Thick line: observed FORGE median.",
         F / "fig_ringing_test.png")
if CV is not None:
    d.p(f"Grid refinement and sponge thickness change the late energy little (Supporting Information, Table S10). "
        f"Relative to the 5-m reference, the median coda-to-S ratio at 40–80 Hz changes by "
        f"{cvs('C25','40-80 Hz'):+.2f} log units on the 2.5-m grid and by {cvs('CSP','40-80 Hz'):+.2f} with the "
        f"doubled sponge ({cvs('C25','20-40 Hz'):+.2f} and {cvs('CSP','20-40 Hz'):+.2f} at 20–40 Hz). Enlarging the "
        f"domain by 150 m on every side, however, raises it by {cvs('CBIG','40-80 Hz'):+.2f} at 40–80 Hz "
        f"({cvs('CBIG','20-40 Hz'):+.2f} at 20–40 Hz) and lowers the coda decay rate from "
        f"{cvd('C5'):.0f} to {cvd('CBIG'):.0f} s^{{-1}}. Because the doubled sponge has a much smaller effect, this "
        f"is not a boundary artefact. It reflects the larger scattering volume: in the reference domain the event "
        f"cloud nearly fills the interior, so heterogeneity beyond it lies in the absorbing sponge and returns no "
        f"scattered energy. For lapse times up to 100 ms after S, scatterers contribute to the coda window only within "
        f"about *V*_{{S}}*t*/2 ≈ 170 m of the source–receiver paths, so a 150-m extension covers most of this volume and "
        f"the effect should be close to saturation; larger domains were not tested. A smooth medium has no such "
        f"scatterers, so the smooth-model deficit is unaffected. The "
        f"slower decay in the larger domain is closer to the observed decay. For the stochastic media, the reference "
        f"domain probably underestimates the late energy, and the effective scattering strength derived below is probably an upper "
        f"bound. With the sensitivity of the median ratio to σ found below ({rslope:.2f} log units per unit "
        f"change of ln σ), this single test implies that σ_{{eff}} from the reference domain may be overestimated by a "
        f"factor of about {fdom:.1f}. Whether the same factor applies to other families and to sensor A was not tested.")
else:
    d.p(pend("convergence results (S26)"))
d.h("Effective scattering strength", 2)
if SE is not None:
    d.p(f"For each family the σ at which the median predicted late energy equals the observed one (σ_{{eff}}) is "
        f"interpolated in ln σ between the simulated levels of the reference domain (Fig. 6; Supporting Information, "
        f"Table S11). At sensors A and B it is {sigeff('pink15','A')} and {sigeff('pink15','B')} for pink media with "
        f"*p* = 1.5, {sigeff('pink18','A')} and {sigeff('pink18','B')} with *p* = 1.8, {sigeff('exp15','A')} and "
        f"{sigeff('exp15','B')} for exponential media with *a* = 15 m, {sigeff('exp50','A')} and {sigeff('exp50','B')} "
        f"with *a* = 50 m, and {sigeff('gau50','A')} and {sigeff('gau50','B')} for Gaussian media with *a* = 50 m "
        f"(\"< 0.09\": the smallest simulated σ already over-predicts). σ_{{eff}} thus differs between sensors and "
        f"families, so it is a model-dependent summary of "
        f"the late energy rather than a unique rock property. Media with less small-scale content need larger σ "
        f"(Gaussian *a* = 50 m), and media richest in structure near the S wavelength need less (Gaussian *a* = 15 m: "
        f"{sigeff('gau15','A')} at both sensors). Realization scatter adds to this. The two seeds at one σ differ "
        f"in the median ratio by up to {SPR['pink']:.2f} log units for pink media and {SPR['exp15']:.2f} for "
        f"exponential media with *a* = 15 m, i.e. a factor of about {np.exp(rspread/rslope):.1f} in σ_{{eff}}, but by "
        f"up to {SPR['gau50']:.2f} for Gaussian media with *a* = 50 m, whose few large features make each "
        f"realization idiosyncratic in a 0.5-km box.")
else:
    d.p(pend("sigma_eff (S29)"))
if has("VSD") and has("VSO"):
    d.p(f"Relaxing the constant *V*_{{P}}/*V*_{{S}} ratio changes the late energy only moderately. With independent "
        f"*V*_{{P}} and *V*_{{S}} fields the median ratio is {cpa('VSD','med_log10_pred_obs_A'):+.2f}/"
        f"{cpa('VSD','med_log10_pred_obs_B'):+.2f} (A/B), against {cpa('P6','med_log10_pred_obs_A'):+.2f}/"
        f"{cpa('P6','med_log10_pred_obs_B'):+.2f} for the coupled medium. The *V*_{{S}}-dominated medium, with half "
        f"the *V*_{{P}} perturbation, gives {cpa('VSO','med_log10_pred_obs_A'):+.2f}/"
        f"{cpa('VSO','med_log10_pred_obs_B'):+.2f}. The S coda is controlled mainly by the *V*_{{S}} perturbation, "
        f"so σ_{{eff}} should be read as a *V*_{{S}} scattering strength.")
else:
    d.p(pend("VSD/VSO results"))
d.figure("Median predicted/observed late S energy versus σ for each spectral family",
         "a) Sensor A. b) Sensor B. Median over events of log_{10}(predicted/observed), averaged over realizations "
         "where several exist. The zero crossing defines σ_{eff}; dotted lines: factor of three; grey band: range of "
         "the smooth media; purple dash-dotted and dashed lines: 6-m layering (L6) and fracture zones (FZ).", F / "fig_sigma_eff.png")
d.h("Inter-event coda coherence", 2)
dpert = lambda w, g, sp: D2P[(D2P.model == "data") & (D2P.window == w) & (D2P.sensor == g) & (D2P.perturb_m == sp)
                              & (D2P.d_lo_m >= 20)].median_coh
dmax = max(abs(dpert(w, g, 50).values - D2W[(D2W.model == "data") & (D2W.window == w) & (D2W.sensor == g)
                                             & (D2W.d_lo_m >= 20)].median_coh.values).max()
           for w in (W15, W40, W50) for g in "AB")
d.p(f"Fig. 7 shows D2 in the three windows. In the early coda (S + 15 to S + 100 ms) the observed coherence falls "
    f"from {d2('data','A',0):.2f} at 0–5 m to {d2('data','A',80):.2f} at 80–160 m (sensor A; {d2('data','B',0):.2f} to "
    f"{d2('data','B',80):.2f} at B). The smooth media stay coherent (all pairs within 160 m: {dar(D2SM, W15, 'A')} "
    f"at A, {dar(D2SM, W15, 'B')} at B, against "
    + (f"{da('data', W15, 'A'):.2f} at both sensors for the data)" if f"{da('data', W15, 'A'):.2f}" == f"{da('data', W15, 'B'):.2f}"
       else f"{da('data', W15, 'A'):.2f} and {da('data', W15, 'B'):.2f} for the data)") + f", and so do the two deterministic small-scale media ({dar(D2DT, W15, 'A')} and {dar(D2DT, W15, 'B')}). "
    f"The stochastic media decorrelate as the data do ({dar(D2ST, W15, 'A')} and {dar(D2ST, W15, 'B')}). This "
    f"window, however, contains the tail of the S pulse and is close to the window of the source fit.")
d.p(f"In the late windows the picture changes. All coherences are lower and approach the coherence of noise "
    f"alone ({dnoise('A'):.2f} at A and {dnoise('B'):.2f} at B). The observed late coda is at or below this level "
    f"({da('data', W40, 'A'):.2f} and {da('data', W40, 'B'):.2f} for S + 40 to S + 100 ms; "
    f"{da('data', W50, 'A'):.2f} and {da('data', W50, 'B'):.2f} for S + 50 to S + 110 ms). Because the synthetics "
    f"carry the same noise, D2 in these windows measures how much coherent late energy a medium adds. The smooth "
    f"media ({dar(D2SM, W40, 'A')} at A, {dar(D2SM, W40, 'B')} at B) and the fracture zones "
    f"({da('FZ', W40, 'A'):.2f}, {da('FZ', W40, 'B'):.2f}) remain more coherent than the data in both late windows; "
    f"their 95 % intervals lie above the data interval in {nsep} of {ntot} medium–sensor–window cases "
    f"(Supporting Information, Table S17). The stochastic media "
    f"({dar(D2ST, W40, 'A')}, {dar(D2ST, W40, 'B')}) and the 6-m layering ({da('L6', W40, 'A'):.2f}, "
    f"{da('L6', W40, 'B'):.2f}) come closest; in S + 50 to S + 110 ms the 6-m layering is as close to the data as "
    f"the stochastic media at sensor A ({da('L6', W50, 'A'):.2f} against {dar(D2ST, W50, 'A')} and "
    f"{da('data', W50, 'A'):.2f} observed). In the pure late coda D2 therefore separates the data from the smooth "
    f"media and from the planar fracture zones, but not from fine 1-D layering. The separation of distributed 3-D "
    f"heterogeneity from fine layering rests on the early-coda window and is correspondingly less secure.")
d.p(f"Location error does not change these conclusions. Perturbing the hypocentres of data and media by up to "
    f"50 m flattens the binned curves, as expected, but changes the observed medians at separations of 20 m or more "
    f"by at most {dmax:.2f} (Supporting Information, Fig. S10), and the all-pairs levels, which do not depend on the "
    f"binning, are unaffected. We therefore do not use the 0–5 and 5–10 m bins for inference (the 0–5 m bin "
    f"contains only {int(d2('data','A',0,'n_pairs'))} pairs from {int(d2('data','A',0,'n_events'))} events). D2 is "
    f"measured on the same records as the late energy, so it is a second descriptor of the same data, not an "
    f"independent experiment.")
d.figure("Inter-event coda coherence in early and late coda windows",
         "Median coherence of the 40–80 Hz coda versus inter-event separation, with 95 % cluster-bootstrap "
         "intervals over events, for sensor A (a–c) and sensor B (d–f), in the windows S + 15 to S + 100 ms "
         "(a, d), S + 40 to S + 100 ms (b, e) and S + 50 to S + 110 ms (c, f). Black: FORGE; dotted: FORGE "
         "noise-only window (S + 550 to S + 610 ms); grey: smooth media (H, L, GRAD3, A); purple: 6-m layering and "
         "fracture zones; orange: pink (σ = 0.09, 0.13); blue: exponential (*a* = 15 m); green: Gaussian "
         "(*a* = 50 m), σ = 0.13. Symbols are offset horizontally for clarity.", F / "fig_d2_windows.png")
d.h("Statistical comparison and the spectral family", 2)
d.p(f"In the synthetic tests the statistical comparison selected the generating class in both controlled examples. "
    f"With the pink truth T1 the best medium was a pink medium with the true σ "
    f"({M1.PHI.idxmin()}, Φ = {f2(M1.PHI.min())}), ahead of H, L and A "
    f"(Φ = {f2(M1.loc[['H','L','A']].PHI.min())}"
    + ("" if f2(M1.loc[['H','L','A']].PHI.min()) == f2(M1.loc[['H','L','A']].PHI.max())
       else f"–{f2(M1.loc[['H','L','A']].PHI.max())}") + " each). With the layered-VTI truth T2, "
    f"L and A ranked first ({f2(M2.loc['L','PHI'])}, {f2(M2.loc['A','PHI'])}) and every pink medium was rejected. "
    f"Two examples do not establish that the procedure is unbiased. The deterministic moment-tensor fit behaved "
    f"differently. For the pink truth it preferred the smoothest media (median VR {f2(MT1.loc['H','VR_med'])} for H "
    f"against {f2(MT1.loc['P1','VR_med'])} for a pink medium with the correct statistics but another realization), "
    f"because a wrong realization places scattered energy at wrong times. This applies to selecting a model class "
    f"from a set of arbitrary stochastic realizations. It does not imply that full-waveform inversion, which "
    f"updates the structure itself, cannot recover heterogeneity when the data coverage allows it.")
d.p(f"The equal-σ family test selected the correct family for {n_ok} of 9 truths. Under a chance model with "
    f"*p* = 1/3 the probability of at least {n_ok} successes is {p_binom:.3f}, so this is not significant. Pink and "
    f"exponential truths were assigned to the correct member of that pair in {n_pe} of 6 cases (chance level), "
    f"whereas small-scale-rich (pink or exponential) and Gaussian truths were separated in {n_pl} of 9. With so few "
    f"realizations, even the latter is only indicative.")
FAMLAB = {"pink": "pink", "exp15": "exponential *a* = 15 m", "exp50": "exponential *a* = 50 m",
          "gau15": "Gaussian *a* = 15 m", "gau50": "Gaussian *a* = 50 m", "smooth": "smooth",
          "deterministic": "deterministic", "pink-VpVs": "pink, variable *V*_{P}/*V*_{S}"}


def seed_of(m):
    if m.startswith("X_"):
        return int(m.split("_")[-1])
    return 33 if (m.endswith("b") or m.startswith("R")) else 11


if PB is not None and FB is not None and HS is not None:
    fb = FB.set_index("family"); hs = HS.set_index("family")
    SF = ["pink", "exp15", "exp50", "gau15", "gau50"]
    PBf = PB[PB.family.isin(SF)].copy(); PBf["seed"] = PBf.model.map(seed_of)
    GP = (PBf.groupby(["family", "sigma", "shape"]).agg(PHI_mean=("PHI", "mean"), n=("PHI", "size"),
                                                         PHI_min=("PHI", "min"), PHI_max=("PHI", "max")).reset_index())
    GP = GP[GP.n >= 2]; GP["dseed"] = GP.PHI_max - GP.PHI_min
    BA = GP.loc[GP.groupby("family").PHI_mean.idxmin()].sort_values("PHI_mean")
    win = fb.pct_boot_wins.idxmax()
    ba_txt = "; ".join(f"{FAMLAB[r.family]} {r.PHI_mean:.2f} (σ = {r.sigma:.2f}"
                       + (f", *p* = {r.shape:g})" if r.family == "pink" else ")") for r in BA.itertuples())
    sm_min = min(fb.loc["smooth", "best_PHI"], fb.loc["deterministic", "best_PHI"])
    st_max = BA.PHI_mean.max()
    d.p(f"On the FORGE data (Fig. 8), the smooth media (best Φ = {fb.loc['smooth','best_PHI']:.2f}) and the "
        f"deterministic small-scale media ({fb.loc['deterministic','best_PHI']:.2f}) are separated from every "
        f"stochastic family (best members {FB[FB.family.isin(SF)].best_PHI.min():.2f}–"
        f"{FB[FB.family.isin(SF)].best_PHI.max():.2f}) by far more than the bootstrap intervals. Between the "
        f"stochastic families, the best single member belongs to the {FAMLAB[win]} family (Φ = "
        f"{fb.loc[win,'best_PHI']:.2f} [{fb.loc[win,'lo95']:.2f}, {fb.loc[win,'hi95']:.2f}]), which has the lowest "
        f"Φ in {fb.loc[win,'pct_boot_wins']:.0f} % of event-bootstrap resamples, ahead of the pink family "
        f"({fb.loc['pink','best_PHI']:.2f} [{fb.loc['pink','lo95']:.2f}, {fb.loc['pink','hi95']:.2f}]). In the "
        f"holdout test the best training member scores a median held-out Φ of {hs.loc[win,'median_PHI_test']:.2f} "
        f"({FAMLAB[win]}) and {hs.loc['pink','median_PHI_test']:.2f} (pink). Both tests, however, resample events "
        f"with the realization fixed, and the holdout selected the same realization in almost every split. They "
        f"show that the ranking generalizes across events, not across realizations. The realization matters more: "
        f"the two seeds at the same (σ, shape) differ by up to {GP.dseed.max():.2f} in Φ (median "
        f"{GP.dseed.median():.2f}), whereas the best members of the leading families differ by "
        f"{abs(fb.loc['pink','best_PHI'] - fb.loc[win,'best_PHI']):.2f}. Averaged over both seeds, the best grid point "
        f"of each family gives {ba_txt}. With two realizations per grid point, the 20–80 Hz statistics therefore do "
        f"not select a spectral family; Gaussian media with *a* = 15 m fit least well, but even they are separated "
        f"from the smooth and deterministic media by a margin ({sm_min - st_max:.2f} in Φ) larger than any "
        f"realization effect.")
else:
    d.p(pend("family bootstrap / holdout (S20)"))
d.figure("Statistical misfit Φ with event-bootstrap intervals and holdout scores",
         "a) Φ for all media with the full event set (circles) and 95 % event-bootstrap intervals (bars), sorted by "
         "Φ. Labels give family, σ, *p* or *a* (m) and white-noise seed (s11, s33). Grey: smooth media; purple: "
         "deterministic small-scale media; orange: pink; light orange: variable *V*_{P}/*V*_{S}; blue: exponential "
         "(dark *a* = 15 m, light *a* = 50 m); green: Gaussian (dark *a* = 15 m, light *a* = 50 m). The intervals "
         "reflect event resampling only, not realization scatter. b) Held-out Φ of each family's best training member "
         "over 21 splits in both directions (box: interquartile range; whiskers: range).", F / "fig_phi_family.png")
d.h("Sonic-log spectrum", 2)
d.p("Fig. 9 shows the log spectrum with four fitted spectral models. The spectrum is consistent with a power law "
    "over 12–200 m, but this range spans only 1.2 decades. Short-correlation-length spectra fit almost as well, "
    "and the estimated slope depends on estimator and segment. The log therefore does not select a spectral "
    "family independently either.")
d.figure("Spectral models of the 56-32 sonic-log fluctuations",
         "Multitaper power spectrum of detrended ln *V* (grey) and log-binned values (dots), with least-squares fits "
         "of a power law, exponential and Gaussian autocorrelation spectra and a power law with a corner (dashed). "
         "a) 12–200 m; the corner model converges to a pure power law and overlies it. b) 1–200 m. AIC values in the "
         "legend.",
         F / "fig_log_spectrum_uncertainty.png")

# ================================================================= DISCUSSION
d.h("Discussion", 1)
import scipy.io as sio
dom = cvs("CBIG", "40-80 Hz")
p1A = cpa("P1", "med_log10_pred_obs_A") + dom; p1B = cpa("P1", "med_log10_pred_obs_B") + dom
# intrinsic attenuation: energy loss of the coda window (mean delay 70 ms after the S peak) at 60 Hz
qloss = lambda Q: 2 * np.pi * 60 * 0.070 / Q / np.log(10)
qfac = np.exp(qloss(50) / rslope)
# distance dependence of observed and predicted late energy (S19 files)
S0 = sio.loadmat(str(R / "s19_H.mat")); R0 = S0["R"]; OB = np.log10(S0["obsratio"][:, :, 0])
def _slope(model, g):
    S = sio.loadmat(str(R / f"s19_{model}.mat"))
    pr = S["lr"][:, g, 0, 0] + np.log10(S["obsratio"][:, g, 0]); rr = S["R"][:, g]
    ok = np.isfinite(pr) & np.isfinite(rr)
    return np.polyfit(np.log10(rr[ok]), pr[ok], 1)[0]
def _drop(g):
    q = np.nanquantile(R0[:, g], [1 / 3, 2 / 3])
    return np.nanmedian(OB[R0[:, g] <= q[0], g]) - np.nanmedian(OB[R0[:, g] > q[1], g])
STM = ["P1", "P3", "P6", "E15", "G50"]
sens = [max(abs(_slope(m, g)) for m in STM) * np.log10(1.15) for g in (0, 1)]
sigf = np.exp(max(sens) / rslope)
RA, RB = np.nanmedian(R0[:, 0]), np.nanmedian(R0[:, 1])
obsA, obsB = np.nanmedian(OB[:, 0]), np.nanmedian(OB[:, 1])
allpred = [_slope(m, g) for m in STM + ["H"] for g in (0, 1)]
def _oslope(g):
    ok = np.isfinite(OB[:, g]) & np.isfinite(R0[:, g])
    return np.polyfit(np.log10(R0[ok, g]), OB[ok, g], 1)[0]

d.h("What the data show", 2)
d.p("Three results are robust. Media without structure below ~50 m produce one to two orders of magnitude less "
    "late S energy than recorded. Media with deterministic structure at the 5–10 m scale supply much of the missing "
    "energy, but not the decorrelation of the early coda between neighbouring events. Random media reproduce both, "
    "at a fluctuation strength of the same order as that of the local sonic log. None of this depends on the "
    "power-law (\"pink-noise\") hypothesis, which the data neither require nor exclude. The main alternative to "
    "small-scale scattering is an unmodelled site or coupling response of comparable size, which only a measured "
    "sensor transfer function or a third, independently installed sensor could rule out.")
d.h("Relation to earlier scattering studies", 2)
d.p("Two-dimensional FD studies showed that coda strength and apparent attenuation are controlled by the "
    "fluctuation strength and by the part of the heterogeneity spectrum near the wavelength [@frankel1986; "
    "@roth1993; @shapiro1993; @hestholm1994], and envelope theory for power-law media relates coda excitation to "
    "the short-wavelength part of the spectrum [@saito2003; @sato2012]. Our 3-D elastic simulations show the same "
    "behaviour in a borehole geometry at 40–80 Hz, where S wavelengths are 40–85 m. At fixed σ, the media richest in "
    "structure near the S wavelength (Gaussian, *a* = 15 m) produce the most late energy and need the smallest "
    "σ_{eff}, and the media whose heterogeneity is concentrated at large scales (Gaussian, *a* = 50 m) produce the "
    "least. This is why spectral family and fluctuation strength trade off. Three-dimensional simulations of "
    "regional high-frequency wavefields have examined how scattering and intrinsic attenuation together shape "
    "envelopes [@takemura2017]; our data constrain the two only jointly (Section 5.3). The required fluctuation "
    "strength is of the same order as the velocity fluctuations measured in sonic logs [@wu1994; @holliger1996], "
    f"including the 56-32 log ({K.log_sig6:.3f}).")
d.h("Intrinsic attenuation, density and the source", 2)
d.p(f"The FD media are elastic. Intrinsic attenuation reduces the energy of the coda window relative to the S peak "
    f"by exp(−2π*f*Δ*t*/*Q*_{{S}}), where Δ*t* ≈ 70 ms is the mean delay of the coda window. At *f* = 60 Hz this "
    f"amounts to {qloss(50):.2f}, {qloss(100):.2f} and {qloss(300):.2f} log units for *Q*_{{S}} = 50, 100 and 300. "
    f"Intrinsic attenuation therefore increases, rather than reduces, the smooth-media deficit, and no value of "
    f"*Q*_{{S}} explains the late energy without scattering. In random media it raises the σ_{{eff}} needed to match "
    f"the data by at most a factor of {qfac:.1f} for *Q*_{{S}} ≥ 50, comparable to the realization scatter. "
    f"Separating scattering from intrinsic attenuation, as done regionally from S-wave energy versus distance "
    f"[@fehler1992], would need more sensors and a wider distance range than two sensors in one well provide. "
    f"Density is constant in all media. If it co-varied with velocity as in Gardner's relation, ρ ∝ *V*_{{P}}^{{1/4}} "
    f"[@gardner1974], impedance fluctuations would be about 25 per cent larger than the velocity fluctuations, and a "
    f"correspondingly smaller σ_{{eff}} (by up to ~20 per cent) would produce the same scattering. Finally, the "
    f"moment tensor is fitted only to the direct waves. Its non-double-couple parts cannot create late energy that "
    f"the double-couple solutions lack (Table 4), and a finite source duration down to a 60-Hz Brune corner raises "
    f"the late energy of the homogeneous medium by only {sd('H',60,'A')-sd('H',np.inf,'A'):.2f} log units "
    f"(Section 3.3).")
d.h("What would resolve the spectral family", 2)
d.p("At 20–80 Hz the per-event test, the common-grid statistics and the equal-σ truth tests do not separate "
    "power-law from exponential or Gaussian spectra, and the sonic log does not separate them either (Section 4.6). "
    "The families differ mainly at wavenumbers well above 1/*a*; for *a* = 15 m this means S wavelengths below "
    "~15–30 m, i.e. frequencies above ~110–220 Hz. An exploratory comparison at 100 Hz on a 2.5-m grid suggested "
    "that a large-scale Gaussian medium loses its coda at that frequency (Supporting Information, Fig. S9), but it "
    "used a fixed σ, two realizations and noise-free synthetics, and the coda attenuation Q_{c}^{-1} is only a "
    "finite-window diagnostic at 25 Hz (Figs S7 and S8). Resolving the family would require (i) simulations to at "
    "least 200 Hz on grids of 2.5 m or finer, in domains that extend at least ~170 m beyond the event cloud; "
    "(ii) ensembles of 10 or more realizations per parameter set, because realization scatter now exceeds the "
    "family differences; (iii) sensors at several azimuths and depths, so that realization-specific features average "
    "out; and (iv) measured sensor transfer functions, because fixed-frequency resonances near 330, 480 and 920 Hz "
    "prevent the use of the band above ~200 Hz. The 4-kHz FORGE records carry S energy above noise up to ~1.8 kHz "
    "(Figs S5 and S6), so the bandwidth exists; the limits are instrument calibration and computational cost.")
d.h("The effective scattering strength and the two sensors", 2)
d.p(f"σ_{{eff}} is not a unique measure of rock heterogeneity. It absorbs every unmodelled contribution to the "
    f"late energy: density contrasts, anisotropic or fluid-filled fracture scattering, intrinsic attenuation, and "
    f"site or coupling effects that the ringing tests cannot exclude. It also depends on the simulation domain. One "
    f"domain test (one medium, sensor B) suggests that the reference-domain values ({se_txt}) may be about "
    f"{fdom:.1f} times too large. With that correction the pink medium with the sonic-log value (σ = 0.045) would "
    f"give median ratios of {p1A:+.2f} (A) and {p1B:+.2f} (B), within a factor of about 3 of the data at both "
    f"sensors. This is consistency, not a measurement.")
d.p(f"σ_{{eff}} at the shallower sensor A is about twice that at the deeper sensor B for the pink and exponential "
    f"families. A sensor-depth error cannot explain this. An error of 50 m changes hypocentral distances by about "
    f"15 per cent (median distances {RA:.0f} m at A and {RB:.0f} m at B), and the predicted late-energy ratio of the "
    f"random media varies with distance slowly enough that this shifts it by at most {sens[0]:.2f} (A) and "
    f"{sens[1]:.2f} (B) log units, a factor of at most {sigf:.1f} in σ_{{eff}}. Distance itself does not explain it "
    f"either. In the data the late energy decreases with distance at both sensors, by {_drop(0):.2f} (A) and "
    f"{_drop(1):.2f} (B) log units from the nearest to the farthest third of the events (log–log slopes of "
    f"{_oslope(0):.1f} and {_oslope(1):.1f}), whereas the simulated media predict slopes between {min(allpred):.1f} "
    f"and {max(allpred):+.1f}, i.e. a much weaker decrease or an increase. Moreover, sensor A records the events at larger "
    f"distances but has the stronger late energy (median log_{{10}} coda/S {obsA:.2f} against {obsB:.2f}). The "
    f"statistically uniform random media therefore reproduce neither the distance dependence nor the sensor "
    f"dependence of the late energy. Both point to spatially non-uniform scattering, for example stronger "
    f"heterogeneity along the paths to sensor A, which cross the upper part of the stimulated volume, or late energy "
    f"generated close to the sources, or to a sensor-specific site response. Two sensors in one well cannot separate "
    f"these possibilities.")
d.p("Stimulation might raise the scattering strength by opening fluid-filled fractures, which lower *V*_{S} "
    "preferentially. This is a hypothesis. No pre-stimulation waveforms exist for these paths, and a stage "
    "comparison is not possible, because 402 of the 416 events belong to stage 3 (Supporting Information, Table S8). "
    "The *V*_{S}-dominated medium shows only that a *V*_{S}-weighted perturbation produces similar late energy, and "
    "the fracture-zone medium that twelve planar zones are not enough. Denser 3-D fracture networks, rough "
    "interfaces, irregular damage zones and anisotropic compliance were not tested, and a sufficiently dense "
    "deterministic network would approach a random medium.")
d.h("Implications for induced-seismicity monitoring", 2)
d.p("Velocity models adequate for locating events at ~1-ms travel-time precision leave most of the late wavefield "
    "unexplained. This matters for EGS monitoring and induced-seismicity assessment [@majer2007; @ellsworth2013] in "
    "three ways. First, waveform-based source studies that assume smooth media, such as moment-tensor inversion, "
    "source spectra and coda-based magnitudes, attribute the late energy to noise or to the source; with few sensors, "
    "non-double-couple components and apparent corner frequencies are particularly exposed. Second, high-frequency "
    "ground-motion and attenuation estimates for induced events depend on the balance between scattering and "
    "intrinsic attenuation [@fehler1992; @takemura2017], which these data constrain only jointly. Third, the coda "
    "carries information on the small-scale structure of the stimulated volume, and its decorrelation between "
    "events and over time can be followed with coda-wave methods [@snieder2002; @pradhan2025]. Such comparisons "
    "must be statistical, because a deterministic comparison with an arbitrary random realization favours smooth "
    "media (Section 4.5). For models of crustal heterogeneity, the results are consistent with the power-law "
    "fluctuations inferred from well logs [@leary1991; @wu1994; @holliger1996], but they do not discriminate them "
    "from media with a short correlation length or, in the late coda, from fine layering.")
d.h("Limitations", 2)
d.bullets([
    "Two sensors in one well, with unoriented horizontals and unreliable vertical channels; no calibration shots, "
    "perforation shots or independent sensors for 2022, so the sensor–borehole transfer function is not measured.",
    "No pre-stimulation records exist for these paths, so changes caused by the stimulation cannot be assessed.",
    "Sensor positions, clock offsets and *V*_{P} trade off with the catalogue locations (~50 m, ~8 per cent).",
    "The random media are isotropic, lognormal and elastic, with constant density; intrinsic attenuation, "
    "anisotropic fracture compliance and rough interfaces are not modelled.",
    "The domain-size effect on the random-media coda was measured for one medium and one sensor; the correction of "
    "σ_{eff} is approximate.",
    f"The common grid has two realizations per grid point; realization scatter limits σ_{{eff}} to about a factor of "
    f"{np.exp(rspread/rslope):.1f} for pink and short-*a* media and much more for large-*a* Gaussian media.",
    "The analysis band is 20–80 Hz, where the spectral families differ least, and in the pure late coda the "
    "inter-event coherence is limited by noise.",
])

# ================================================================= CONCLUSIONS
d.h("Conclusions", 1)
d.p("(1) In 416 downhole records of the 2022 FORGE stimulation, the S-wave energy 40–100 ms after the S arrival "
    "exceeds that of smooth elastic media by one to two orders of magnitude. The excess is not an artefact of the "
    "analysis windows, the source parametrization, event timing or a resonant sensor response.")
d.p("(2) A model that reproduces this energy needs structure at the 5–10 m scale. Fine layering and thin planar "
    "zones provide much of the energy; the early-coda decorrelation between neighbouring events is reproduced only "
    "by distributed 3-D heterogeneity, whereas in the noise-limited late coda only smooth and planar-zone media are "
    "excluded.")
d.p(f"(3) In isotropic lognormal media the required fluctuation strength is {se_txt}. It depends on sensor and "
    f"spectral family, is probably an upper bound, and is of the same order as the 56-32 sonic-log value.")
d.p("(4) Neither the 20–80 Hz waveforms nor the sonic log resolve the spectral family. Doing so requires "
    "frequencies above ~200 Hz, calibrated sensors, several stations and realization ensembles.")
d.p("(5) Selecting a medium class by matching waveforms to one random realization is biased towards smooth media; "
    "statistical descriptors of the wavefield avoid this bias.")

# ================================================================= BACK MATTER (GJI)
d.h("Acknowledgements", 1, numbered=False)
d.p("We thank Heiner Igel for the acoustic finite-difference MATLAB codes from which the elastic solver used here "
    "was developed, Geo-Energie Suisse and the Utah FORGE project of the U.S. Department of Energy for the event "
    "catalogue and the downhole records, and an internal reviewer for critical readings of earlier versions. This "
    "research received no specific grant from any funding agency in the public, commercial or not-for-profit "
    "sectors. Author contributions: WS designed the waveform tests, processed the data, developed the software, "
    "carried out the simulations and analyses, prepared the figures and wrote the original draft; PL conceived the "
    "lognormal power-law crustal model, contributed to the methodology and reviewed and edited the manuscript.")
d.h("Data availability", 1, numbered=False)
d.p("The event catalogue [@dyer2022] and the downhole geophone records of wells 56-32, 58-32 and 78B-32 from the "
    "2022 stimulation [@pankow2023] are publicly available from the U.S. Department of Energy Geothermal Data "
    "Repository. The SEG-2 versions of the 56-32 records used here were provided by Geo-Energie Suisse. All code "
    "(finite-difference solver, reciprocity, feature extraction, tests and figure scripts, with fixed random seeds) "
    "and all derived result tables are available at https://github.com/maswiet/Paper2_PinkNoise_WaveformInversion, "
    f"release {RELEASE}.")
d.h("Conflict of interest", 1, numbered=False)
d.p("The authors declare that they have no conflict of interest.")
d.h("Supporting information", 1, numbered=False)
d.p("Supplementary data are available at GJI online: Figures S1–S11 and Tables S1–S18 (reciprocity validation, sonic "
    "log, earlier statistical tests, bandwidth and multi-octave diagnostics, location stress test of the coda "
    "coherence, window and moment-tensor sensitivity, S-timing subsets, calibration bootstrap, source spectra, "
    "ringing tests, stage comparison, field statistics, numerical convergence, σ_{eff} curves, complete Φ and "
    "holdout tables, equal-σ truth tests, log-spectrum variants, per-event results for every medium and seed, and "
    "coda-coherence tables). Please note: Oxford University Press is not responsible for the content or "
    "functionality of any supporting materials supplied by the authors. Any queries (other than missing material) "
    "should be directed to the corresponding author for the paper.")
unused = d.references()
print("unused refs:", unused)
print(f"summary words: {NW_ABS}; figures: {d.fig_n}; tables: {d.tab_n}")
print("PENDING:", sorted(set(PENDING)))
d.save(P / "manuscript" / ("Paper2_GJI_review.docx" if d.embed else "Paper2_GJI.docx"))
