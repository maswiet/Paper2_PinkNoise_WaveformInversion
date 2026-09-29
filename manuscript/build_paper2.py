"""Paper 2 manuscript (Geothermal Energy / SpringerOpen style), built from results/*.csv.
Run:  python build_paper2.py            -> Paper2_manuscript.docx (figures uploaded separately)
      EMBED=1 python build_paper2.py    -> Paper2_manuscript_review.docx (figures embedded)
"""
import os, sys
from pathlib import Path
import pandas as pd

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))  # local copy of springer_doc.py (from Paper 1)
from springer_doc import SDoc  # noqa: E402

R = P / "results"; F = P / "figs"
rd = lambda f: pd.read_csv(R / f)
K = rd("key_numbers.csv").set_index("key").value
MD = rd("table_misfit_data_forge.csv").set_index("model")
MR = rd("table_misfit_data_random.csv").set_index("model")
M1 = rd("table_misfit_T1_forge.csv").set_index("model")
M2 = rd("table_misfit_T2_forge.csv").set_index("model")
MT = {o: rd(f"table_mt_{o}.csv").set_index("model") for o in ("T1", "T2", "data")}
TT = rd("table_tt.csv").set_index("model")
TC = rd("table_ttcorr.csv")
FE = rd("table_features.csv")

f2, f3 = (lambda x: f"{x:.2f}"), (lambda x: f"{x:.3f}")
CR = rd("table_coda_ratio.csv").set_index("model")
LV = rd("table_log_vs_media.csv").set_index("medium")
det = [m for m in ("H", "L", "A") if m in MD.index]
comp = [m for m in ("E15", "E50", "G15", "G50") if m in MD.index]
pink = [m for m in MD.index if m not in det + comp]
best = MD.loc[pink].PHI.idxmin(); bestd = MD.loc[det].PHI.idxmin(); bestc = MD.loc[comp].PHI.idxmin()
sig_best = MD.loc[best, "sigma"]
real045 = MD.loc[[m for m in ("P1", "T1", "R1", "R2") if m in MD.index]].PHI
real09 = MD.loc[[m for m in ("P3", "R3") if m in MD.index]].PHI
real18 = MD.loc[[m for m in ("P9", "R5") if m in MD.index]].PHI
fe = lambda s, b, m, c: FE[(FE.sensor == s) & (FE.band_Hz == b) & (FE.model == m)][c].iloc[0]
cr = lambda m, c: CR.loc[m, c]
lv = lambda name: LV.loc[name, "rms_log10_vs_log"]
ST = rd("table_shape_test.csv")
st = lambda tr, m, c="PHI": ST[(ST.truth == tr) & (ST.model == m)][c].iloc[0]
stbest = lambda tr, c="PHI": ST[(ST.truth == tr) & ~ST.model.isin(["H", "A"])].sort_values(c).model.iloc[0]
CO = rd("table_coherence.csv").set_index("model")
CF = rd("table_confusion_runs.csv")
DFm = rd("table_data_families.csv").set_index("family")
FC = rd("table_field_corr.csv").set_index("medium").corr_with_pink_seed33
rexp = FC["exponential a=15 m (seed 33)"]; rgau = FC["Gaussian a=50 m (seed 33)"]
n_ok = int((CF.truth_family == CF.selected).sum())
pl = CF.truth_family.isin(["pink", "exponential"])
n_pl = int((pl & CF.selected.isin(["pink", "exponential"])).sum() + (~pl & (CF.selected == "Gaussian")).sum())
n_pe = int(((CF.truth_family == CF.selected) & pl).sum())
BW = rd("table_bandwidth.csv").set_index("sensor")
CMO = rd("table_coda_multioctave.csv")
CCB = rd("table_coda_cleanband.csv")
QS = rd("table_qc_slope.csv")
ccb = lambda s, m, c: CCB[(CCB.sensor == s) & (CCB.measure == m)][c].iloc[0]
cmo = lambda s, f, c: CMO[(CMO.sensor == s) & (CMO.fc_Hz == f)][c].iloc[0]
qs_fam = lambda fam: QS[QS.family == fam].slope_25_50
qs_d = QS[QS.model == "data"].iloc[0]
BT = rd("table_break_test.csv").set_index("model")
bt = lambda m, c: BT.loc[m, c]
btx = BT[BT.index.str.startswith("EX")]                      # exponential correlation-length series
fam_spread = (CF[["score_pink", "score_exponential", "score_Gaussian"]].max(axis=1)
              - CF[["score_pink", "score_exponential", "score_Gaussian"]].min(axis=1))
TITLE = ("Utah FORGE microearthquake waveforms reject smooth crustal models: 3-D elastic tests of a pink-noise "
         "crust against layered, anisotropic and correlation-length random media")

d = SDoc()
d.embed = bool(os.environ.get("EMBED"))
d.p(f"**{TITLE}**", align="center", size=14)
d.p("Wiwit Suryanto^{1*}, Peter Leary^{2}, [CO-AUTHORS TO BE CONFIRMED]", align="center")
for a in ["^{1} Department of Physics, Faculty of Mathematics and Natural Sciences, Universitas Gadjah Mada, Sekip Utara BLS 21, Yogyakarta 55281, Indonesia",
          "^{2} Geoflow Imaging, Auckland 1010, New Zealand"]:
    d.p(a, align="left", size=11)
d.p("^{*} Corresponding author: Wiwit Suryanto, ws@ugm.ac.id", align="left", size=11)

# ================================================================= ABSTRACT
d.h("Abstract", 1)
d.p(f"Well logs, cores and flow data from sedimentary, crystalline and volcanic crust share a spatial-fluctuation "
    f"signature: power spectra scaling as *S*(*k*) ∝ 1/*k* over cm–km scales and lognormal permeability, "
    f"*κ* ∝ exp(*αφ*). If the seismogenic crust is a pink-noise medium in this sense, the microearthquakes (MEQs) "
    f"induced by geothermal stimulation should carry its imprint. We test this with 3-D elastic finite-difference "
    f"(FD) wavefields and the April 2022 stimulation of well 16A(78)-32 at Utah FORGE. A pink-noise crust is compared "
    f"with homogeneous, layered and transversely isotropic (VTI) crusts, and with random media that have an "
    f"exponential or Gaussian correlation function. {int(K.n_events_box)} catalogued MEQs recorded by two downhole "
    f"sensors in well 56-32 were extracted from continuous SEG-2 records. The sensors were self-calibrated from about "
    f"1000 P onsets (depths {K.sensorA_depth_m:.0f} and {K.sensorB_depth_m:.0f} m; *V*_{{P}} = {K.vp_m_s:.0f} m/s). "
    f"Strain Green's tensors computed by reciprocity give synthetics for every event–sensor path from six FD runs per "
    f"model, with realistic mechanisms and real noise. A deterministic moment-tensor waveform inversion prefers the "
    f"smoothest model even when the synthetic truth is pink noise, so it cannot identify the medium class. A "
    f"statistical inversion of S-coda features and envelopes recovers a pink truth and a layered-VTI truth without "
    f"bias. For FORGE, all smooth crusts fail. With moment tensors fitted to the direct waves of 150 observed events, "
    f"the VTI crust predicts S coda {10**-cr('A','med_log10_pred_over_obs_A'):.0f} times too weak (median, sensor A), "
    f"and only {cr('A','frac_within_x3_A')*100:.0f} % of events fall within a factor of three. Smooth crusts also keep "
    f"the coda of neighbouring events coherent (median {CO.loc['H','cohA_0_5m']:.2f}→{CO.loc['H','cohA_80_160m']:.2f} "
    f"from 0–5 to 80–160 m separation), whereas the observed coda decorrelates "
    f"({CO.loc['data','cohA_0_5m']:.2f}→{CO.loc['data','cohA_80_160m']:.2f}) as in scattering media. The statistical misfit is "
    f"{f2(MD.loc[bestd,'PHI'])}–{f2(MD.loc[det].PHI.max())} for smooth crusts and {f2(MD.loc[best,'PHI'])} for the best "
    f"pink crust (σ = {sig_best:.2f}, bracketed by runs at σ = 0.13 and 0.25). The seismic data do not, however, "
    f"establish spectral shape. An equal-σ truth test with three realisations per family classifies a pink, "
    f"exponential or Gaussian truth correctly in {n_ok} of 9 cases. It separates small-scale-rich (pink or "
    f"exponential) from Gaussian media in {n_pl} of 9 cases, but pink from exponential only at chance "
    f"({n_pe} of 6), because on the simulated scales the two fields are {rexp*100:.0f} % correlated. With two "
    f"realisations per family, the FORGE misfits (pink {f2(DFm.loc['pink','mean'])}, exponential "
    f"{f2(DFm.loc['exponential','mean'])}, Gaussian {f2(DFm.loc['Gaussian','mean'])}) differ by less than the "
    f"realisation scatter, so the data prefer no family. Shape is constrained instead by the 56-32 sonic log. Its "
    f"12–200 m spectrum is scale-free and is matched by pink noise (rms {lv('pink p=1.5'):.2f} in log power), less "
    f"well by exponential media ({lv('exponential a=15 m'):.2f}–{lv('exponential a=50 m'):.2f}), and not at all by "
    f"Gaussian media (≥ {lv('Gaussian a=15 m'):.1f}). The 4-kHz downhole records are broadband: the S wave "
    f"exceeds noise by a factor of three up to about 1.8 kHz, far above the 100-Hz Nyquist frequency of 200-sample/s "
    f"recording. Over the resonance-free band 25–200 Hz, the single-scattering coda attenuation follows one power law, "
    f"Q_{{c}}^{{-1}} ∝ *f*^{{{ccb('A','QI','slope'):.2f}}} (sensor A) and *f*^{{{ccb('B','QI','slope'):.2f}}} (sensor B). "
    f"Simulations with correlation lengths of 5–134 m show, however, that such a power law is not diagnostic of "
    f"scale-free heterogeneity. A 2.5-m-grid simulation series extended to 100 Hz is more selective. At 100 Hz the "
    f"observed coda level ({bt('data','CL100_B'):.2f} in log_{{10}} of coda/S energy) is close to that of pink media "
    f"({min(bt('R6','CL100_B'), bt('P6','CL100_B')):.2f} to {max(bt('R6','CL100_B'), bt('P6','CL100_B')):.2f}) and is "
    f"bracketed by exponential media ({btx.CL100_B.min():.2f} to {btx.CL100_B.max():.2f}, decreasing with *a*), "
    f"whereas the Gaussian medium with *a* = 50 m "
    f"({min(bt('G50','CL100_B'), bt('G50b','CL100_B')):.2f} to {max(bt('G50','CL100_B'), bt('G50b','CL100_B')):.2f}) and "
    f"the homogeneous crust ({bt('H','CL100_B'):.2f}) lose their coda. The MEQ "
    f"wavefields thus reject the smooth crustal models of routine microseismic practice and require a strongly "
    f"scattering crust. The power-law (pink) character of that crust is supported by the borehole data, not "
    f"independently by the seismic wavefield.")
d.p("**Keywords:** EGS; Utah FORGE; microseismicity; seismic scattering; pink noise; 3-D finite difference; "
    "coda; statistical waveform inversion", align="left")

# ================================================================= INTRODUCTION
d.h("Introduction", 1)
d.p("Engineered geothermal system (EGS) stimulation is monitored almost entirely through induced microseismicity, "
    "and the interpretation of that microseismicity rests on a crustal model. Standard practice treats the "
    "reservoir as homogeneous, layered or weakly anisotropic: the model is smooth enough that ray theory holds, "
    "travel times locate events and waveforms are fitted by a few source parameters. A different view emerges from "
    "borehole data. Well-log fluctuation spectra in sedimentary, crystalline and volcanic settings scale as "
    "*S*(*k*) ∝ 1/*k*^{β} with β ≈ 1 over cm–km scales, and well-core porosity *φ* and permeability *κ* follow "
    "*κ* ∝ exp(*αφ*) with large *α* [@leary1997; @leary2002; @holliger1996; @leary_abg]. The ambient crust is then a "
    "spatially correlated, lognormally distributed \"pink-noise\" medium without a characteristic scale. This is "
    "the physical basis of geoflow imaging (GFI), in which MEQ emissions map the flow structure of that medium "
    "[@leary_abg; @leary2026].")
d.p("Paper 1 of this series quantified the α, β and γ descriptors of flow heterogeneity at Utah FORGE from logs, "
    "cores, flow and catalogue data [@suryanto2026]. The present paper asks whether the *seismic wavefield* of "
    "FORGE stimulation MEQs carries the imprint of a pink-noise crust, and whether that imprint can be "
    "distinguished from conventional smooth crustal models. Short-period seismic scattering by random media is well "
    "established [@aki1975; @frankel1986; @sato2012]. However, to our knowledge no study has used FORGE downhole MEQ "
    "waveforms to test a well-log-constrained pink-noise model against homogeneous, layered and VTI alternatives "
    "within a single 3-D elastic framework and with synthetic truth tests.")
d.p("We make four contributions: (i) a reproducible extraction and self-calibration of the April 2022 56-32 "
    "downhole records; (ii) a reciprocity-based 3-D elastic FD framework that yields full-waveform synthetics for "
    "hundreds of event–sensor paths per model; (iii) a demonstration, with synthetic truths, that deterministic "
    "waveform inversion cannot identify the medium class whereas statistical (coda/envelope) waveform inversion "
    "can; and (iv) a direct comparison of observed and predicted FORGE seismograms for smooth, correlation-length "
    "and pink-noise crusts, combined with the independent constraint of the sonic log.")

# ================================================================= CONCEPT
d.h("The pink-noise crust", 1)
d.p("Following the GFI recipe [@leary_abg], a 3-D pink-noise field *f*(**x**) is generated by filtering Gaussian "
    "white noise in the wavenumber domain,")
d.eq(r"\hat{f}(\mathbf{k}) = \frac{\hat{w}(\mathbf{k})}{\left(1+|\mathbf{k}|\right)^{p}}")
d.p("normalised to zero mean and unit variance. Seismic velocities follow the lognormal form of the poro-permeability "
    "relation,")
d.eq(r"V_{P}(\mathbf{x}) = \bar{V}_{P}\,e^{\sigma f(\mathbf{x})},\qquad V_{S}(\mathbf{x}) = \bar{V}_{S}\,e^{\sigma f(\mathbf{x})}")
d.p(f"with σ the standard deviation of ln *V* at grid scale. One-dimensional lines through the field reproduce "
    f"well-log spectra: *p* = 1.5 gives *S*(*k*) ∝ *k*^{{-1.05}} (the canonical 1/*k*) and *p* = 1.8 gives "
    f"*k*^{{-1.48}}, the slope measured on the 56-32 sonic log (β = {K.log_beta:.2f}). The conventional "
    f"alternatives are a homogeneous crust (H), a layered crust built from 50-m blocks of the 56-32 sonic log (L), "
    f"and a homogeneous VTI crust (A) whose Thomsen parameters are inverted from P travel times [@thomsen1986].")
d.p("To test whether the *spectral shape* matters, and not only the strength of heterogeneity, we also use "
    "conventional random media with a correlation length *a* [@sato2012]. Their fields are generated from the same "
    "white noise with amplitude spectra")
d.eq(r"\hat{f}_{\mathrm{exp}} \propto \frac{\hat{w}}{1+k^{2}a^{2}},\qquad \hat{f}_{\mathrm{Gauss}} \propto \hat{w}\,e^{-k^{2}a^{2}/8}")
d.p("corresponding to exponential and Gaussian autocorrelation functions. We use *a* = 15 and 50 m and set σ = 0.13 "
    "with the same ±3σ truncation as the pink model P6, so that the two families differ only in spectral shape.")

# ================================================================= DATA
d.h("Data", 1)
d.h("Stimulation MEQs and continuous downhole records", 2)
d.p(f"Geo-Energie Suisse recorded the April 2022 stimulation of well 16A(78)-32 (stages 1–3) in 60-s SEG-2 files "
    f"(4 kHz, six channels, GPS-timed) from monitoring well 56-32, and published a catalogue of 36 641 located "
    f"events [@dyer2022]. We read the acquisition time of 3556 files (18–24 April 2022) and cut 1-s windows at "
    f"catalogue origin times for events with *M* ≥ −1, which gave 800 events. Of these, {int(K.n_events_box)} "
    f"with a clean P onset lie inside the simulation volume (Fig. 1).")
d.h("Sensor self-calibration", 2)
d.p("The file headers carry no receiver coordinates, and the 2024 56-32 sensor position predicts P arrivals ~0.1 s "
    "too late. Automatic AIC onsets [@maeda1985] on the two three-component sensors (channels 1–3 and 4–6) were "
    "therefore inverted jointly with the catalogue hypocentres for sensor position, clock offset and *V*_{P}. "
    "An L1 misfit was used, the sensors were constrained to one vertical well and outliers were rejected. The "
    f"solution places both sensors at E = {K.sensor_E_m:.0f} m, N = {K.sensor_N_m:.0f} m (16A wellhead frame), "
    f"depths {K.sensorA_depth_m:.0f} m (A) and {K.sensorB_depth_m:.0f} m (B), with *V*_{{P}} = {K.vp_m_s:.0f} m/s "
    f"and residual MAD {TT.loc['H','MAD_ms']:.2f} ms. The S–P time grows with distance at 1.31 × 10^{{-4}} s/m, "
    f"consistent with *V*_{{P}}/*V*_{{S}} ≈ 1.73–1.76, which confirms the phase identification. Channels 1 and 4 "
    f"carry only ~3 % of the P energy and have much lower SNR, whereas channels 2–3 and 5–6 behave as horizontal "
    f"pairs (P rectilinearity 0.82 and 0.98). The analysis therefore uses rotation-invariant horizontal energy.")
d.h("Sonic log", 2)
d.p(f"The 56-32 ThruBit sonic porosity (1950–2780 m, calliper-screened) converts to a mean *V*_{{P}} of 5895 m/s, "
    f"close to the travel-time value. The detrended ln *V* has σ = {K.log_sig6:.3f} at 6-m and "
    f"{K.log_sig50:.3f} at 50-m averaging and a spectrum *S*(*k*) ∝ *k*^{{-{K.log_beta:.2f}}} over 1–200 m "
    f"(Fig. 2). These values set the reference pink-noise model P1 (σ = 0.045) and bracket the model grid.")
d.figure("Utah FORGE 2022 geometry: stimulation MEQs, 56-32 sensors and FD volume",
         f"a) {int(K.n_events_box)} catalogued MEQs (colour = stimulation stage), the self-calibrated sensors A and B "
         "(red triangles) and the interior of the FD volume (dotted box). b) Distribution of hypocentral distances "
         "to sensors A and B.", F / "fig_geometry.png")
d.figure("Velocity fluctuation statistics of the 56-32 sonic log",
         f"a) Sonic-derived *V*_{{P}} (grey) and its 50-m block average (black); dashed lines mark the depths of "
         f"sensors A and B. b) Distribution of the detrended, 6-m-averaged δln *V* (σ = {K.log_sig6:.3f}). "
         f"c) Power spectrum of δln *V* with the power-law fit *S*(*k*) ∝ *k*^{{-{K.log_beta:.2f}}} (red).", F / "fig_log5632.png")

# ================================================================= METHODS
d.h("Methods", 1)
d.h("3-D elastic finite differences", 2)
d.p(f"We solve the velocity–stress equations on a staggered grid [@virieux1986] with fourth-order spatial and "
    f"second-order temporal differences [@levander1988] for a VTI medium (stiffnesses from *V*_{{P}}, *V*_{{S}}, "
    f"ρ and Thomsen ε, δ, γ). The treatment of anisotropy on the staggered grid follows Igel, Mora and Riollet "
    f"[@igel1995], and the MATLAB implementation grew out of the acoustic finite-difference codes of Igel "
    f"[@igel2016], which we extended to 3-D elastic VTI media, moment-tensor and point-force sources and "
    f"strain-rate recording. A Cerjan sponge [@cerjan1985] absorbs outgoing energy. The grid is "
    f"{int(K.nx)}×{int(K.ny)}×{int(K.nz)} nodes at {K.dx_m:.0f} m (0.46 × 0.58 × 0.64 km interior plus 20-node sponges), "
    f"with Δ*t* = {K.dt_s*1e3:.2f} ms (the data sampling), {int(K.nt)} steps and a {K.f0_Hz:.0f}-Hz Ricker wavelet; "
    f"the analysis band is 20–80 Hz (≥ 5 grid points per minimum S wavelength). In a homogeneous test the solver "
    f"reproduces the P velocity to 0.6 %.")
d.h("Reciprocity and strain Green's tensors", 2)
d.p("Instead of one run per event, three point-force runs per sensor record the strain-rate tensor at every event "
    "node. By reciprocity the velocity at sensor component *n* from a moment tensor **M** at event **x**_{e} is")
d.eq(r"v_{n}(\mathbf{x}_{s},t) = M_{pq}\int_{0}^{t}\dot{\varepsilon}^{(n)}_{pq}(\mathbf{x}_{e},\tau)\,d\tau")
d.p("[@zhao2006]. A validation in a heterogeneous VTI pink-noise medium reproduces direct moment-tensor "
    "simulations with correlation 0.995 and amplitude ratio 0.995 (Fig. 3). Six runs per model thus "
    "give synthetics for all ~1100 event–sensor paths with any mechanism.")
d.figure("Reciprocity validation in a heterogeneous VTI pink-noise medium",
         "a–c) East, north and vertical velocity at a sensor from a direct double-couple simulation (black) and from "
         "the strain-Green's-tensor reciprocity synthetic (red dashed), in a pink-noise medium with σ = 0.045, "
         "ε = 0.05, δ = 0.02 and γ = 0.05.",
         F / "val_reciprocity.png")
d.h("Model suite", 2)
rows = [["H", "homogeneous", "*V*_{P} 5808 m/s, *V*_{P}/*V*_{S} 1.73"],
        ["L", "layered", "56-32 sonic, 50-m blocks"],
        ["A", "VTI", "travel-time inverted ε, δ"],
        ["P4, P1, P3, P6, P9, P12, P10", "pink, *p* = 1.5", "σ = 0.02, 0.045, 0.09, 0.13, 0.18, 0.22, 0.25"],
        ["P5, P2", "pink, σ = 0.045", "*p* = 1.2, 1.8"],
        ["P7, P8", "pink, σ = 0.09", "*p* = 1.2, 1.8"],
        ["P11", "pink, σ = 0.18", "*p* = 1.8"],
        ["T1, R1, R2", "pink σ = 0.045, *p* = 1.5", "independent realisations (T1 = synthetic truth)"],
        ["R3; R5", "pink, *p* = 1.5", "independent realisations of σ = 0.09; 0.18"],
        ["E15, E50", "exponential ACF", "σ = 0.13, *a* = 15, 50 m"],
        ["G15, G50", "Gaussian ACF", "σ = 0.13, *a* = 15, 50 m"],
        ["T2", "layered + VTI", "synthetic truth: ε = 0.06, δ = 0.03, γ = 0.06"]]
d.table("Crustal models simulated", ["Model", "Class", "Parameters"], rows, [4.2, 4.0, 7.8],
        legend="Fields with σ > 0.1 are truncated at ±3σ (±2.5σ for σ ≥ 0.18). For σ ≥ 0.22 the time step is "
               "halved for stability; at σ = 0.25 the slowest S waves have only 4.5 grid points per wavelength at "
               "80 Hz.")
d.p("Synthetic records use double-couple mechanisms consistent with FORGE stress (strike 0–50°, dip 55–85°, "
    "rake −90 ± 30°); a second set uses fully random mechanisms. Real 56-32 noise is added to each record at "
    "that event's observed SNR, and the data are convolved with the same Ricker wavelet, so data and synthetics "
    "share a source spectrum (FORGE MEQ corner frequencies lie well above 80 Hz).")
d.h("Deterministic waveform inversion", 2)
d.p("For each of 150 events we invert the six moment-tensor elements linearly from the four horizontal traces, "
    "with 20–80 Hz filtering, a window from P − 10 ms to S + 60 ms and a per-sensor time shift of ±6 ms. For the "
    "data, the tool azimuth and handedness of each sensor are grid-searched. Fit quality is the variance "
    "reduction (VR).")
d.h("Statistical waveform inversion", 2)
d.p("Because the realisation of a random medium is unknowable from two sensors, we invert the *statistics* of "
    "the wavefield. On the horizontal energy envelope we measure three mechanism-robust, S-normalised features: "
    "coda level log_{10}(*E*_{coda}/*E*_{S}) 30–100 ms after the S peak, coda decay rate, and rms S-pulse width. "
    "The misfit between observed and model features is")
d.eq(r"\Phi = \overline{D_{\mathrm{KS}}} + \overline{\left\| \log_{10}\tilde{E}_{\mathrm{obs}} - \log_{10}\tilde{E}_{\mathrm{mod}} \right\|_{\mathrm{rms}}}")
d.p("where *D*_{KS} is the two-sample Kolmogorov–Smirnov distance per feature, sensor and band, and *Ẽ* is the "
    "median post-peak S envelope. Sensor B at 20–40 Hz is excluded because it carries a persistent non-event "
    "plateau, likely a borehole resonance. P travel-time residual statistics are reported separately, since "
    "catalogue location error dominates them in the data.")
d.h("Coda prediction test", 2)
d.p("The moment tensor fitted to the direct waves (P − 10 ms to S + 60 ms) of each of the 150 events is used to "
    "predict the complete record, including the S coda 30–100 ms after the S peak. That coda lies outside the fitting "
    "window and is therefore a pure prediction of the medium. For each event and model we compare the predicted and "
    "observed coda-to-S energy ratios at 40–80 Hz.")
d.h("Shape-sensitive diagnostics", 2)
d.p("Two further observables are chosen because scale-free and correlation-length media should differ in them. "
    "D1 is the frequency dependence of the coda: per event, the coda level at 40–80 Hz minus that at 20–40 Hz "
    "(sensor A). D2 is the inter-event coherence of the coda: for pairs of events recorded on the same sensor, the "
    "maximum normalised two-component cross-correlation (lag ±8 ms) of the 40–80 Hz window S + 15 … S + 100 ms, "
    "summarised as the median in bins of inter-event separation (0–5 to 80–160 m). To test whether these "
    "observables and Φ can resolve spectral shape, nine further synthetic truths were simulated at equal "
    "σ = 0.13: pink, exponential with *a* = 15 m and Gaussian with *a* = 50 m, each with seeds 33, 44 and 55. "
    "Candidates are two realisations per family (seeds 11 and 33). A candidate that shares the truth's seed is "
    "excluded, in which case only seed 11 is used for every family. The family score is the mean Φ over its "
    "candidate realisations, and the truth is assigned to the family with the lowest score.")
d.h("Multi-octave coda analysis of the raw records", 2)
d.p("The raw 4-kHz records (no source-wavelet convolution) are analysed in octave bands centred at 25, 50, 100, "
    "200, 400 and 800 Hz, plus 1000–1700 Hz. For each event, sensor and band we measure four quantities on the "
    "horizontal-energy envelope. The first is the coda level, defined as the energy 30–100 ms after the S peak "
    "divided by the S-peak energy. The second is the energy decay rate, taken from 20 to 120 ms after the S peak. "
    "The third is the single-scattering coda attenuation Q_{c}^{-1} [@aki1975], from the slope of "
    "ln(*E* *t*^{2}) against lapse time *t*, where *E* *t*^{2} ∝ exp(−2π*f**t*/Q_{c}). The fourth is the rms "
    "width of the S pulse. A band is used only when the coda energy exceeds five times the pre-event noise "
    "energy in that band. All four measures are energy ratios within one band, so a linear instrument or coupling "
    "response cancels. The frequency dependence is tested by comparing a single power law in *f* with a continuous "
    "two-segment law, using BIC. A medium with a correlation length *a* should show a change of behaviour near "
    "*ka* ≈ 1 (about 36 Hz for *a* = 15 m and 11 Hz for *a* = 50 m at *V*_{S} ≈ 3.4 km/s), whereas a scale-free "
    "medium should not. The same measures are computed on the FD synthetics in the 25 and 50 Hz bands, which the "
    "5-m grid resolves (sensors A and B).")
d.h("High-frequency simulations and correlation-length series", 2)
d.p("To reach 100 Hz the FD grid is refined to 2.5 m (169×209×153 nodes, Δ*t* = 0.125 ms, 2200 steps, Ricker "
    "90 Hz; at least 5.7 grid points per minimum S wavelength at 160 Hz). To keep the cost manageable the volume is "
    "reduced to sensor B and the core of the event cloud, which retains 542 events at 92–392 m. The records are "
    "decimated to the 4-kHz data sampling. This grid is run for the homogeneous crust, pink σ = 0.13 (two "
    "realisations), exponential *a* = 15 m (two), Gaussian *a* = 50 m (two) and an exponential correlation-length "
    "series with *a* = 5, 9, 17, 34, 67 and 134 m (σ = 0.13). The series *a* = 9–134 m is also run on the 5-m grid. "
    "Together the two grids give Q_{c}^{-1} and coda level at 25 and 50 Hz (sensor A, 5-m grid) and at 50 and 100 Hz "
    "(sensor B, 2.5-m grid). The curvature, slope(50–100 Hz) − slope(25–50 Hz), measures a change of slope of the "
    "kind expected near *ka* ≈ 1. Synthetic records on the 2.5-m grid are noise-free.")

# ================================================================= RESULTS
d.h("Results", 1)
d.h("Travel times", 2)
d.p(f"P onsets fit a homogeneous crust with MAD {TT.loc['H','MAD_ms']:.2f} ms. A vertical gradient "
    f"(ΔAIC = {TT.loc['G','AIC']-TT.loc['H','AIC']:.1f}, ΔBIC = {TT.loc['G','BIC']-TT.loc['H','BIC']:.1f}) and VTI "
    f"(ΔAIC = {TT.loc['VTI','AIC']-TT.loc['H','AIC']:.1f}, ΔBIC = {TT.loc['VTI','BIC']-TT.loc['H','BIC']:.1f}) barely "
    f"improve the fit. However, the residuals are spatially correlated between neighbouring events "
    f"(ρ = {TC.rho_A[0]:.2f}/{TC.rho_B[0]:.2f} at < 10 m and {TC.rho_A[1]:.2f}/{TC.rho_B[1]:.2f} at 10–20 m for "
    f"sensors A/B, falling to zero beyond 40 m; Fig. 4). Such correlation is not predicted by any smooth model; it could "
    f"arise from structure or from correlated location error.")
d.figure("FORGE 56-32 records, S envelopes and travel-time residual correlation",
         "a, d) Horizontal records (H1, 40–80 Hz) of 40 events ordered by hypocentral distance for sensors A and B; the "
         "dashed line is the predicted S time. b, e) Median (solid) and interquartile range (dotted) of the horizontal "
         "S envelope, 40–80 Hz, for sensors A and B. c) Correlation of P travel-time residuals (high-frequency picks, "
         "homogeneous model) versus inter-event separation.",
         F / "fig_data_overview.png")
d.h("Synthetic truth tests", 2)
d.p(f"With a pink-noise truth (T1), the deterministic inversion prefers the smooth models: VR "
    f"{f2(MT['T1'].loc['H','VR_med'])} (H) and {f2(MT['T1'].loc['L','VR_med'])} (L), against "
    f"{f2(MT['T1'].loc['P1','VR_med'])} for P1, which has the correct statistics but a different realisation. VR "
    f"falls monotonically with σ ({f2(MT['T1'].loc['P3','VR_med'])} for σ = 0.09, "
    f"{f2(MT['T1'].loc['P6','VR_med'])} for 0.13), because a wrong realisation places scattered energy at wrong "
    f"times. The statistical inversion instead ranks a pink model with the true σ first (best {M1.PHI.idxmin()}, "
    f"an independent realisation of the truth statistics, Φ = {f2(M1.PHI.min())}). Pink models with σ = 0.02–0.045 "
    f"and *p* ≤ 1.5 follow (Φ ≤ {f2(M1.loc[['P4','P1','R2','P5']].PHI.max())}), ahead of H, L and A "
    f"({f3(M1.loc[['H','L','A']].PHI.min())}–{f3(M1.loc[['H','L','A']].PHI.max())}). All models with σ ≥ 0.09 are "
    f"rejected (Φ ≥ {f2(M1.loc[M1.sigma >= 0.09].PHI.min())}). Realisation scatter at the true parameters "
    f"(Φ = {f2(M1.loc['R1','PHI'])}–{f2(M1.loc[['P1','R2']].PHI.max())}) limits the σ resolution to about a factor of two. With a layered-VTI truth (T2), L and A are ranked first "
    f"(Φ = {f2(M2.loc['L','PHI'])}, {f2(M2.loc['A','PHI'])}) and every pink model is rejected "
    f"(Φ ≥ {f2(M2.loc[[m for m in M2.index if m not in ('H','L','A')]].PHI.min())}). The method therefore resolves the "
    f"medium class and σ, is not biased toward pink noise, and resolves the spectral exponent *p* only weakly (the correct-σ model with "
    f"*p* = 1.8 ties with the smooth models, Φ = {f2(M1.loc['P2','PHI'])}). The correlation-length media were "
    f"simulated only at σ = 0.13 and are rejected for both truths (Φ ≥ {f2(M1.loc[comp].PHI.min())} for T1, "
    f"≥ {f2(M2.loc[comp].PHI.min())} for T2) because of their σ. Resolution of spectral shape at equal σ is tested "
    f"separately below.")
d.figure("Model ranking by statistical waveform misfit for synthetic truths and FORGE",
         "a) Synthetic truth T1 (pink, σ = 0.045, *p* = 1.5). b) Synthetic truth T2 (layered + VTI). c) Utah FORGE "
         "2022 data (56-32 sensors). Grey: smooth crusts (H, L, A); orange: pink-noise crusts; blue: exponential and "
         "Gaussian random media (σ = 0.13). Lower Φ is better.",
         F / "fig_misfit_ranking.png")
d.h("Utah FORGE", 2)
d.p(f"On the FORGE data the deterministic inversion again prefers the smooth crust (VR "
    f"{f2(MT['data'].loc['H','VR_med'])} for H to {f2(MT['data'].loc['P6','VR_med'])} for P6), yet it leaves "
    f"coda energy equal to {MT['data'].loc['H','unexpl_coda_over_S']*100:.0f} % of the S-window energy unexplained. The same "
    f"smooth model leaves only {MT['T1'].loc['H','unexpl_coda_over_S']*100:.0f}–{MT['T2'].loc['H','unexpl_coda_over_S']*100:.0f} % "
    f"unexplained when fitted to the synthetic truths. The statistical inversion ranks every pink crust with σ ≥ 0.045 ahead of every "
    f"conventional crust (Figs. 5 and 7). Φ decreases from {f2(MD.loc['H','PHI'])} (H), {f2(MD.loc['L','PHI'])} (L) and "
    f"{f2(MD.loc['A','PHI'])} (A) through {f2(real045.min())}–{f2(real045.max())} (σ = 0.045, "
    f"{len(real045)} realisations) and {f2(real09.min())}–{f2(real09.max())} (σ = 0.09) to "
    f"{f2(MD.loc['P6','PHI'])} (σ = 0.13) and {f2(real18.min())}–{f2(real18.max())} (σ = 0.18, two realisations). "
    f"It then rises to {f2(MD.loc['P12','PHI'])} (σ = 0.22) and {f2(MD.loc['P10','PHI'])} (σ = 0.25), which brackets "
    f"the optimum near σ ≈ 0.18. At σ = 0.18 the exponent *p* = 1.8 fits worse than 1.5 ({f2(MD.loc['P11','PHI'])}). "
    f"Realisation spread is up to {max(real045.max()-real045.min(), real09.max()-real09.min()):.2f} at σ ≤ 0.09 but "
    f"only {real18.max()-real18.min():.3f} at σ = 0.18. At σ = 0.045 the worst realisation comes within "
    f"{f2(MD.loc[det].PHI.min()-real045.max())} of the best smooth crust, so the evidence is decisive only for "
    f"σ ≥ 0.09. With random mechanisms the smooth crusts again rank last (Φ ≥ {f2(MR.loc[det].PHI.min())}).")
d.p(f"The physics is visible in the envelopes (Fig. 6). Smooth models predict an S pulse followed by a coda two "
    f"orders of magnitude below the peak within 50 ms. The data retain 10^{{-2}}–10^{{-1}} of the peak energy for "
    f"≥ 100 ms. At 40–80 Hz the median coda level is {fe('A','40-80','data','coda_level_log10'):.2f} (A) and "
    f"{fe('B','40-80','data','coda_level_log10'):.2f} (B), against {fe('A','40-80','H','coda_level_log10'):.2f} and "
    f"{fe('B','40-80','H','coda_level_log10'):.2f} for H. The S pulse is also broader in the data "
    f"({fe('A','40-80','data','S_width_ms'):.1f} ms versus {fe('A','40-80','H','S_width_ms'):.1f} ms for H at "
    f"sensor A). Scattering media with σ ≈ 0.09–0.18 reproduce both. Residual differences are a slower coda decay "
    f"in the data and excess coda at 20–40 Hz; these point to a spectral exponent or intrinsic attenuation not "
    f"captured by the grid.")
d.figure("Median S-coda envelopes: FORGE data versus crustal model classes",
         "Envelopes of horizontal energy aligned on the S peak and normalised by it. a) Sensor A, 20–40 Hz. b) Sensor A, "
         "40–80 Hz. c) Sensor B, 40–80 Hz. Black: FORGE data; H homogeneous, A VTI, E50/G50 exponential/Gaussian "
         "(*a* = 50 m, σ = 0.13), P3/P9 pink (σ = 0.09/0.18).", F / "fig_envelopes.png")
d.figure("Statistical waveform misfit over pink-noise parameters (σ, p)",
         f"Colour and labels: Φ for FORGE data at each (σ, *p*), with mean and range where several realisations exist. "
         f"For comparison, the best smooth crust is {bestd} (Φ = {f2(MD.loc[bestd,'PHI'])}) and the best "
         f"correlation-length medium is {bestc} (Φ = {f2(MD.loc[bestc,'PHI'])}).",
         F / "fig_sigma_p_map.png")
d.h("Observed versus predicted seismograms", 2)
d.p(f"Fig. 8 shows sensor-A records of four events at 359–419 m with the predictions of five crusts; the moment "
    f"tensor was fitted to the direct waves only. In the grey coda window, homogeneous and VTI predictions are flat, "
    f"whereas the observed traces keep oscillating at 10–60 % of the S peak. Pink and correlation-length crusts "
    f"produce coda of comparable amplitude, though not the same wiggles, since the realisation is unknown. On the "
    f"direct waves the smooth models often reach higher VR (median {f2(cr('H','VR_med'))} for H against "
    f"{f2(cr('P6','VR_med'))} for P6), as in the truth tests. Fig. 9 shows the same events as envelopes for both "
    f"sensors. Smooth-model coda falls to 10^{{-3}}–10^{{-4}} of the S peak within 50–80 ms, whereas the observed and "
    f"scattering-model coda stays at 10^{{-2}}–10^{{-1}}.")
rows = []
for m, name in (("H", "homogeneous"), ("L", "layered"), ("A", "VTI"), ("E15", "exponential, a = 15 m"),
                ("E50", "exponential, a = 50 m"), ("G15", "Gaussian, a = 15 m"), ("G50", "Gaussian, a = 50 m"),
                ("P3", "pink, σ = 0.09"), ("P6", "pink, σ = 0.13"), ("P9", "pink, σ = 0.18")):
    if m in CR.index:
        rows.append([f"{m} ({name})", f"{cr(m,'med_log10_pred_over_obs_A'):+.2f}",
                     f"{cr(m,'med_log10_pred_over_obs_B'):+.2f}", f"{cr(m,'frac_within_x3_A')*100:.0f}",
                     f"{cr(m,'frac_within_x3_B')*100:.0f}"])
d.table("Predicted versus observed S-coda ratio for 150 FORGE events",
        ["Model", "Median log_{10}(pred/obs), A", "Median, B", "% within ×3, A", "% within ×3, B"], rows,
        [5.2, 3.2, 2.2, 2.6, 2.6],
        legend="Coda ratio = mean energy 30–100 ms after the S peak divided by S-peak energy, 40–80 Hz horizontals. "
               "Moment tensors were fitted to P − 10 ms … S + 60 ms, so the coda is a prediction.")
d.p(f"Across all 150 events (Table 2, Fig. 10) the smooth crusts under-predict the coda by a factor of "
    f"{10**-cr('L','med_log10_pred_over_obs_A'):.0f}–{10**-cr('H','med_log10_pred_over_obs_A'):.0f} at sensor A and "
    f"{10**-cr('H','med_log10_pred_over_obs_B'):.0f}–{10**-cr('A','med_log10_pred_over_obs_B'):.0f} at sensor B, and only "
    f"{min(cr(m,'frac_within_x3_A') for m in det)*100:.0f}–{max(cr(m,'frac_within_x3_B') for m in det)*100:.0f} % of "
    f"events fall within a factor of three. The VTI crust is indistinguishable from the homogeneous one. The pink "
    f"crust with σ = 0.09 is unbiased (median {cr('P3','med_log10_pred_over_obs_A'):+.2f}/"
    f"{cr('P3','med_log10_pred_over_obs_B'):+.2f}), while σ = 0.13–0.18 over-predicts by a factor of 2–4. This "
    f"per-event test places the effective σ near 0.09–0.13, slightly below the statistical optimum. The "
    f"correlation-length media perform comparably; the Gaussian medium with *a* = 50 m has the most events within a "
    f"factor of three ({cr('G50','frac_within_x3_A')*100:.0f}/{cr('G50','frac_within_x3_B')*100:.0f} %).")
d.figure("Observed FORGE seismograms versus predictions of five crustal models",
         "Sensor A, horizontal H1, 40–80 Hz. Rows: four events at increasing distance (event number and hypocentral "
         "distance on the left). Columns: a, f, k, p) homogeneous; b, g, l, q) VTI; c, h, m, r) exponential, *a* = 15 m; "
         "d, i, n, s) Gaussian, *a* = 15 m; e, j, o, t) pink, σ = 0.13. Black: observed; red: predicted with the moment "
         "tensor fitted to P − 10 ms … S + 60 ms; grey: coda prediction window; dotted: S. The VR of the fit is given "
         "in each panel.", F / "fig_seis_overlay.png")
d.figure("Observed and predicted S envelopes for the events of Fig. 8",
         "Horizontal-energy envelopes normalised by the S peak for the four events of Fig. 8: a–d) sensor A, "
         "e–h) sensor B, events in the same order as the rows of Fig. 8. Black: observed; colours: model predictions "
         "(legend).", F / "fig_seis_envelopes.png")
d.figure("Per-event predicted/observed S-coda ratio for 150 FORGE events",
         "a) Sensor A. b) Sensor B. 40–80 Hz, 150 events. Boxes: interquartile range; whiskers: 10–90 %; diamonds: median; "
         "dotted lines: factor of three. Grey: smooth crusts; blue: correlation-length media (σ = 0.13); orange: pink "
         "crusts.", F / "fig_coda_ratio.png")
d.h("Inter-event coda coherence and frequency dependence", 2)
d.p(f"Fig. 11 shows D1 and D2 for FORGE and the models. Neighbouring FORGE events have coda coherence "
    f"{CO.loc['data','cohA_0_5m']:.2f} (sensor A) and {CO.loc['data','cohB_0_5m']:.2f} (B) at 0–5 m separation, "
    f"falling to {CO.loc['data','cohA_80_160m']:.2f} and {CO.loc['data','cohB_80_160m']:.2f} at 80–160 m. In the "
    f"homogeneous and VTI crusts the coda is a deterministic tail of the direct S wave and stays coherent "
    f"({CO.loc['H','cohA_0_5m']:.2f}→{CO.loc['H','cohA_80_160m']:.2f} for H, "
    f"{CO.loc['A','cohA_0_5m']:.2f}→{CO.loc['A','cohA_80_160m']:.2f} for VTI at sensor A). Every scattering medium "
    f"reproduces the observed decorrelation ({CO.loc['P9','cohA_0_5m']:.2f}→{CO.loc['P9','cohA_80_160m']:.2f} for "
    f"pink σ = 0.18, {CO.loc['G50','cohA_0_5m']:.2f}→{CO.loc['G50','cohA_80_160m']:.2f} for Gaussian *a* = 50 m). "
    f"This independent observable rejects smooth crust for the same reason as the coda level. The observed D1 "
    f"(median {CO.loc['data','D1_median']:+.2f}) is likewise far from the smooth crusts "
    f"({CO.loc['H','D1_median']:+.2f}) and just below the range of the scattering media "
    f"({min(CO.loc[m,'D1_median'] for m in ('P3','P6','P9','E15','E50','G15','G50')):+.2f} to "
    f"{max(CO.loc[m,'D1_median'] for m in ('P3','P6','P9','E15','E50','G15','G50')):+.2f}). D1 is, however, also sensitive "
    f"to intrinsic attenuation, which the FD model omits.")
d.figure("Frequency dependence and inter-event coherence of the S coda",
         "a) D1, coda level at 40–80 Hz minus that at 20–40 Hz (sensor A); bars are interquartile ranges, the black line "
         "and grey band the FORGE median and interquartile range. b, c) D2, median coda coherence versus inter-event "
         "separation for sensors A and B. Black: FORGE; grey: smooth crusts; blue: correlation-length media; orange: "
         "pink crusts.", F / "fig_shape_data.png")
d.h("Can spectral shape be resolved? Equal-σ truth test", 2)
rows = []
for _, r in CF.iterrows():
    rows.append([f"{r.truth_family} ({r.truth})", f2(r.score_pink), f2(r.score_exponential), f2(r.score_Gaussian),
                 r.selected + (" ✓" if r.selected == r.truth_family else " ✗")])
rows.append(["FORGE data", f2(DFm.loc["pink", "mean"]), f2(DFm.loc["exponential", "mean"]),
             f2(DFm.loc["Gaussian", "mean"]), "—"])
d.table("Equal-σ shape-resolution test: family scores for nine synthetic truths and FORGE",
        ["Truth", "Pink", "Exponential (a = 15 m)", "Gaussian (a = 50 m)", "Selected"], rows,
        [4.2, 2.0, 3.2, 3.2, 3.2],
        legend="Family score = mean Φ over candidate realisations with a seed different from the truth; σ = 0.13 "
               "throughout. FORGE row: mean over seeds 11 and 33 "
               f"(pink {f2(DFm.loc['pink','PHI_seed11'])}/{f2(DFm.loc['pink','PHI_seed33'])}, exponential "
               f"{f2(DFm.loc['exponential','PHI_seed11'])}/{f2(DFm.loc['exponential','PHI_seed33'])}, Gaussian "
               f"{f2(DFm.loc['Gaussian','PHI_seed11'])}/{f2(DFm.loc['Gaussian','PHI_seed33'])}).")
d.p(f"Table 3 lists the nine truth tests. The correct family is selected in {n_ok} of 9 cases, against 3 expected "
    f"by chance. The errors are structured. Pink and exponential truths are confused with each other: {n_pe} of 6 "
    f"are assigned to the right member of the pair, which is chance level. However, none of the six is assigned to "
    f"the Gaussian family, and small-scale-rich and Gaussian media are separated in {n_pl} of 9 cases overall. The "
    f"pink–exponential confusion has a simple cause. Built from the same white noise, the pink field and the "
    f"exponential field with *a* = 15 m are {rexp*100:.0f} % correlated over the simulated scales (10–700 m), "
    f"whereas the Gaussian field with *a* = 50 m is {rgau*100:.0f} % correlated with it. Over this range the two are "
    f"almost the same medium, and they would separate only at scales outside the FD band. The D1 and D2 "
    f"diagnostics add no shape information: the combined D1–D2 score selected the Gaussian candidate for every "
    f"truth in the single-realisation test.")
d.p(f"The second realisation changes the reading of the FORGE comparison. With seed 11 the Gaussian medium fits "
    f"best (Φ = {f2(DFm.loc['Gaussian','PHI_seed11'])}), but with seed 33 it scores "
    f"{f2(DFm.loc['Gaussian','PHI_seed33'])}, behind the exponential medium. The family means "
    f"({f2(DFm.loc['pink','mean'])}, {f2(DFm.loc['exponential','mean'])} and {f2(DFm.loc['Gaussian','mean'])} for "
    f"pink, exponential and Gaussian) differ by {f2(DFm['mean'].max()-DFm['mean'].min())}. This is less than the "
    f"median spread between family scores in the truth tests ({f2(fam_spread.median())}) and comparable to the "
    f"realisation scatter of a single family ({f2(abs(DFm.loc['Gaussian','PHI_seed11']-DFm.loc['Gaussian','PHI_seed33']))} "
    f"for the Gaussian medium). The seismic data therefore do not prefer any of the three spectral families.")
d.h("Competing random media and the well log", 2)
d.p(f"At equal σ = 0.13, the correlation-length media fit the FORGE coda statistics as well as or better than the "
    f"pink crust: Φ = {f2(MD.loc['E15','PHI'])} and {f2(MD.loc['E50','PHI'])} (exponential), {f2(MD.loc['G15','PHI'])} "
    f"and {f2(MD.loc['G50','PHI'])} (Gaussian, *a* = 15 and 50 m), against {f2(MD.loc['P6','PHI'])} for P6 and "
    f"{f2(MD.loc[best,'PHI'])} for the best pink crust. The same holds with random mechanisms. In the 20–80 Hz band "
    f"and at 100–600 m range, the coda therefore constrains the strength of heterogeneity but not its spectral shape.")
d.p(f"The shape is constrained by the sonic log (Fig. 12). Over 12–200 m the 56-32 log spectrum is a straight power "
    f"law with no corner. Vertical lines through the pink field reproduce it best (rms difference "
    f"{lv('pink p=1.5'):.2f} in log_{{10}} power for *p* = 1.5; {lv('pink p=1.8'):.2f} for *p* = 1.8). The exponential "
    f"media do so less well ({lv('exponential a=15 m'):.2f} and {lv('exponential a=50 m'):.2f}), because their spectra "
    f"steepen toward *k*^{{-2}} above the corner. The Gaussian media, with their exponential cut-off, are incompatible "
    f"with the log ({lv('Gaussian a=15 m'):.1f} and {lv('Gaussian a=50 m'):.0f}). The log therefore excludes the one "
    f"family that the wavefield can separate from pink noise, and leaves pink and exponential media, which the "
    f"wavefield cannot separate, with a modest preference for pink.")
d.figure("Vertical fluctuation spectra of the 56-32 sonic log and of the model media",
         "Spectra of detrended ln *V* along vertical lines, normalised at 50-m wavelength, over the common 12–200 m "
         "band, for the 56-32 sonic log (black) and the model media. The legend gives the fitted slope and the rms "
         "log_{10} difference from the log.",
         F / "fig_log_vs_media.png")
d.h("Bandwidth of the FORGE records", 2)
d.p(f"The FD analysis above is confined to 20–80 Hz, but the 56-32 records extend much higher. Fig. 13 shows four "
    f"events of magnitude −0.6 to +0.2 at 300–370 m from sensor B. Their spectrograms carry S and coda energy up to "
    f"about 1.5 kHz, and the S-window spectra stand one to three orders of magnitude above pre-event noise from "
    f"about 20 Hz to the anti-alias limit. Over all events (Fig. 14), the median amplitude SNR rises from about 4 at "
    f"5 Hz to several hundred at 300–1000 Hz, and every event exceeds SNR 3 up to about 1.8 kHz "
    f"({int(BW.loc['A','n_events'])} events at sensor A, {int(BW.loc['B','n_events'])} at B). A recorder sampling "
    f"at 200 samples/s (Nyquist 100 Hz) would capture only the low-frequency tail of these S waves. Nearly all of "
    f"the S-wave velocity energy lies above 100 Hz (median {BW.loc['A','pct_S_energy_above100_med']:.1f} % at A, "
    f"{BW.loc['B','pct_S_energy_above100_med']:.1f} % at B). This percentage is inflated, however, by sharp spectral "
    f"peaks near 330, 480 and 920 Hz, which occur in nearly every event. Peak frequencies cluster at these fixed "
    f"values, so the peaks are most likely resonances of the tool or its coupling rather than source or path "
    f"properties.")
d.figure("Example FORGE 2022 microearthquakes recorded at 4 kHz in well 56-32",
         "Four events (rows; magnitude −0.61 to +0.24) at sensor B, 4 kHz. a, d, g, j) Horizontal velocity (> 5 Hz) "
         "with P (blue) and S (red) marked. b, e, h, k) Spectrograms to 2 kHz on a common dB scale relative to each "
         "event's maximum; white dotted and dashed lines mark the Nyquist frequencies of 100 and 200 sample/s "
         "recording, and vertical lines the P and S times. c, f, i, l) S-window spectrum (black) and pre-event noise "
         "(grey); shaded, the FD analysis band 20–160 Hz.", F / "fig_event_spectrograms.png")
d.figure("Signal-to-noise ratio and spectral peak frequency of all FORGE records",
         f"a, b) Median (solid) and interquartile range (dotted) of the S-window amplitude SNR spectrum for sensors A "
         f"and B; vertical lines at 50 and 100 Hz, horizontal line at SNR = 3. c) Peak frequency of the S-window "
         f"velocity spectrum for all events. The clustering at about 330, 480 and 920 Hz indicates instrument or "
         f"coupling resonances. Median share of S-wave velocity energy above 100 Hz: "
         f"{BW.loc['A','pct_S_energy_above100_med']:.1f} % (A) and {BW.loc['B','pct_S_energy_above100_med']:.1f} % (B).",
         F / "fig_bandwidth.png")
d.h("Frequency dependence of the coda over two decades", 2)
d.p(f"Fig. 15 and Table 4 give the coda measures in seven octave bands. Q_{{c}}^{{-1}} falls steadily from "
    f"{cmo('A',25,'Qc_inv_med'):.2f} at 25 Hz to {cmo('A',1304,'Qc_inv_med'):.4f} at 1.3 kHz (sensor A), so the "
    f"energy decay rate varies little with frequency ({CMO[CMO.sensor=='A'].decay_per_s_med.min():.0f}–"
    f"{CMO[CMO.sensor=='A'].decay_per_s_med.max():.0f} s^{{-1}} over all bands). The coda level and the S-pulse width decrease up to 200 Hz and "
    f"rise again above it. The rise coincides with the suspected resonance bands. A high-Q resonance does not change "
    f"within-band energy ratios in the linear limit, but its ringing lengthens the S pulse and raises the coda. We "
    f"therefore interpret only the resonance-free band, 25–200 Hz (S wavelengths of about 17–134 m; 50–200 Hz at "
    f"sensor B, whose 25-Hz band carries the borehole plateau).")
rows = []
for s in ("A", "B"):
    for m, lab in (("QI", "Q_{c}^{-1}"), ("CL", "coda level"), ("SW", "S-pulse width")):
        rows.append([s, lab, ccb(s, m, "band_Hz"),
                     f"{ccb(s,m,'slope'):.2f} [{ccb(s,m,'slope_lo95'):.2f}, {ccb(s,m,'slope_hi95'):.2f}]",
                     f"{ccb(s,m,'dBIC_break'):+.0f}",
                     f"{ccb(s,m,'f_break'):.0f}" if ccb(s, m, 'dBIC_break') < 0 else "—"])
d.table("Frequency dependence of FORGE coda measures in the resonance-free band",
        ["Sensor", "Measure", "Band (Hz)", "Power-law slope [95 %]", "ΔBIC (break − single)", "Break (Hz)"],
        rows, [1.5, 3.0, 2.2, 4.2, 3.0, 2.0],
        legend="Slopes are d log_{10}(measure)/d log_{10} *f* (for the coda level, d log_{10}(E_{coda}/E_{S})/d "
               "log_{10} *f*). Confidence intervals are from 500 bootstrap resamples of events. Negative ΔBIC favours "
               "a break at the listed frequency.")
d.p(f"In this band Q_{{c}}^{{-1}} follows a single power law on both sensors, *f*^{{{ccb('A','QI','slope'):.2f}}} "
    f"[{ccb('A','QI','slope_lo95'):.2f}, {ccb('A','QI','slope_hi95'):.2f}] at A and "
    f"*f*^{{{ccb('B','QI','slope'):.2f}}} [{ccb('B','QI','slope_lo95'):.2f}, {ccb('B','QI','slope_hi95'):.2f}] at B. "
    f"A break is disfavoured (ΔBIC = {ccb('A','QI','dBIC_break'):+.0f} and {ccb('B','QI','dBIC_break'):+.0f}), so "
    f"Q_{{c}}^{{-1}} shows no characteristic frequency over S wavelengths of about 17–134 m. For a medium with a "
    f"correlation length *a*, a change of behaviour would be expected near *ka* ≈ 1. For 25–200 Hz this corresponds "
    f"to *a* ≈ 3–21 m, not to the wavelengths themselves. Whether this test can detect a correlation length is "
    f"examined with simulations below (Table 5); it turns out that it cannot. The coda level, by contrast, flattens "
    f"above 50–100 Hz "
    f"(ΔBIC = {ccb('A','CL','dBIC_break'):.0f} and {ccb('B','CL','dBIC_break'):.0f}). This may reflect the growing "
    f"influence of the resonances toward 200 Hz, and we do not use it for inference.")
d.figure("FORGE S-coda measures in seven octave bands from 25 Hz to 1.3 kHz",
         "Median and interquartile range for sensors A and B of a) the coda level, b) the single-scattering "
         "Q_{c}^{-1} and c) the rms S-pulse width. The dotted line marks *ka* = 1 for *a* = 15 m. Thin grey lines "
         "mark the suspected resonance frequencies, and the shaded band (> 200 Hz) is affected by them.", F / "fig_coda_multioctave.png")
qp, qe, qg5 = qs_fam('pink s=.09-.18'), qs_fam('exponential a=15'), qs_fam('Gaussian a=50')
d.p(f"The FD synthetics resolve the 25 and 50 Hz bands, which gives a direct test of the Q_{{c}}^{{-1}} slope "
    f"between them (Fig. 16). The FORGE slope at sensor A is {qs_d.slope_25_50:.2f} (95 % bootstrap interval "
    f"{qs_d.slope_lo95:.2f} to {qs_d.slope_hi95:.2f}). Smooth crusts give {qs_fam('smooth').min():.2f} to "
    f"{qs_fam('smooth').max():.2f} and are excluded. Six pink realisations give {qp.min():.2f} to {qp.max():.2f} "
    f"(median {qp.median():.2f}), and four exponential (*a* = 15 m) realisations give {qe.min():.2f} to {qe.max():.2f} "
    f"(median {qe.median():.2f}); both bracket the data. The Gaussian medium with *a* = 50 m gives {qg5.min():.2f} to "
    f"{qg5.max():.2f} (median {qg5.median():.2f}), and only {int((qg5.between(qs_d.slope_lo95, qs_d.slope_hi95)).sum())} "
    f"of {len(qg5)} realisations fall inside the data interval. The Gaussian medium with *a* = 15 m is also steep "
    f"({qs_fam('Gaussian a=15').iloc[0]:.2f}). The slope is therefore sensitive to heterogeneity at scales of about "
    f"15 m and below rather than to spectral form as such. It favours media that contain such small-scale structure, "
    f"which is the property the log also shows. In the 25-Hz band the coda window spans only one to three periods, and "
    f"even the homogeneous crust gives a non-zero Q_{{c}}^{{-1}}. Data and synthetics are measured identically, so "
    f"the comparison is fair, but the absolute Q_{{c}}^{{-1}} values at 25 Hz should not be read as medium "
    f"properties.")
d.figure("Q_{c}^{-1} frequency slope between 25 and 50 Hz: FORGE versus simulated realisations",
         "Q_{c}^{-1} frequency slope log_{2}[Q_{c}^{-1}(50 Hz)/Q_{c}^{-1}(25 Hz)] at sensor A. Each dot is one "
         "simulated medium (σ = 0.13 unless stated; pink σ = 0.09–0.18). Black line and grey band: FORGE with its "
         "95 % bootstrap interval.", F / "fig_qc_slope.png")
d.h("Coda at 100 Hz and the correlation-length series", 2)
rows = []
lab = {"H": "homogeneous", "P6": "pink (seed 11)", "R6": "pink (seed 33)", "E15": "exponential, a = 15 m",
       "E15b": "exponential, a = 15 m (seed 33)", "G50": "Gaussian, a = 50 m", "G50b": "Gaussian, a = 50 m (seed 33)",
       "EX5": "exponential, a = 5 m", "EX9": "exponential, a = 9 m", "EX17": "exponential, a = 17 m",
       "EX34": "exponential, a = 34 m", "EX67": "exponential, a = 67 m", "EX134": "exponential, a = 134 m"}
fm = lambda v: "—" if pd.isna(v) else f"{v:+.2f}"
for m in lab:
    rows.append([lab[m], fm(bt(m, 'slope_25_50')), fm(bt(m, 'slope_50_100')), fm(bt(m, 'curvature')),
                 f"{bt(m, 'CL100_B'):.2f}"])
rows.append(["FORGE data", f"{bt('data','slope_25_50'):+.2f} [{bt('data','s1_lo'):+.2f}, {bt('data','s1_hi'):+.2f}]",
             f"{bt('data','slope_50_100'):+.2f} [{bt('data','s2_lo'):+.2f}, {bt('data','s2_hi'):+.2f}]",
             f"{bt('data','curvature'):+.2f} [{bt('data','curv_lo'):+.2f}, {bt('data','curv_hi'):+.2f}]",
             f"{bt('data','CL100_B'):.2f}"])
d.table("Q_{c}^{-1} slopes, their curvature and the 100-Hz coda level",
        ["Medium (σ = 0.13)", "Slope 25–50 Hz (A)", "Slope 50–100 Hz (B)", "Curvature", "Coda level 100 Hz (B)"],
        rows, [5.2, 3.0, 3.0, 3.0, 2.6],
        legend="Slopes are log_{2} ratios of median Q_{c}^{-1} between octave bands; curvature = slope(50–100) − "
               "slope(25–50). Brackets: 95 % bootstrap intervals over events. Coda level = log_{10}(E_{coda}/E_{S}) at "
               "100 Hz. The 25–50 Hz slope uses the 5-m grid, the 50–100 Hz slope and 100-Hz level the 2.5-m grid; "
               "a = 5 m is not resolved on the 5-m grid.")
cv = BT.drop(index="data").curvature.dropna()
d.p(f"Table 5 and Fig. 17 give the results of the correlation-length series. The Q_{{c}}^{{-1}} curvature does not "
    f"identify a correlation length. The simulated media, including the pink, Gaussian and homogeneous crusts, give "
    f"curvatures of {cv.min():+.2f} to {cv.max():+.2f}, positive in all but {['none','one','two','three'][int((cv < 0).sum())]}, with no systematic "
    f"dependence on *a*, while "
    f"FORGE gives {bt('data','curvature'):+.2f} [{bt('data','curv_lo'):+.2f}, {bt('data','curv_hi'):+.2f}]. The "
    f"single power law of the observed Q_{{c}}^{{-1}} is therefore not evidence for scale-free heterogeneity. The "
    f"common positive curvature of the models probably reflects the short coda window in the 25-Hz band rather than "
    f"the medium. The coda level at 100 Hz is far more selective. It falls steadily with correlation length in the "
    f"exponential series ({bt('EX9','CL100_B'):.2f} for *a* = 9 m to {bt('EX134','CL100_B'):.2f} for *a* = 134 m), "
    f"because large-*a* media lack the small-scale structure that scatters 30-m wavelengths. The Gaussian media with "
    f"*a* = 50 m ({bt('G50','CL100_B'):.2f}, {bt('G50b','CL100_B'):.2f}) and the homogeneous crust "
    f"({bt('H','CL100_B'):.2f}) fall 0.6–1.6 log units below FORGE ({bt('data','CL100_B'):.2f}). The pink realisations "
    f"({bt('P6','CL100_B'):.2f}, {bt('R6','CL100_B'):.2f}) and the exponential media with *a* ≳ 30 m lie within about "
    f"0.4. The Gaussian medium that matched the 20–80 Hz coda statistics as well as pink noise is thus rejected once "
    f"the analysis reaches 100 Hz. The observed 50–100 Hz slope ({bt('data','slope_50_100'):+.2f}) is steeper than in "
    f"every simulated medium except the smallest-scale exponential ones (*a* = 5 and 9 m: {bt('EX5','slope_50_100'):+.2f}, "
    f"{bt('EX9','slope_50_100'):+.2f}); the pink realisations ({bt('P6','slope_50_100'):+.2f}, "
    f"{bt('R6','slope_50_100'):+.2f}) lie at the edge of its interval. Among the simulated media, only pink noise comes "
    f"close to both the 100-Hz coda level and the 50–100 Hz slope. Media rich in very small scales match the slope "
    f"but give too much 100-Hz coda at σ = 0.13, and media with large *a* match the level but not the slope. With "
    f"one realisation per correlation length and σ fixed, this is suggestive rather than decisive.")
d.figure("Q_{c}^{-1} frequency dependence and 100-Hz coda level versus correlation length",
         "a) Q_{c}^{-1} normalised at 50 Hz, from the 25–50 Hz slope (sensor A, 5-m grid, circles) and the 50–100 Hz "
         "slope (sensor B, 2.5-m grid, squares); black: FORGE; blue: exponential series (darker = larger *a*); orange: "
         "pink; green: Gaussian *a* = 50 m; grey: homogeneous. b) Curvature, slope(50–100 Hz) − slope(25–50 Hz), "
         "versus correlation length; black line and grey band: FORGE with 95 % bootstrap interval; pink (orange) and "
         "homogeneous (grey) media are plotted at the right and left edges. c) Coda level at 100 Hz (sensor B) versus "
         "correlation length; black line: FORGE; symbols as in b).", F / "fig_break_test.png")

# ================================================================= DISCUSSION
d.h("Discussion", 1)
d.h("Why deterministic waveform fitting favours smooth crust", 2)
d.p("With two sensors, the realisation of a scale-free random medium cannot be recovered, and any assumed "
    "realisation adds coherent scattered energy at wrong times. Deterministic waveform fitting, whether by "
    "moment-tensor or full-waveform inversion, therefore penalises heterogeneity: it favours the smoothest model "
    "that fits the direct waves and assigns the coda to noise. Our truth test shows that this preference persists "
    "when the truth is a pink-noise crust. Conventional practice is thus biased toward the smooth models it "
    "assumes [@virieux2009]. A scale-free crust has to be tested through wavefield statistics.")
d.h("Heterogeneity strength: log versus wavefield", 2)
d.p(f"The FORGE coda requires σ ≈ 0.09–0.18 at 5-m scale (per-event test near 0.09–0.13, statistical optimum near "
    f"0.18), about {0.09/K.log_sig6:.0f}–{0.18/K.log_sig6:.0f} times the "
    f"sonic-log value. Several effects may contribute. (i) The log samples one vertical line in intact rock outside "
    f"the stimulated volume, whereas the MEQ paths cross the fracture volume being pressurised by injection, where "
    f"fluid-filled fractures lower *V*_{{S}} strongly; in the GFI framework this is the poro-connectivity "
    f"enhancement of stimulation. (ii) Fluid-filled fractures scatter S waves more strongly than an isotropic "
    f"lognormal perturbation with coupled *V*_{{P}} and *V*_{{S}}. (iii) Near-sensor effects such as cement "
    f"coupling and borehole modes add coda that our model lacks; the sensor-B 20–40 Hz plateau is one example. "
    f"(iv) Intrinsic attenuation is absent from the FD model, although it would *shorten* the coda and so cannot "
    f"explain the excess. The inferred σ should therefore be read as an effective scattering strength of the "
    f"stimulated volume, with the log value as a reference for intact ambient crust.")
d.h("What the data do and do not constrain", 2)
d.p("Three results should be kept apart. First, every smooth crust (homogeneous, layered, VTI) fails, on three "
    "independent observables. It under-predicts the observed S coda by more than an order of magnitude, event by "
    "event; it keeps the coda of neighbouring events coherent where the data decorrelate; and it fails the "
    "statistical test. VTI anisotropy changes the direct waves but creates neither coda nor decorrelation. Second, "
    "the coda requires strong heterogeneity (σ of roughly 0.1–0.2 at 5-m scale) but does not establish its spectral "
    "shape. In repeated truth tests the wavefield separates small-scale-rich media from Gaussian media reliably, but "
    "cannot separate pink from exponential media, which are nearly the same field over 10–700 m. On FORGE, with two "
    "realisations per family, it prefers no family in the 20–80 Hz coda statistics. Higher frequencies are more "
    "selective. The 25–50 Hz Q_{c}^{-1} slope excludes smooth crusts and most large-scale Gaussian realisations, "
    "and the 100-Hz coda level rejects the large-scale Gaussian medium outright while matching pink noise. Third, the "
    "borehole data supply the shape. The sonic log is "
    "scale-free over 12–200 m, which is consistent with pink noise, marginally with an exponential medium and not with "
    "a Gaussian one. The two lines of evidence are complementary rather than in conflict: the seismic data fix the "
    "strength of heterogeneity in the stimulated volume, and the log fixes its scaling in the intact rock. Taken "
    "together they support a strongly scattering crust with power-law heterogeneity, as assumed in GFI. The MEQ "
    "wavefield alone is not independent proof of the pink 1/*k* spectrum. Testing the shape seismically would "
    "require wider-band records, up to the kHz corner frequencies of FORGE MEQs, and a finer FD grid, because the "
    "pink and exponential media differ mainly below the present 10-m resolution.")
d.h("What kind of heterogeneity makes the coda?", 2)
d.p("That coda waves are scattered by crustal heterogeneity is long established [@aki1975; @sato2012], and "
    "power-law heterogeneity spectra have been inferred from well logs [@holliger1996; @leary1997]. What has "
    "remained open is which heterogeneity produces the coda of a given reservoir. Stochastic models are often "
    "chosen for convenience and fitted with a free correlation length. The FORGE data allow a more specific "
    "answer at the reservoir scale. The heterogeneity must be strong (σ ≈ 0.1–0.2 at 5-m scale), and it must "
    "retain enough small-scale structure to sustain the coda at 100 Hz. That requirement rejects the smooth crusts "
    "and the large-scale Gaussian medium. Of the media tested, pink noise alone comes close to both the 100-Hz coda "
    "level and the 50–100 Hz Q_{c}^{-1} slope. The single power law of Q_{c}^{-1} over 25–200 Hz is not, by itself, "
    "evidence of scale-free heterogeneity, because media with correlation lengths of 9–134 m reproduce it equally "
    "well in our simulations. The borehole data show a scale-free spectrum over 12–200 m. A lognormal pink-noise "
    "field, the same statistical object that describes poro-permeability in the GFI framework, is consistent with "
    "all of these observations. If this holds more widely, the coda is not a nuisance "
    "tail but a direct observable of the flow-controlling heterogeneity. This link can be tested only with "
    "broadband near-source records such as those at FORGE. At 200 samples/s the band above 100 Hz, which carries "
    "most of the S-wave energy of these events, is not recorded at all. Our own results show why the band matters: "
    "at 20–80 Hz a Gaussian medium fits the coda as well as pink noise, and it fails only at 100 Hz. The present "
    "data do not, however, exclude an exponential medium, which differs from pink noise mainly below 10 m.")
d.h("Implications for microseismic flow imaging", 2)
d.p("A strongly heterogeneous crust with σ of several percent is invisible to travel-time location (residuals of ~1 ms) "
    "but dominates the MEQ coda. For EGS monitoring this means that (a) travel-time residuals of ~1 ms "
    "are compatible with smooth-model locations, whereas (b) the scattered wavefield carries information on the heterogeneity "
    "that controls flow. Coda and envelope statistics, rather than ray-theoretical travel times, are the natural "
    "observables for GFI/seismic emission tomography [@leary_abg].")
d.h("Limitations", 2)
d.bullets([
    "Only two downhole sensors, with unoriented horizontals and questionable vertical channels, are available for 2022.",
    "The FD band (20–80 Hz) covers only the low end of the FORGE MEQ spectrum; the data-model comparison of the "
    "Q_{c}^{-1} slope uses two octave bands only, and the 25-Hz band is short relative to the coda window.",
    "Spectral peaks near 330, 480 and 920 Hz, attributed to the tool or its coupling, prevent interpretation of the "
    "coda above about 200 Hz.",
    "The model grid in (σ, *p*) is coarse and *p* is weakly resolved; intrinsic *Q* and fracture-specific (non-lognormal) "
    "scatterers were not included.",
    "Correlation-length competitors were simulated only at σ = 0.13, with at most two candidate and three truth "
    "realisations per family and a single realisation per correlation length in the *a* = 5–134 m series; a full "
    "(σ, *a*) search was not performed, so the trade-off between σ and *a* in the 100-Hz coda level is not resolved.",
    "The 25–50 Hz and 50–100 Hz slopes come from different sensors and grids, and the 2.5-m-grid synthetics are "
    "noise-free.",
    "Catalogue hypocentres were computed by the operator with a smooth velocity model, so the travel-time residual "
    "statistics mix structure with location error.",
])

# ================================================================= CONCLUSIONS
d.h("Conclusions", 1)
d.bullets([
    f"Downhole MEQ records of the 2022 FORGE stimulation were matched to the catalogue and the 56-32 sensors "
    f"self-calibrated (depths {K.sensorA_depth_m:.0f} and {K.sensorB_depth_m:.0f} m).",
    "Reciprocity-based 3-D elastic FD gives full-waveform synthetics for hundreds of paths per model, validated "
    "to 0.995 correlation in heterogeneous VTI media.",
    "Deterministic waveform inversion favours the smoothest crust even for a pink-noise truth and cannot "
    "identify the medium class.",
    "Statistical waveform inversion recovers a pink and a layered-VTI synthetic truth without bias.",
    "For FORGE, homogeneous, layered and VTI crusts under-predict the observed S coda of 150 events by a factor of "
    "about 10–30, keep the coda of neighbouring events coherent where the data decorrelate, and rank last in the "
    "statistical inversion. Smooth crustal models are rejected.",
    f"The coda requires strong heterogeneity: the pink optimum is σ ≈ {sig_best:.2f}, bracketed by runs at σ = 0.13 "
    f"and 0.25, and the per-event test gives σ ≈ 0.09–0.13. This is several times the sonic-log value, which "
    f"suggests stimulation-enhanced scattering.",
    f"At equal σ, nine synthetic truth tests show that the wavefield separates small-scale-rich from Gaussian media "
    f"({n_pl} of 9) but not pink from exponential media ({n_pe} of 6), which are {rexp*100:.0f} % correlated on the "
    f"simulated scales. At 20–80 Hz the FORGE data prefer no spectral family; at 100 Hz they reject the large-scale "
    f"Gaussian medium. The scale-free 56-32 log spectrum favours pink "
    f"noise and excludes Gaussian media. The seismic data establish heterogeneity strength; the pink spectral form "
    f"rests on the borehole data.",
    f"The 4-kHz records carry S-wave energy above noise to about 1.8 kHz. Over 25–200 Hz the coda attenuation "
    f"follows a single power law, Q_{{c}}^{{-1}} ∝ *f*^{{{ccb('A','QI','slope'):.1f}}}. Simulations with correlation "
    f"lengths of 5–134 m show that this is not diagnostic of scale-free heterogeneity.",
    f"Raising the simulated band to 100 Hz (2.5-m grid) separates the media that 20–80 Hz could not. The 100-Hz coda "
    f"level rejects the large-scale Gaussian medium and the smooth crust by 0.6–1.6 log units. Pink noise comes "
    f"closest to both the 100-Hz coda level and the 50–100 Hz Q_{{c}}^{{-1}} slope; with one realisation per medium "
    f"this is suggestive rather than decisive. Everything above 100 Hz is inaccessible at 200 samples/s.",
])

# ================================================================= BACK MATTER
d.h("Availability of data and materials", 1)
d.p("The GES 2022 catalogue is public [@dyer2022]. The 56-32 SEG-2 records were provided to the authors; their "
    "public release is being confirmed. All MATLAB code (FD solver, reciprocity, features, inversion) and derived "
    "results are in the Paper2_PinkNoise_WaveformInversion repository [URL TO BE ADDED].")
d.h("Competing interests", 1)
d.p("[TO BE COMPLETED BY AUTHORS]")
d.h("Funding", 1)
d.p("[TO BE COMPLETED BY AUTHORS]")
d.h("Authors' contributions", 1)
d.p("[TO BE COMPLETED BY AUTHORS]")
d.h("Acknowledgements", 1)
d.p("We thank Heiner Igel for the acoustic finite-difference MATLAB codes from which the elastic solver used here "
    "was developed. [FURTHER ACKNOWLEDGEMENTS TO BE COMPLETED BY AUTHORS]")

REFS = {
    "leary1997": "Leary P. Rock as a critical-point system and the inherent implausibility of reliable earthquake prediction. Geophys J Int. 1997;131:451–66.",
    "leary2002": "Leary PC, Al-Kindy F. Power-law scaling of spatially correlated porosity and log(permeability) sequences from north-central North Sea Brae oilfield well core. Geophys J Int. 2002;148:426–42.",
    "holliger1996": "Holliger K. Upper-crustal seismic velocity heterogeneity as derived from a variety of P-wave sonic logs. Geophys J Int. 1996;125:813–29.",
    "leary_abg": "Leary P. The αβγ of EGS crustal heat extraction. J Energy Power Technol. 2026; in press.",
    "leary2026": "Leary P. EGS sustainability: deconstructing UtahForge engineered geothermal system flow data. Sustainability. 2026;18:5308. https://doi.org/10.3390/su18115308.",
    "suryanto2026": "Suryanto W, Leary P, Saunders G, Onto, Fleure T, Pramono B, et al. Data-constrained multiscale descriptors of heterogeneous flow and microseismicity in enhanced geothermal systems. Geotherm Energy. 2026; submitted.",
    "aki1975": "Aki K, Chouet B. Origin of coda waves: source, attenuation, and scattering effects. J Geophys Res. 1975;80:3322–42.",
    "frankel1986": "Frankel A, Clayton RW. Finite difference simulations of seismic scattering: implications for the propagation of short-period seismic waves in the crust and models of crustal heterogeneity. J Geophys Res. 1986;91:6465–89.",
    "sato2012": "Sato H, Fehler MC, Maeda T. Seismic wave propagation and scattering in the heterogeneous Earth. 2nd ed. Berlin: Springer; 2012.",
    "thomsen1986": "Thomsen L. Weak elastic anisotropy. Geophysics. 1986;51:1954–66.",
    "dyer2022": "Dyer B, Karvounis D. Utah FORGE: updated seismic event catalogue from the April 2022 stimulation of well 16A(78)-32. DOE Data Explorer. 2022. https://www.osti.gov/dataexplorer/biblio/dataset/1908927. Accessed 26 Sep 2026.",
    "maeda1985": "Maeda N. A method for reading and checking phase times in autoprocessing system of seismic wave data. Zisin. 1985;38:365–79.",
    "virieux1986": "Virieux J. P-SV wave propagation in heterogeneous media: velocity-stress finite-difference method. Geophysics. 1986;51:889–901.",
    "levander1988": "Levander AR. Fourth-order finite-difference P-SV seismograms. Geophysics. 1988;53:1425–36.",
    "igel1995": "Igel H, Mora P, Riollet B. Anisotropic wave propagation through finite-difference grids. Geophysics. 1995;60:1203–16.",
    "igel2016": "Igel H. Computational seismology: a practical introduction. Oxford: Oxford University Press; 2016.",
    "cerjan1985": "Cerjan C, Kosloff D, Kosloff R, Reshef M. A nonreflecting boundary condition for discrete acoustic and elastic wave equations. Geophysics. 1985;50:705–8.",
    "zhao2006": "Zhao L, Chen P, Jordan TH. Strain Green's tensors, reciprocity, and their applications to seismic source and structure studies. Bull Seismol Soc Am. 2006;96:1753–63.",
    "virieux2009": "Virieux J, Operto S. An overview of full-waveform inversion in exploration geophysics. Geophysics. 2009;74:WCC1–26.",
}
d.h("References", 1)
missing = [k for k in d.cite_order if k not in REFS]
assert not missing, missing
for i, k in enumerate(d.cite_order, 1):
    d.p(f"{i}. {REFS[k]}", align="left")
print("unused:", [k for k in REFS if k not in d.cite_order])
d.save(P / "manuscript" / ("Paper2_manuscript_review.docx" if d.embed else "Paper2_manuscript.docx"))
