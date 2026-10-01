"""S18  Uncertainty of the 56-32 sonic-log fluctuation spectrum (reviewer point 16).

Same log, interval and conversion as s04_log_statistics.m (SPHI -> V with DTma 55.5, DTf 189 us/ft,
1950-2780 m, calliper < 10.5 in).  Varies: detrending (mean, linear, quadratic), PSD estimator
(periodogram, Welch, multitaper), sub-interval (thirds), and fits power-law, exponential-ACF
(Lorentzian) and Gaussian-ACF spectra over 12-200 m and 1-200 m, compared by AIC on log-binned PSD.
Outputs: results/table_log_slope_variants.csv, results/table_log_spectral_models.csv,
         figs/fig_log_spectrum_uncertainty.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import welch, periodogram
from scipy.signal.windows import dpss
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LAS = ROOT.parent / "56-32" / "A.1008700.02.01_UniversityUtah_Forge-56-32_ThruBit_Main.las.dat"

L = np.loadtxt(LAS)
L[L == -999.25] = np.nan
z = L[:, 0] * 0.3048
sphi, cal = L[:, 13], L[:, 3]
v = 304800.0 / (55.5 + sphi * (189 - 55.5))
sel = (z > 1950) & (z < 2780) & np.isfinite(v) & (cal < 10.5)
z, lv = z[sel], np.log(v[sel])
dz = np.median(np.diff(z))


def detrend(y, x, order):
    if order < 0:
        return y - y.mean()
    return y - np.polyval(np.polyfit(x, y, order), x)


def psd(y, method):
    y = y - y.mean()
    if method == "periodogram":
        f, p = periodogram(y, fs=1 / dz, window="boxcar", detrend=False)
    elif method == "welch":
        f, p = welch(y, fs=1 / dz, nperseg=min(len(y), 2048), detrend=False)
    else:  # multitaper, NW = 4, K = 7
        n = len(y); tap = dpss(n, 4, 7)
        P = np.abs(np.fft.rfft(tap * y[None, :], axis=1)) ** 2
        f = np.fft.rfftfreq(n, dz); p = P.mean(0) * dz
        return f[1:], p[1:], P[:, 1:] * dz
    return f[1:], p[1:], None


def slope(f, p, lo, hi):
    m = (f >= 1 / hi) & (f <= 1 / lo)
    c = np.polyfit(np.log(f[m]), np.log(p[m]), 1)
    return -c[0]


rows = []
for det, order in (("mean", -1), ("linear", 1), ("quadratic", 2)):
    y = detrend(lv, z, order)
    for meth in ("periodogram", "welch", "multitaper"):
        f, p, _ = psd(y, meth)
        rows.append(dict(detrend=det, psd=meth, segment="full", beta_12_200=slope(f, p, 12, 200),
                         beta_1_200=slope(f, p, 1, 200)))
    # thirds
    edges = np.linspace(z.min(), z.max(), 4)
    for s in range(3):
        m = (z >= edges[s]) & (z < edges[s + 1])
        f, p, _ = psd(detrend(lv[m], z[m], order), "multitaper")
        rows.append(dict(detrend=det, psd="multitaper", segment=f"third {s + 1} ({edges[s]:.0f}-{edges[s + 1]:.0f} m)",
                         beta_12_200=slope(f, p, 12, min(200, (edges[1] - edges[0]) / 2)), beta_1_200=slope(f, p, 1, 100)))
T = pd.DataFrame(rows)
# jackknife over tapers for the reference estimate (linear detrend, multitaper)
y = detrend(lv, z, 1)
f, p, P = psd(y, "multitaper")
jk = [slope(f, np.delete(P, k, 0).mean(0), 12, 200) for k in range(P.shape[0])]
K = P.shape[0]
se_jk = np.sqrt((K - 1) / K * np.sum((np.array(jk) - np.mean(jk)) ** 2))
T.to_csv(ROOT / "results" / "table_log_slope_variants.csv", index=False)
print(T.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
print(f"reference beta(12-200 m) = {slope(f, p, 12, 200):.2f} +- {se_jk:.2f} (taper jackknife)")

# ---------------- spectral-model comparison on log-binned multitaper PSD ----------------
def logbin(f, p, lo, hi, nb=20):
    e = np.logspace(np.log10(1 / hi), np.log10(1 / lo), nb + 1)
    fc, pc = [], []
    for a, b in zip(e[:-1], e[1:]):
        m = (f >= a) & (f < b)
        if m.sum() >= 2:
            fc.append(np.exp(np.mean(np.log(f[m])))); pc.append(np.exp(np.mean(np.log(p[m]))))
    return np.array(fc), np.array(pc)


models = {
    "power law": (lambda q, k: q[0] - q[1] * np.log(k), [0.0, 1.0]),
    "exponential ACF (Lorentzian)": (lambda q, k: q[0] - np.log1p((2 * np.pi * k * np.exp(q[1])) ** 2), [0.0, np.log(20.0)]),
    "Gaussian ACF": (lambda q, k: q[0] - (np.pi * k * np.exp(q[1])) ** 2, [0.0, np.log(20.0)]),
    "power law with corner": (lambda q, k: q[0] - q[1] / 2 * np.log1p((2 * np.pi * k * np.exp(q[2])) ** 2), [0.0, 1.0, np.log(50.0)]),
}
mrows = []
plt.rcParams.update({"font.size": 8, "font.family": "Arial", "axes.linewidth": 0.6})
fig, ax = plt.subplots(1, 2, figsize=(17.4 / 2.54, 7.0 / 2.54))
for ib, (lo, hi) in enumerate(((12, 200), (1, 200))):
    fc, pc = logbin(f, p, lo, hi)
    yv = np.log(pc); n = len(yv)
    ax[ib].loglog(f, p, color="0.8", lw=0.5); ax[ib].loglog(fc, pc, "ko", ms=2.5, label="binned log PSD")
    for name, (fun, q0) in models.items():
        r = least_squares(lambda q: fun(q, fc) - yv, q0)
        rss = np.sum(r.fun ** 2); kpar = len(q0)
        aic = n * np.log(rss / n) + 2 * kpar
        extra = ""
        if name != "power law" and name != "Gaussian ACF":
            extra = f"a = {np.exp(r.x[-1]):.1f} m"
        if name == "Gaussian ACF":
            extra = f"a = {np.exp(r.x[1]):.1f} m"
        mrows.append(dict(band=f"{lo}-{hi} m", model=name, n_bins=n, rms_log=np.sqrt(rss / n), AIC=aic,
                          params=np.array2string(r.x, precision=3), note=extra))
        kk = np.logspace(np.log10(1 / hi), np.log10(1 / lo), 100)
        ax[ib].loglog(kk, np.exp(fun(r.x, kk)), lw=1.0, ls=("--" if name == "power law with corner" else "-"), label=f"{name} (AIC {aic:.1f})")
    ax[ib].set_xlim(1 / hi * 0.8, 1 / lo * 1.2); ax[ib].set_xlabel("vertical wavenumber (cycles/m)")
    ax[ib].set_ylabel("PSD of ln V"); ax[ib].legend(fontsize=6, loc="lower left", frameon=False); ax[ib].grid(True, which="both", alpha=0.3, lw=0.4)
    ax[ib].text(0.97, 0.95, "ab"[ib] + ")", transform=ax[ib].transAxes, ha="right", va="top", fontweight="bold", fontsize=9,
                bbox=dict(facecolor="w", edgecolor="none", pad=1))
M = pd.DataFrame(mrows)
for b in M.band.unique():
    M.loc[M.band == b, "dAIC"] = M.loc[M.band == b, "AIC"] - M.loc[M.band == b, "AIC"].min()
M.to_csv(ROOT / "results" / "table_log_spectral_models.csv", index=False)
print(M.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
fig.tight_layout(); fig.savefig(ROOT / "figs" / "fig_log_spectrum_uncertainty.png", dpi=300)
fig.savefig(ROOT / "figs" / "fig_log_spectrum_uncertainty.pdf")
pd.DataFrame([dict(beta_ref=slope(f, p, 12, 200), beta_se_jackknife=se_jk,
                   beta_min=T.beta_12_200.min(), beta_max=T.beta_12_200.max())]).to_csv(
    ROOT / "results" / "table_log_slope_summary.csv", index=False)
