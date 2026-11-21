"""Illustrative narrowband patterns of the proposed NavIC L5 8-element (7+1) CRPA.

Polar sky plots (azimuth around the circle, zenith at the centre) at NavIC L5 for
(a) a single antenna, (b) power inversion (PI), (c) one MVDR beam steered to PRN 6 and
shared by all channels, and (d) the per-satellite MVDR beam of channel A3, one of the
per-satellite beams adopted in the paper. Writes the figure used in the manuscript and
the null depths, satellite gains and output SINR (crpa8_metrics.txt), including the
common-beam LCMV with all seven constraints and scenario J6.

Constellation: seven satellites, i.e., the four NavIC satellites at the start of the
deep-tracking evaluation (Nashik, 1 January 2025) and three artificial satellites at the
vacant nominal NavIC slots (navic_constellation.py).

Assumptions (stated in the paper): isotropic elements, ideal calibration,
no mutual coupling, narrowband model at the L5 carrier, spatially white noise.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from navic_constellation import REAL, ARTIFICIAL

C = 299792458.0
F_L5 = 1176.45e6
LAM = C / F_L5
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.environ.get("FIG_OUT") or os.path.join(REPO, "paper", "R1")   # figure PDFs/PNGs
os.makedirs(OUT, exist_ok=True)

# 7 elements on a ring of radius lambda_L5/2 plus one centre (reference) element
R = C / F_L5 / 2
ring = np.array([[R * np.cos(2 * np.pi * k / 7), R * np.sin(2 * np.pi * k / 7), 0.0] for k in range(7)])
POS = np.vstack([[0.0, 0.0, 0.0], ring])  # element 0 = centre, ENU metres
M = POS.shape[0]


def unit(az_deg, el_deg):
    az, el = np.radians(az_deg), np.radians(el_deg)
    return np.stack([np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)], axis=-1)


def steer(az_deg, el_deg):
    return np.exp(1j * 2 * np.pi / LAM * (unit(az_deg, el_deg) @ POS.T))


# Seven-satellite constellation (azimuth, elevation) in degrees
SATS = {**REAL, **{k.split()[0]: v for k, v in ARTIFICIAL.items()}}

# Three broadband jammers, JNR = 50 dB each (per element, over the receiver bandwidth)
JAM = [(40.0, 5.0), (165.0, 12.0), (290.0, 8.0)]
JAM_J6 = [JAM[0], (SATS["PRN 10"][0], SATS["PRN 10"][1] - 10.0), JAM[2]]  # jammer 10 deg below PRN 10
JNR_DB = 50.0
SNR_DB = -30.0  # per-element pre-correlation SNR of one NavIC L5 SPS signal

e0 = np.zeros(M); e0[0] = 1.0


def covariance(jams):
    return np.eye(M, dtype=complex) + sum(10 ** (JNR_DB / 10) * np.outer(steer(*j), steer(*j).conj()) for j in jams)


def weights(Rm):
    Ri = np.linalg.inv(Rm)
    w = {"PI": Ri @ e0 / (e0 @ Ri @ e0)}
    a6 = steer(*SATS["PRN 6"])
    w["MVDR common (PRN 6)"] = Ri @ a6 / (a6.conj() @ Ri @ a6)
    A = np.stack([steer(*v) for v in SATS.values()], axis=1)
    w["LCMV common (7 constraints)"] = Ri @ A @ np.linalg.solve(A.conj().T @ Ri @ A, np.ones(A.shape[1]))
    for n, v in SATS.items():
        a = steer(*v)
        w[f"MVDR per-satellite {n}"] = Ri @ a / (a.conj() @ Ri @ a)
    return w


def sinr_db(w, a, Rm):
    return 10 * np.log10(10 ** (SNR_DB / 10) * abs(w.conj() @ a) ** 2 / np.real(w.conj() @ Rm @ w))


Rm = covariance(JAM)
W = weights(Rm)

az = np.linspace(0, 360, 721)
el = np.linspace(0, 90, 181)
AZ, EL = np.meshgrid(az, el)
S = steer(AZ, EL)

PANELS = [("(a) Single antenna", e0.astype(complex)), ("(b) Power inversion", W["PI"]),
          ("(c) Shared MVDR beam (PRN 6)", W["MVDR common (PRN 6)"]),
          ("(d) Beam of channel A3", W["MVDR per-satellite A3"])]

fig, axs = plt.subplots(1, 4, figsize=(7.16, 2.35), subplot_kw={"projection": "polar"}, constrained_layout=True)
for ax, (title, w) in zip(axs, PANELS):
    G = np.abs(S @ w.conj()) ** 2
    G_db = np.clip(10 * np.log10(G / G.max() + 1e-12), -60, 0)
    cs = ax.contourf(np.radians(AZ), 90 - EL, G_db, levels=np.arange(-60, 1, 3), cmap="viridis")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_ylim(0, 90)
    ax.set_yticks([30, 60]); ax.set_yticklabels([])
    ax.set_xticks(np.radians([0, 90, 180, 270])); ax.set_xticklabels(["N", "E", "S", "W"], fontsize=6)
    ax.tick_params(pad=-2)
    ax.grid(color="w", alpha=0.4, lw=0.5)
    for k, (ja, je) in enumerate(JAM):
        ax.plot(np.radians(ja), 90 - je, "rx", ms=6, mew=1.8, label="Jammer" if k == 0 else None)
    for k, (n, (sa, se)) in enumerate(SATS.items()):
        real = n.startswith("PRN")
        ax.plot(np.radians(sa), 90 - se, "^" if real else "s", ms=5, mfc="w" if real else "c", mec="k",
                label=("NavIC satellite" if n == "PRN 2" else ("Artificial satellite" if n == "A1" else None)))
        ax.annotate(n.split()[-1], (np.radians(sa), 90 - se), xytext=(3, 3), textcoords="offset points",
                    fontsize=5, color="k", bbox=dict(boxstyle="round,pad=0.1", fc="w", ec="none", alpha=0.8))
    ax.set_title(title, fontsize=7, pad=6)
fig.legend(*axs[0].get_legend_handles_labels(), fontsize=6, loc="lower left", ncol=3, framealpha=0.8)
cb = fig.colorbar(cs, ax=axs, shrink=0.8, pad=0.01)
cb.set_label("Normalized gain (dB)", fontsize=6); cb.ax.tick_params(labelsize=6)
fig.savefig(os.path.join(OUT, "fig_crpa8_pattern.pdf"))
fig.savefig(os.path.join(OUT, "fig_crpa8_pattern.png"), dpi=300)


def report(f, tag, Rm, W, jams):
    f.write(f"== {tag}\n")
    for name, w in W.items():
        if name.startswith("MVDR per-satellite"):
            continue
        g = (np.abs(S @ w.conj()) ** 2).max()
        nulls = [round(float(10 * np.log10(abs(w.conj() @ steer(*j)) ** 2 / g)), 1) for j in jams]
        gs = {n: round(float(10 * np.log10(abs(w.conj() @ steer(*v)) ** 2 / g)), 1) for n, v in SATS.items()}
        sinr = {n: round(float(sinr_db(w, steer(*v), Rm)), 1) for n, v in SATS.items()}
        f.write(f"{name}: null depths (dB) = {nulls}; satellite gains (dB) = {gs}; output SINR (dB) = {sinr}\n")
    per = {n: round(float(sinr_db(W[f'MVDR per-satellite {n}'], steer(*v), Rm)), 1) for n, v in SATS.items()}
    f.write(f"Per-satellite MVDR beams: output SINR (dB) = {per}\n")


with open(os.path.join(HERE, "crpa8_metrics.txt"), "w") as f:
    f.write(f"lambda_L5 = {LAM:.4f} m, ring radius = {R:.4f} m, ring spacing = {2*R*np.sin(np.pi/7):.4f} m\n")
    f.write(f"satellites (az, el) = {SATS}\n")
    report(f, "base scenario", Rm, W, JAM)
    R6 = covariance(JAM_J6)
    report(f, "scenario J6 (jammer 10 deg below PRN 10)", R6, weights(R6), JAM_J6)
    f.write(f"Single element SINR = {SNR_DB - 10*np.log10(1 + len(JAM)*10**(JNR_DB/10)):.1f} dB\n")
print(open(os.path.join(HERE, "crpa8_metrics.txt")).read())
