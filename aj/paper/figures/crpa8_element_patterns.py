"""Element and array radiation patterns of the 7+1 NavIC L5 CRPA in one vertical cut.

(a) Embedded patterns of the eight patch elements (isolated patch model plus a simple
    mutual-coupling model), compared with the isolated patch.
(b) Combined pattern with uniform weights (zenith beam) and with the beam steered to
    NavIC PRN 10 (conventional beam steering, no interference).
(c) Null steering by power inversion (PI) and the per-satellite MVDR beam of channel PRN 10 adopted in
    the paper, for scenario J6: a jammer 10 deg below PRN 10 plus the two other jammers
    of the illustrative scenario.

Cut: the vertical plane through azimuth 101.2 deg (PRN 10) and 281.2 deg. Angles are
zenith angles, positive toward azimuth 101.2 deg. Array gains are referred to the
noise gain of the weights, |w^H a|^2 / (w^H w), so that uniform weights give the
classical 10 log10(8) = 9 dB of array gain over one element.

Assumptions (stated in the paper): RHCP patch amplitude pattern
G(el) = G0 ((1 + sin el)/2)^n with G0 = 5 dBic at zenith and 10 dB roll-off at the
horizon; illustrative coupling coefficient -20 dB between adjacent elements, decaying
with the square of distance, with phase -(k d + pi/2); coupling known (calibrated array);
narrowband model at the L5 carrier; JNR 50 dB per jammer.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 299792458.0
F_L5 = 1176.45e6
LAM = C / F_L5
K = 2 * np.pi / LAM
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.environ.get("FIG_OUT") or os.path.join(REPO, "paper", "R1")   # figure PDFs/PNGs
os.makedirs(OUT, exist_ok=True)

R = LAM / 2
ring = np.array([[R * np.cos(2 * np.pi * k / 7), R * np.sin(2 * np.pi * k / 7), 0.0] for k in range(7)])
POS = np.vstack([[0.0, 0.0, 0.0], ring])  # ENU metres, element 0 = centre (reference)
M = POS.shape[0]

G0_DBI = 5.0
N_ROLL = np.log10(10 ** (10 / 10)) / np.log10(2)  # 10 dB roll-off at the horizon
FLOOR_DBI = -25.0

# mutual coupling: c_mn = kappa (d_adj / d_mn)^2 exp(-j (k d_mn + pi/2)), kappa = -20 dB
D = np.linalg.norm(POS[:, None, :] - POS[None, :, :], axis=-1)
D_ADJ = 2 * R * np.sin(np.pi / 7)
KAPPA = 10 ** (-20 / 20)
CPL = np.eye(M, dtype=complex)
off = ~np.eye(M, dtype=bool)
CPL[off] = KAPPA * (D_ADJ / D[off]) ** 2 * np.exp(-1j * (K * D[off] + np.pi / 2))

AZ_CUT = 101.2
from navic_constellation import REAL, ARTIFICIAL
SATS = {**REAL, **{k.split()[0]: v for k, v in ARTIFICIAL.items()}}  # seven-satellite constellation
JAM_J6 = [(40.0, 5.0), (101.2, 14.8), (290.0, 8.0)]
JNR_DB = 50.0
SNR_DB = -30.0


def unit(az_deg, el_deg):
    az, el = np.radians(az_deg), np.radians(el_deg)
    return np.stack([np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)], axis=-1)


def elem_amp(el_deg):
    el = np.radians(el_deg)
    g = G0_DBI + 10 * N_ROLL * np.log10(np.clip((1 + np.sin(el)) / 2, 1e-12, None))
    return 10 ** (np.maximum(g, FLOOR_DBI) / 20)


def steer(az_deg, el_deg):
    """Coupled array manifold including the element pattern, shape (..., M)."""
    ideal = elem_amp(np.asarray(el_deg))[..., None] * np.exp(1j * K * (unit(az_deg, el_deg) @ POS.T))
    return ideal @ CPL.T


def cut_dirs(zen_deg):
    """Zenith angle in the cut (positive toward AZ_CUT) -> (az, el)."""
    az = np.where(zen_deg >= 0, AZ_CUT, (AZ_CUT + 180) % 360)
    return az, 90 - np.abs(zen_deg)


def gain_db(w, A):
    return 10 * np.log10(np.abs(A @ w.conj()) ** 2 / np.real(np.vdot(w, w)) + 1e-15)


zen = np.linspace(-90, 90, 18001)
AZc, ELc = cut_dirs(zen)
A_cut = steer(AZc, ELc)

# (b) uniform and conventional beam steering toward PRN 10
w_uni = np.ones(M, dtype=complex)
w_bs = steer(*SATS["PRN 10"])

# (c) J6 interference covariance, PI and per-satellite MVDR beam of channel PRN 10
Rm = np.eye(M, dtype=complex) + sum(10 ** (JNR_DB / 10) * np.outer(steer(*j), steer(*j).conj()) for j in JAM_J6)
Ri = np.linalg.inv(Rm)
e0 = np.zeros(M); e0[0] = 1.0
w_pi = Ri @ e0 / (e0 @ Ri @ e0)
a10_ = steer(*SATS["PRN 10"])
w_lc = Ri @ a10_ / (a10_.conj() @ Ri @ a10_)  # per-satellite MVDR beam of channel PRN 10

th = np.radians(zen)
fig = plt.figure(figsize=(7.16, 5.0), constrained_layout=True)
axs = [fig.add_subplot(2, 2, k + 1, projection="polar") for k in range(3)] + [fig.add_subplot(2, 2, 4)]
RMIN, RMAX = -40, 15


def setup(ax, title, rmin=RMIN, rmax=RMAX, rticks=(-30, -15, 0, 10)):
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_thetamin(-90); ax.set_thetamax(90)
    ax.set_rlim(rmin, rmax)
    ax.set_rticks(list(rticks)); ax.set_rlabel_position(80); ax.tick_params(labelsize=6)
    ax.set_xticks(np.radians([-90, -60, -30, 0, 30, 60, 90]))
    ax.set_xticklabels(["H", "60°", "30°", "Z", "30°", "60°", "H"], fontsize=6)
    ax.set_title(title, fontsize=8, pad=2)


def mark(ax, zen_deg, style, label=None, rmin=RMIN, rmax=RMAX):
    ax.plot([np.radians(zen_deg)] * 2, [rmin, rmax], style, lw=0.9, label=label)


iso = 10 * np.log10(elem_amp(ELc) ** 2)
setup(axs[0], "(a) Embedded element patterns", -15, 10, (-10, -5, 0, 5))
for m in range(M):
    g = 10 * np.log10(np.abs(A_cut[:, m]) ** 2)
    axs[0].plot(th, np.clip(g, -15, None), lw=0.9, color="0.15" if m == 0 else plt.cm.viridis(m / 8),
                label="Element 0 (center)" if m == 0 else f"Element {m}")
axs[0].plot(th, np.clip(iso, -15, None), "k--", lw=1.1, label="Isolated patch")

setup(axs[1], "(b) Combined pattern: uniform and beam steering")
axs[1].plot(th, np.clip(gain_db(w_uni, A_cut), RMIN, None), color="tab:blue", lw=1.1, label="Uniform weights")
axs[1].plot(th, np.clip(gain_db(w_bs, A_cut), RMIN, None), color="tab:orange", lw=1.1, label="Beam steered to PRN 10")
axs[1].plot(th, np.clip(iso, RMIN, None), "k--", lw=0.9, label="Single element")
mark(axs[1], 90 - 24.8, "g:", "PRN 10")

setup(axs[2], "(c) PI vs. per-satellite beam of PRN 10, J6")
axs[2].plot(th, np.clip(gain_db(w_pi, A_cut), RMIN, None), color="tab:red", lw=1.1, label="PI (null steering)")
axs[2].plot(th, np.clip(gain_db(w_lc, A_cut), RMIN, None), color="tab:purple", lw=1.1, label="Per-satellite beam (PRN 10)")
mark(axs[2], 90 - 24.8, "g:", "PRN 10")
mark(axs[2], 90 - 14.8, "r:", "Jammer")

ax = axs[3]
sel = zen >= 0
el_cut = 90 - zen[sel]
ax.plot(el_cut, gain_db(w_pi, A_cut[sel]), color="tab:red", lw=1.1, label="PI (null steering)")
ax.plot(el_cut, gain_db(w_lc, A_cut[sel]), color="tab:purple", lw=1.1, label="Per-satellite beam (PRN 10)")
ax.plot(el_cut, gain_db(w_bs, A_cut[sel]), color="tab:orange", lw=0.9, ls="--", label="Beam to PRN 10 (no jamming)")
ax.axvline(24.8, color="g", ls=":", lw=1.0, label="PRN 10")
ax.axvline(14.8, color="r", ls=":", lw=1.0, label="Jammer")
ax.set_xlim(0, 60); ax.set_ylim(-120, 15)
ax.set_xlabel("Elevation toward azimuth 101.2° (deg)", fontsize=7)
ax.set_ylabel("Gain (dB)", fontsize=7)
ax.set_title("(d) Detail of (c) near PRN 10 and the jammer", fontsize=8)
ax.grid(alpha=0.3); ax.tick_params(labelsize=6)

axs[0].legend(fontsize=5, loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=3, framealpha=0.8)
for a in axs[1:3]:
    a.legend(fontsize=5.5, loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2, framealpha=0.8)
ax.legend(fontsize=5.5, loc="lower left", framealpha=0.8)
fig.savefig(os.path.join(OUT, "fig_crpa8_elements.pdf"))
fig.savefig(os.path.join(OUT, "fig_crpa8_elements.png"), dpi=300)

# metrics
with open(os.path.join(HERE, "crpa8_element_metrics.txt"), "w") as f:
    emb = np.array([10 * np.log10(np.abs(A_cut[:, m]) ** 2) for m in range(M)])
    upper = np.abs(zen) <= 80
    f.write(f"ring radius = {R:.4f} m, adjacent spacing = {D_ADJ:.4f} m, coupling adjacent = -20 dB\n")
    f.write(f"embedded element gain spread over |zenith| <= 80 deg: max {np.max(emb[:, upper].max(0) - emb[:, upper].min(0)):.2f} dB\n")
    f.write(f"embedded zenith gains (dBic): {[round(float(x), 2) for x in emb[:, np.argmin(np.abs(zen))]]}\n")
    a10 = steer(*SATS["PRN 10"]); aj = steer(*JAM_J6[1])
    for name, w in (("uniform", w_uni), ("beam PRN 10", w_bs), ("PI", w_pi), ("per-sat beam PRN 10", w_lc)):
        z = gain_db(w, A_cut)
        f.write(f"{name}: peak {z.max():.1f} dB at zenith angle {zen[np.argmax(z)]:.1f}; "
                f"gain to PRN 10 {gain_db(w, a10[None, :])[0]:.1f} dB; gain to J6 jammer {gain_db(w, aj[None, :])[0]:.1f} dB\n")
    for name, w in (("PI", w_pi), ("per-sat beam PRN 10", w_lc)):
        s = 10 ** (SNR_DB / 10) * abs(np.vdot(w, a10)) ** 2 / np.real(w.conj() @ Rm @ w)
        f.write(f"{name}: output SINR PRN 10 = {10 * np.log10(s):.1f} dB\n")
print(open(os.path.join(HERE, "crpa8_element_metrics.txt")).read())
