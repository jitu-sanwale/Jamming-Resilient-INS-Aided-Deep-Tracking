"""Effect of aircraft roll rate on CRPA null depth, with and without INS attitude aiding.

A single ground jammer is fixed in the local-level frame while the aircraft
rolls about its longitudinal axis. The 7+1 array is fixed in the body frame,
so the jammer direction in the body frame rotates at the roll rate.

Conventional processing: power-inversion (PI) weights from the block-averaged
covariance of the previous block, held for the next block (block length T_b).
INS-propagated processing: the jammer direction estimated from a 10-ms
covariance block (with an estimation error) is held in the navigation frame,
rotated into the body frame with the INS attitude (with an attitude error),
and PI weights are recomputed every T_u from the nominal array model. The true
array has per-element gain/phase errors; blind PI adapts to them, whereas the
model-based propagation does not, so its performance depends on calibration.

Output: output jammer-to-noise ratio (INR) after PI versus roll rate
(fig_roll_nulling.pdf/.png and crpa8_roll_metrics.txt).
Assumptions: isotropic elements, exact covariances, narrowband model.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 299792458.0
F_L5 = 1176.45e6
LAM_L5 = C / F_L5
R_RING = LAM_L5 / 2  # NavIC L5 design: ring radius = lambda_L5 / 2
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.environ.get("FIG_OUT") or os.path.join(REPO, "paper", "R1")   # figure PDFs/PNGs
os.makedirs(OUT, exist_ok=True)

ring = np.array([[R_RING * np.cos(2 * np.pi * k / 7), R_RING * np.sin(2 * np.pi * k / 7), 0.0] for k in range(7)])
POS = np.vstack([[0.0, 0.0, 0.0], ring])  # body frame (x fwd, y right, z down), element 0 = centre
M = POS.shape[0]
E0 = np.zeros(M); E0[0] = 1.0

JNR = 10 ** (60 / 10)
JAM_AZ, JAM_EL = 60.0, 10.0      # jammer direction in the local-level frame (deg)
DOA_ERR_DEG = 0.2                # nav-frame jammer DOA estimation error (INS-aided case)
ATT_ERR_DEG = 0.01               # INS attitude error (as in the UT paper's initialization)
RNG = np.random.default_rng(1)
CAL_ERR = {"uncal": (0.5, 3.0), "cal": (0.05, 0.3)}  # rms gain (dB), phase (deg) per element
ERR = {}
for k, (g_db, p_deg) in CAL_ERR.items():
    r = np.random.default_rng(7)
    e = 10 ** (g_db * r.standard_normal(M) / 20) * np.exp(1j * np.radians(p_deg) * r.standard_normal(M))
    e[0] = 1.0
    ERR[k] = e


def ned_unit(az, el):
    a, e = np.radians(az), np.radians(el)
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), -np.sin(e)])


def rot_x(phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])  # NED -> body for a pure roll phi


def steer_body(u_b, lam):
    return np.exp(1j * 2 * np.pi / lam * (POS @ u_b))


def pi_weights(R):
    Ri = np.linalg.inv(R)
    return Ri @ E0 / (E0 @ Ri @ E0)


def out_inr(w, u_list, lam, err):
    # mean output jammer power over the block divided by output noise power
    j = np.mean([abs(np.vdot(w, err * steer_body(u, lam))) ** 2 for u in u_list])
    return JNR * j / np.real(np.vdot(w, w))


def block_cov(u_list, lam, err):
    R = np.eye(M, dtype=complex)
    for u in u_list:
        a = err * steer_body(u, lam)
        R += JNR / len(u_list) * np.outer(a, a.conj())
    return R


def run(rate_dps, lam, mode, Tb=1e-3, Tu=1e-4, cal="uncal", n_start=36, n_sub=20):
    err = ERR[cal]
    u_n = ned_unit(JAM_AZ, JAM_EL)
    w_rate = np.radians(rate_dps)
    out = []
    for phi0 in np.linspace(0, 2 * np.pi, n_start, endpoint=False):
        if mode == "conv":
            t_prev = np.linspace(-Tb, 0, n_sub)
            t_next = np.linspace(0, Tb, n_sub)
            w = pi_weights(block_cov([rot_x(phi0 + w_rate * t) @ u_n for t in t_prev], lam, err))
            out.append(out_inr(w, [rot_x(phi0 + w_rate * t) @ u_n for t in t_next], lam, err))
        else:
            # nav-frame DOA with estimation error, body prediction with attitude error
            d = np.radians(DOA_ERR_DEG) * RNG.standard_normal(2)
            u_hat_n = ned_unit(JAM_AZ + np.degrees(d[0]) / max(np.cos(np.radians(JAM_EL)), 1e-3), JAM_EL + np.degrees(d[1]))
            dphi = np.radians(ATT_ERR_DEG) * RNG.standard_normal()
            u_hat_b = rot_x(phi0 + dphi) @ u_hat_n
            a = steer_body(u_hat_b, lam)
            w = pi_weights(np.eye(M) + JNR * np.outer(a, a.conj()))
            t_next = np.linspace(0, Tu, max(2, n_sub // 4))
            out.append(out_inr(w, [rot_x(phi0 + w_rate * t) @ u_n for t in t_next], lam, err))
    return 10 * np.log10(np.mean(out))


rates = np.array([0, 25, 50, 100, 150, 200, 300, 450])
cases = [("Conventional PI, $T_b$ = 10 ms", "conv", dict(Tb=10e-3)),
         ("Conventional PI, $T_b$ = 1 ms", "conv", dict(Tb=1e-3)),
         ("INS-propagated PI, uncalibrated", "ins", dict(Tu=1e-4, cal="uncal")),
         ("INS-propagated PI, calibrated", "ins", dict(Tu=1e-4, cal="cal"))]

fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.6), constrained_layout=True)
lines = []
for (label, mode, kw), mk in zip(cases, ("o", "s", "^", "D")):
    y = [run(r, LAM_L5, mode, **kw) for r in rates]
    ax.plot(rates, y, marker=mk, ms=4, label=label)
    lines.append(f"NavIC L5 | {label}: " + ", ".join(f"{r}:{v:.1f}" for r, v in zip(rates, y)))
ax.axvline(447, color="gray", ls=":", lw=1)
ax.set_xlabel("Roll rate (deg/s)", fontsize=8)
ax.set_ylabel("Output INR (dB)", fontsize=8)
ax.grid(alpha=0.3)
ax.tick_params(labelsize=7)
ax.legend(fontsize=6, loc="upper left")
fig.savefig(os.path.join(OUT, "fig_roll_nulling.pdf"))
fig.savefig(os.path.join(OUT, "fig_roll_nulling.png"), dpi=300)

with open(os.path.join(HERE, "crpa8_roll_metrics.txt"), "w") as f:
    f.write(f"ring radius = {R_RING:.4f} m; ring spacing = {2*R_RING*np.sin(np.pi/7):.4f} m "
            f"= {2*R_RING*np.sin(np.pi/7)/LAM_L5:.3f} lambda_L5\n")
    f.write("output INR (dB) vs roll rate (deg/s); input JNR 60 dB\n")
    f.write("\n".join(lines) + "\n")
print(open(os.path.join(HERE, "crpa8_roll_metrics.txt")).read())
