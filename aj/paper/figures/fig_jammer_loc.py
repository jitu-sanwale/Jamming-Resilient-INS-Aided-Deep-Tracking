"""Figure: jammer localization error versus observation time (R1 Section VI-C, Fig. 8).

Reads the results of aj/jammer_localization.py
(<AJ_DATA_ROOT>/campaign/loc/locres_<id>.json) for the single ground BBN
jammer of scenario S1 (seed 1) on the level-cruise and break profiles, for INS
attitude errors of 0 (truth), the receiver's initial error, and 0.1/0.5 deg in
heading, and writes fig_jammer_loc.pdf/.png and jammer_loc_metrics.txt.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.environ.get("FIG_OUT") or os.path.join(REPO, "paper", "R1")   # figure PDFs/PNGs
os.makedirs(OUT, exist_ok=True)
LOC = os.path.join(os.environ.get("AJ_DATA_ROOT") or os.path.dirname(REPO), "campaign", "loc")
items = [("S1_baseline_cruise_s01", "Level cruise"), ("S1_break_gonset_s01", "Break turn")]
cases = [("truth", "No attitude error"), ("ins", r"INS error (0.01$^\circ$/0.01$^\circ$/0.02$^\circ$)"),
         ("yaw0.1", r"Heading error 0.1$^\circ$"), ("yaw0.5", r"Heading error 0.5$^\circ$")]
colors = ["C0", "C1", "C2", "C3"]

COLUMN = os.environ.get("FIG_LAYOUT") == "column"   # one-column version (R1): panels stacked, larger text
FS = 8
if COLUMN:
    fig, axs = plt.subplots(2, 1, figsize=(3.5, 3.9), constrained_layout=True, sharex=True, sharey=True)
else:
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6), constrained_layout=True, sharey=True)
lines = []
for ax, (iid, title) in zip(axs, items):
    d = json.load(open(os.path.join(LOC, f"locres_{iid}.json")))
    j0 = d["jammers"][0]
    for (cname, label), col in zip(cases, colors):
        c = d["cases"][cname]
        kept = [t for t in c["tracks"] if t.get("kept")]
        if not kept:
            lines.append(f"{iid} | {cname}: track rejected (bearings inconsistent)")
            continue
        t = kept[0]
        tt = [h["t"] - 15.0 for h in t["history"]]
        ee = [h["err_m"] for h in t["history"]]
        ax.semilogy(tt, ee, color=col, lw=1.2, label=label)
        az = c["doa_rms_az_deg"].get(str(t["jammer"]), c["doa_rms_az_deg"].get(t["jammer"]))
        el = c["doa_rms_el_deg"].get(str(t["jammer"]), c["doa_rms_el_deg"].get(t["jammer"]))
        lines.append(f"{iid} | {cname}: final error {t['final_err_m']:.0f} m, CEP50 {t['cep50_m']:.1f} m, "
                     f"first below 100 m after {t['t_first_below_100m'] - 15.0 if t['t_first_below_100m'] else float('nan'):.1f} s, "
                     f"bearing rms az {az:.4f} deg el {el:.3f} deg, n {t['n']}")
    ax.axhline(100, color="gray", ls=":", lw=1)
    ax.grid(alpha=0.3, which="both")
    if COLUMN:
        ax.text(0.02, 0.05, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=FS,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
        ax.set_ylabel("Position error (m)", fontsize=FS)
        ax.tick_params(labelsize=FS - 0.5)
        ax.set_xlim(0, 35)
    else:
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("Time since jammer onset (s)", fontsize=8)
        ax.tick_params(labelsize=7)
if COLUMN:
    axs[1].set_xlabel("Time since jammer onset (s)", fontsize=FS)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside upper center", ncol=2, fontsize=FS - 1, frameon=False,
               handlelength=1.6, columnspacing=1.0)
    name = "fig_jammer_loc_col"
else:
    axs[0].set_ylabel("Jammer position error (m)", fontsize=8)
    axs[0].legend(fontsize=6, loc="upper right")
    name = "fig_jammer_loc"
fig.savefig(os.path.join(OUT, name + ".pdf"))
fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "jammer_loc_metrics.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
