"""Figures comparing the Salma-12 run with the CN2 control, from md_analysis.py outputs.

Usage:
    python plot_md.py --s12 salma12.npz --cn2 cn2.npz --outdir figures
Writes rmsd_compare.png and rmsf_compare.png.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

TEAL, AMBER, INK, MUTE, PAPER = "#127C74", "#C2662B", "#16262E", "#7C8B90", "#FAFAF8"


def _style(ax, xlabel, ylabel, title):
    ax.set_facecolor(PAPER)
    ax.set_title(title, fontsize=14, color=INK, pad=8)
    ax.set_xlabel(xlabel, fontsize=12, color=INK)
    ax.set_ylabel(ylabel, fontsize=12, color=INK)
    ax.tick_params(colors=INK, labelsize=11)
    ax.grid(axis="y", color=MUTE, alpha=0.15)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTE)


def plot_rmsd(s12, cn2, path):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), dpi=300)
    fig.patch.set_facecolor(PAPER)
    for a, key, title, ylabel in ((ax[0], "ligand_rmsd", "Ligand heavy-atom RMSD", "RMSD (Å)"),
                                  (ax[1], "protein_rmsd", "Protein Cα RMSD", "RMSD (Å)")):
        a.plot(s12["time_ns"], s12[key], lw=1.0, color=TEAL, label="Salma-12")
        a.plot(cn2["time_ns"], cn2[key], lw=1.0, color=AMBER, label="CN2 (control)")
        _style(a, "Time (ns)", ylabel, title)
        a.legend(frameon=False, fontsize=10)
    plt.tight_layout()
    plt.savefig(path, dpi=300, facecolor=PAPER, bbox_inches="tight")
    plt.close(fig)


def plot_rmsf(s12, cn2, path):
    fig, ax = plt.subplots(figsize=(11, 3.8), dpi=300)
    fig.patch.set_facecolor(PAPER)
    for d, color, label in ((s12, TEAL, "Salma-12"), (cn2, AMBER, "CN2 (control)")):
        ax.plot(np.arange(len(d["rmsf"])), d["rmsf"], lw=0.9, color=color, label=label)
    labels = [str(x) for x in s12["rmsf_labels"]]
    first_beta = next((i for i, l in enumerate(labels) if l.startswith("β")), None)
    if first_beta:
        ax.axvline(first_beta, color=MUTE, lw=0.8, ls=":")
        ax.text(first_beta / 2, ax.get_ylim()[1] * 0.95, "α-tubulin", ha="center", color=MUTE)
        ax.text((first_beta + len(labels)) / 2, ax.get_ylim()[1] * 0.95, "β-tubulin", ha="center", color=MUTE)
    _style(ax, "Residue index (chain order, α then β)", "Cα RMSF (Å)", "Per-residue Cα fluctuation")
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    plt.tight_layout()
    plt.savefig(path, dpi=300, facecolor=PAPER, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--s12", required=True)
    p.add_argument("--cn2", required=True)
    p.add_argument("--outdir", default="figures")
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    s12, cn2 = np.load(a.s12), np.load(a.cn2)
    plot_rmsd(s12, cn2, os.path.join(a.outdir, "rmsd_compare.png"))
    plot_rmsf(s12, cn2, os.path.join(a.outdir, "rmsf_compare.png"))
    print("wrote", os.path.join(a.outdir, "rmsd_compare.png"), "and", os.path.join(a.outdir, "rmsf_compare.png"))


if __name__ == "__main__":
    main()
