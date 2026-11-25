"""Plot a completed benchmark's donor estimates, recall and assignment crossings."""

import argparse
import csv
import json
from pathlib import Path

from .cli import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path, help="New PNG, SVG or PDF")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new output")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    receipt = json.loads((args.run / "receipt.json").read_text())
    path = args.run / "donor_summary.tsv"
    if not path.exists() or receipt["outputs"].get(path.name) != sha(path):
        parser.error("No verified donor summary is available")
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    steps = [("high", "mid"), ("mid", "low")]
    fig, axs = plt.subplots(1, 3, figsize=(11, 4.3), layout="constrained")
    for i, (start, goal) in enumerate(steps):
        group = [r for r in rows if (r["start"], r["goal"]) == (start, goal)]
        for j, r in enumerate(group):
            x = i + 0.28 * ((j + 0.5) / len(group) - 0.5)
            vals = [float(r["direct"]), float(r["relative"])]
            axs[0].plot([x, x], vals, color="#a6b5ba", lw=1)
            for v, c, m in zip(vals, ["#137d78", "#ab512b"], ["o", "s"]):
                axs[0].scatter(x, v, c=c, marker=m, s=30)
        denom = sum(int(r["aligned"]) for r in group)
        for j, key in enumerate(["target_cross", "control_cross"]):
            color = ["#137d78", "#84939b"][j]
            x = i + (0.17 if j else -0.17)
            if denom:
                axs[2].bar(
                    x, sum(float(r[key]) for r in group) / denom, 0.28, color=color
                )
            else:
                axs[2].text(x, 0.02, "NA", ha="center", fontsize=8)
        axs[2].text(i, 1.02, f"aligned n={denom}", ha="center", fontsize=8)
    for ax in (axs[0], axs[2]):
        ax.set_xticks([0, 1], ["High → mid", "Mid → low"])
        ax.set_xlim(-0.55, 1.55)
    axs[0].axhline(0, color="#41596a", lw=0.8)
    axs[0].set_title("Shift by donor")
    axs[0].set_ylabel("Control-adjusted cosine difference")
    axs[0].scatter([], [], c="#137d78", label="Direct")
    axs[0].scatter([], [], c="#ab512b", marker="s", label="Relative")
    axs[0].legend(fontsize=8)
    for i, s in enumerate(["high", "mid", "low"]):
        recall = receipt["per_state_recall"][s]
        if recall is not None:
            axs[1].bar(i, recall, color="#33566c")
        else:
            axs[1].text(i, 0.02, "NA", ha="center")
    axs[1].set_xticks(range(3), ["High", "Mid", "Low"])
    axs[1].set_xlim(-0.55, 2.55)
    axs[1].set_title("Unedited state recall")
    axs[2].set_title("Assignment crossings")
    axs[2].legend(
        handles=[
            Patch(facecolor="#137d78", label="Target"),
            Patch(facecolor="#84939b", label="Control"),
        ],
        fontsize=8,
    )
    for ax in axs[1:]:
        ax.set_ylim(0, 1.12)
        ax.set_ylabel("Fraction")
    fig.suptitle(
        f"{receipt['target']}  ·  {receipt['model']}\n{receipt['scope']}  ·  {receipt['eligible_cells']} / {receipt['evaluation_cells']} evaluation cells eligible",
        fontsize=12,
    )
    fig.savefig(args.output, dpi=145)
    plt.close(fig)
