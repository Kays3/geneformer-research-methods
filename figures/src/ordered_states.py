"""Regenerate the synthetic animated state geometry and its static companion."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

OUT = Path(__file__).resolve().parents[1]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "text.color": "#142a39",
        "axes.titlecolor": "#142a39",
    }
)
fig = plt.figure(figsize=(9.0, 4.6), dpi=100, facecolor="#f8faf7")
ax = fig.add_axes([0.065, 0.18, 0.44, 0.62], facecolor="#f8faf7")
text = fig.add_axes([0.56, 0.18, 0.42, 0.62])
text.axis("off")
refs = np.array([12, 48, 84])


def unit(deg):
    rad = np.radians(deg)
    return np.stack([np.cos(rad), np.sin(rad)], axis=-1)


experiments = [
    ("High → intermediate", 19, 21, 27, 0, 1),
    ("Intermediate → low", 52, 54, 74, 1, 2),
]


def draw(frame):
    part = min(frame // 36, 1)
    t = min(1, (frame % 36) / 23)
    title, start, ctrl, target, s, g = experiments[part]
    ax.clear()
    text.clear()
    text.axis("off")
    ax.set(xlim=(0, 1.20), ylim=(0, 1.19), aspect="equal")
    ax.axis("off")
    arc = unit(np.linspace(0, 90, 100))
    ax.plot(*arc.T, c="#9daeb1", lw=1)
    for i, (label, pos) in enumerate(zip(("High", "Mid", "Low"), unit(refs))):
        ax.plot([0, pos[0]], [0, pos[1]], lw=1, color="#c6d4d1")
        ax.scatter(*pos, s=68, color="#193e49", zorder=4)
        ax.annotate(label, pos, xytext=(3, 8), textcoords="offset points", fontsize=10)
    for angle, color, mark in [(ctrl, "#aa5a25", "s"), (target, "#087d76", "o")]:
        active = start + t * (angle - start)
        xy = unit(np.linspace(start, active, 35)) * 0.85
        ax.plot(*xy.T, c=color, lw=2)
        ax.scatter(*xy[-1], c=color, marker=mark, s=65, zorder=5)
    ax.scatter(
        *(unit(start) * 0.85), facecolors="white", edgecolors="#203746", s=34, zorder=6
    )
    text.text(0, 1, title, fontsize=17, weight="bold", va="top")
    text.text(0, 0.84, "Separate starting cohort", fontsize=10, color="#47626c")
    for y, label, color in [
        (0.68, "●  Target edit", "#087d76"),
        (0.56, "■  Matched control", "#aa5a25"),
    ]:
        text.text(0, y, label, color=color, fontsize=12, weight="bold")
    target_score = unit(refs) @ unit(start + t * (target - start))
    ctrl_score = unit(refs) @ unit(start + t * (ctrl - start))
    direct = target_score[g] - ctrl_score[g]
    relative = direct - (target_score[s] - ctrl_score[s])
    text.text(0, 0.37, f"Direct goal shift     {direct:+.3f}", fontsize=12)
    text.text(0, 0.23, f"Goal minus start    {relative:+.3f}", fontsize=12)
    text.text(
        0,
        0.05,
        "Nearest target state    "
        + ["high", "intermediate", "low"][np.argmax(target_score)],
        fontsize=11,
    )


fig.text(
    0.055,
    0.925,
    "How a matched edit changes state similarity",
    fontsize=20,
    weight="bold",
)
fig.text(
    0.055,
    0.06,
    "Synthetic unit-vector geometry   ·   Motion illustrates each comparison; frames do not represent biological time.",
    fontsize=9,
)
ani = FuncAnimation(fig, draw, frames=72, interval=110)
ani.save(OUT / "08-state-shifts.gif", writer=PillowWriter(fps=9), dpi=95)
draw(65)
fig.savefig(OUT / "08-state-shifts.png", dpi=145)
