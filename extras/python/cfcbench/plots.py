"""Figures from the JSON lines of cfc_solve (no solver runs here).

  PYTHONPATH=python python3 -m cfcbench.plots results/*.jsonl --out figures --target 1e-4

Time to target: the fastest run of the tolerance ladder whose common relative residual
max(balance, bilateral, natural map) / velocity scale reaches the target; solver-native statuses are ignored.
"""
import argparse
import json
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Fixed solver identity (color follows the solver in every figure; markers are the second channel).
SOLVERS = {
    "cuopt": ("cuOpt (GPU)", "#2a78d6", "o"),
    "clarabel": ("Clarabel", "#eb6834", "s"),
    "conicxx": ("ConicXX", "#1baf7a", "D"),
    "qoco": ("QOCO", "#eda100", "^"),
    "scs": ("SCS", "#e87ba4", "v"),
    "mosek": ("MOSEK", "#008300", "P"),
}
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
CPU = ("clarabel", "conicxx", "qoco", "scs", "mosek")


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.size": 10, "text.color": INK, "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "lines.linewidth": 2, "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
    })


def contact_ticks(ax):
    """Plain-number ticks on a logarithmic contacts axis."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    lo, hi = ax.get_xlim()
    ticks = [t for t in (100, 200, 500, 1e3, 2e3, 3e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6) if lo <= t <= hi]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))


def footnote(fig, text):
    """Muted footnote wrapped to the figure width; reserves the space below the axes."""
    lines = []
    for part in text.split("\n"):
        lines += textwrap.wrap(part, width=int(fig.get_figwidth() * 17)) or [""]
    fig.text(0.01, 0.01, "\n".join(lines), fontsize=8, color=MUTED, va="bottom")
    fig.tight_layout(rect=(0, 0.035 * len(lines) * 4.2 / fig.get_figheight(), 1, 1))


def load(paths):
    runs = []
    for p in paths:
        for line in Path(p).read_text().splitlines():
            if line.strip():
                runs.append(json.loads(line))
    return runs


def rel(d):
    r = d.get("residuals")
    if not r or r.get("rel_max") is None:
        return np.nan
    return r["rel_max"]


def time_to_target(runs, target):
    """{(family, problem, form, solver): (nc, time or nan)} over the tolerance ladder."""
    best = {}
    for d in runs:
        k = (d.get("family"), d["problem"], d["form"], d["solver"])
        nc = d.get("nc", np.nan)
        t = d.get("t_total")
        reached = rel(d) <= target and t is not None
        prev = best.get(k, (nc, np.nan))
        best[k] = (nc, min(prev[1], t) if reached and not np.isnan(prev[1]) else (t if reached else prev[1]))
    return best


def label_ends(ax, items):
    """Direct labels at the right end of each line, nudged apart only if they collide."""
    items = sorted(items, key=lambda x: x[1])
    placed = []
    for x, y, text in items:
        ly = y
        for py in placed:
            if abs(np.log10(ly) - np.log10(py)) < 0.12:
                ly = py * 10 ** 0.12
        placed.append(ly)
        ax.annotate(text, (x, y), xytext=(x * 1.12, ly), color=INK2, fontsize=9, va="center",
                    arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.6) if ly != y else None)


def work_precision(runs, out, platform_note):
    fam = defaultdict(list)
    for d in runs:
        if d["form"] == "V" and d.get("t_total") is not None:
            fam[d.get("family")].append(d)
    families = sorted(fam)
    fig, axes = plt.subplots(1, len(families), figsize=(5.2 * len(families), 4.2), squeeze=False)
    for ax, f in zip(axes[0], families):
        # The largest instance of the family (most contacts).
        nc_max = max(d["nc"] for d in fam[f])
        prob = sorted({d["problem"] for d in fam[f] if d["nc"] == nc_max})[0]
        for s, (name, color, marker) in SOLVERS.items():
            pts = sorted((d["settings"]["tol_target"], d["t_total"], rel(d)) for d in fam[f] if d["problem"] == prob and d["solver"] == s)
            pts = [(t, r) for _, t, r in sorted(pts, reverse=True) if np.isfinite(r)]
            if not pts:
                continue
            x, y = zip(*pts)
            ax.plot(x, y, color=color, marker=marker, markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.5, label=name, zorder=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time to solution (setup + solve) [s]")
        ax.set_ylabel("achieved relative residual")
        ax.set_title(f"{f}: {prob} ({nc_max:,} contacts)", loc="left", fontsize=11, color=INK)
        ax.margins(x=0.25)
    axes[0][-1].legend(loc="upper right", fontsize=9)
    footnote(fig, "Each point: one native solver tolerance of the ladder 1e-4 ... 1e-10. Residual: max(balance, bilateral, natural map) / largest free contact velocity.\n" + platform_note)
    save(fig, out / "work_precision")


def time_vs_contacts(runs, out, target, platform_note):
    ttt = time_to_target(runs, target)
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    misses = []
    for s, (name, color, marker) in SOLVERS.items():
        pts = [(f, nc, t) for (f, p, form, solver), (nc, t) in ttt.items() if solver == s and form == "V"]
        if not pts:
            continue
        ok = sorted((nc, t) for f, nc, t in pts if np.isfinite(t))
        misses += [(nc, color, marker) for f, nc, t in pts if not np.isfinite(t)]
        if ok:
            x, y = zip(*ok)
            ax.plot(x, y, linestyle="none", color=color, marker=marker, markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.2, label=name, zorder=3)
        for f in sorted({f for f, nc, t in pts}):  # thin connector within one scene family
            fam = sorted((nc, t) for g, nc, t in pts if g == f and np.isfinite(t))
            if len(fam) >= 2 and fam[-1][0] > 1.5 * fam[0][0]:
                ax.plot(*zip(*fam), color=color, lw=1, alpha=0.5, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    if misses:
        lo, hi = ax.get_ylim()
        top = hi * 2
        for nc, color, marker in misses:
            ax.plot([nc], [top], linestyle="none", marker=marker, markerfacecolor="none", color=color, markersize=7, zorder=3)
        ax.axhline(top / 1.45, color=AXIS, lw=0.8)
        ax.text(ax.get_xlim()[0], top * 1.25, " not reached", color=MUTED, fontsize=8, va="bottom")
        ax.set_ylim(lo, top * 2.2)
    ax.set_xlabel("contacts")
    contact_ticks(ax)
    ax.set_ylabel(f"time to relative residual {target:g} [s]")
    ax.set_title(f"Time to a common residual of {target:g}", loc="left", fontsize=11)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=9)
    ax.margins(x=0.15)
    footnote(fig, "Lines connect the steps of one scene. " + platform_note)
    save(fig, out / "time_vs_contacts")


def performance_profile(runs, out, target, platform_note):
    ttt = time_to_target(runs, target)
    probs = sorted({(f, p) for (f, p, form, s) in ttt if form == "V"})
    solvers = [s for s in SOLVERS if any(k[3] == s for k in ttt)]
    T = np.full((len(probs), len(solvers)), np.inf)
    for i, (f, p) in enumerate(probs):
        for j, s in enumerate(solvers):
            nc, t = ttt.get((f, p, "V", s), (0, np.nan))
            if np.isfinite(t):
                T[i, j] = t
    best = T.min(axis=1, keepdims=True)
    keep = np.isfinite(best[:, 0])
    ratio = T[keep] / best[keep]
    taus = np.logspace(0, 2.5, 300)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for j, s in enumerate(solvers):
        name, color, marker = SOLVERS[s]
        rho = [(ratio[:, j] <= tau).mean() for tau in taus]
        ax.step(taus, rho, where="post", color=color, label=name)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("performance ratio  t / t_best")
    ax.set_ylabel("fraction of problems")
    ax.set_title(f"Performance profile, time to residual {target:g} ({keep.sum()} problems)", loc="left", fontsize=11)
    ax.legend(loc="lower right", fontsize=9)
    footnote(fig, platform_note)
    save(fig, out / "performance_profile")


def speedup(runs, out, target, platform_note):
    ttt = time_to_target(runs, target)
    rows = []
    for (f, p, form, s), (nc, t) in ttt.items():
        if s != "cuopt" or form != "V" or not np.isfinite(t):
            continue
        cpu = [ttt[(f, p, "V", c)][1] for c in CPU if (f, p, "V", c) in ttt and np.isfinite(ttt[(f, p, "V", c)][1])]
        if cpu:
            rows.append((nc, min(cpu) / t))
    if not rows:
        print("speedup: no cuOpt results reaching the target; figure skipped")
        return
    x, y = zip(*sorted(rows))
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axhline(1, color=AXIS, lw=1)
    ax.plot(x, y, linestyle="none", marker="o", color=SOLVERS["cuopt"][1], markersize=7, markeredgecolor=SURFACE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("contacts")
    contact_ticks(ax)
    ax.set_ylabel("speedup of cuOpt over the fastest CPU solver")
    ax.set_title(f"cuOpt speedup, time to residual {target:g}", loc="left", fontsize=11)
    footnote(fig, platform_note)
    save(fig, out / "speedup")


def summary(runs, out, target):
    """Markdown table per scene and solver: problems reaching the target, median time to target, median
    iterations and median achieved residual at each native tolerance."""
    ttt = time_to_target(runs, target)
    fams = sorted({k[0] for k in ttt})
    lines = [f"| scene | solver | reached {target:g} | median time to {target:g} [s] | median iterations (tol 1e-8) | median rel. residual at tol 1e-6 / 1e-8 / 1e-10 |",
             "|---|---|---|---|---|---|"]
    for f in fams:
        for s, (name, _, _) in SOLVERS.items():
            keys = [k for k in ttt if k[0] == f and k[3] == s and k[2] == "V"]
            if not keys:
                continue
            times = [ttt[k][1] for k in keys]
            ok = [t for t in times if np.isfinite(t)]
            sel = [d for d in runs if d.get("family") == f and d["solver"] == s and d["form"] == "V"]
            it = [d["iterations"] for d in sel if d["settings"]["tol_target"] == 1e-8 and d.get("iterations", -1) >= 0]
            res = []
            for tol in (1e-6, 1e-8, 1e-10):
                vals = [rel(d) for d in sel if d["settings"]["tol_target"] == tol and np.isfinite(rel(d))]
                res.append(f"{np.median(vals):.1e}" if vals else "–")
            lines.append(f"| {f} | {name} | {len(ok)}/{len(keys)} | {np.median(ok):.3g} | {np.median(it):.0f} | {' / '.join(res)} |" if ok else
                         f"| {f} | {name} | 0/{len(keys)} | – | {np.median(it):.0f} | {' / '.join(res)} |" if it else
                         f"| {f} | {name} | 0/{len(keys)} | – | – | {' / '.join(res)} |")
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def save(fig, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    fig.savefig(stem.with_suffix(".svg"))
    plt.close(fig)
    print("wrote", stem.with_suffix(".png"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", nargs="+")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--target", type=float, default=1e-4)
    ap.add_argument("--note", default="", help="e.g. 'Preliminary, CPU only, 1 thread, <cpu>'")
    args = ap.parse_args(argv)
    style()
    runs = [d for d in load(args.results) if "t_total" in d or not d.get("ok", True)]
    hosts = sorted({(d.get("platform") or {}).get("cpu", "") for d in runs} - {""})
    note = args.note or ("Platform: " + "; ".join(hosts))
    out = Path(args.out)
    work_precision(runs, out, note)
    time_vs_contacts(runs, out, args.target, note)
    performance_profile(runs, out, args.target, note)
    speedup(runs, out, args.target, note)
    summary(runs, out, args.target)


if __name__ == "__main__":
    main()
