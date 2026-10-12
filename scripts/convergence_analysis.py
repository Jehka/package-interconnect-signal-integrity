#!/usr/bin/env python3
"""Adaptive-mesh convergence audit of every solved Q2D variation.

AEDT refines the 2D mesh pass by pass until the change in stored energy between
passes (*delta energy*, %) drops below the target (1 % here) or the pass limit
is reached. This script reads, for every solved variation in the archive,

* the full convergence table (pass, triangles, energy error %, delta energy %)
  for the CG (capacitance) and RL (resistance/inductance) solves, and
* the C11 and R11 values after every pass, from the stored per-pass matrices,

and flags

``NOT_CONVERGED``
    The last pass is still above the delta-energy target. The solve stopped on
    a limit (pass count or mesh size), not on accuracy.
``EARLY_STOP``
    Converged in fewer passes than the design's typical run, and not because a
    minimum-pass setting held it there (MinPasses was below the pass count), so
    the target may have been met by a lucky pair of coarse meshes. Worth a
    re-solve with a higher minimum. A solve that stopped exactly at its
    MinPasses is not flagged.

It also reports the last-pass change of C11 and R11 directly, because a 1 %
energy criterion does not bound every matrix entry to 1 %: the R values in this
project still move a few percent on the final pass.

Outputs, per question folder ``results/q*/analysis/``:
    convergence_summary.csv   one row per variation
    convergence_passes.csv    every adaptive pass of every variation
    q1_convergence.png        delta energy vs mesh size (Q1)

Example
-------
    python scripts/convergence_analysis.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from aedt_results import AedtArchive, default_archive_path, repo_root

ROOT = repo_root()
QUESTION_DIR = {
    "Q1": "q1_parametric_studies",
    "Q2": "q2_coupled_lines",
    "Q3": "q3_soldermask",
    "Q4": "q4_rc_crosstalk",
}
LABEL_VARS = ("W", "H", "er", "Wgnd", "Tmask")


def label(v) -> str:
    parts = []
    for k in LABEL_VARS:
        if k in v.variables:
            x = v.value(k)
            parts.append(f"{k}={x:g}" if k == "er" else f"{k}={x * 1e6:g}um")
    return " ".join(parts)


def rel_change_last(values) -> float | None:
    if len(values) < 2 or values[-2] == 0:
        return None
    return abs(values[-1] / values[-2] - 1) * 100


def audit(arch: AedtArchive, design: str):
    vs = arch.variations(design)
    typical_passes = int(np.median([len(v.convergence["CG"].passes) for v in vs]))
    summary, passes = [], []
    for v in vs:
        c11 = [m[0, 0] for m in v.pass_matrices("C")]
        r11 = [m[0, 0] for m in v.pass_matrices("R")]
        row = {"design": design, "variation": label(v), "link": v.link}
        flags = []
        for solver in ("CG", "RL"):
            t = v.convergence.get(solver)
            if not t:
                continue
            last = t.passes[-1]
            row.update({
                f"{solver}_passes": len(t.passes),
                f"{solver}_min_passes": t.min_passes,
                f"{solver}_max_passes": t.max_passes,
                f"{solver}_target_pct": t.target_pct,
                f"{solver}_final_triangles": last.triangles,
                f"{solver}_final_delta_energy_pct": last.delta_energy_pct,
                f"{solver}_final_energy_error_pct": last.energy_error_pct,
                f"{solver}_converged": t.converged,
            })
            if not t.converged:
                flags.append(f"{solver}_NOT_CONVERGED")
            elif (solver == "CG" and len(t.passes) < typical_passes
                  and (t.min_passes or 0) < len(t.passes)):
                flags.append(f"{solver}_EARLY_STOP")
            for p in t.passes:
                passes.append({
                    "design": design, "variation": label(v), "solver": solver,
                    "pass": p.pass_no, "triangles": p.triangles,
                    "total_energy": p.total_energy,
                    "energy_error_pct": p.energy_error_pct,
                    "delta_energy_pct": p.delta_energy_pct,
                    "C11_pF_per_m" if solver == "CG" else "R11_ohm_per_m":
                        (c11[p.pass_no - 1] * 1e12 if solver == "CG" else r11[p.pass_no - 1])
                        if p.pass_no - 1 < (len(c11) if solver == "CG" else len(r11)) else None,
                })
        row["C11_final_pF_per_m"] = c11[-1] * 1e12
        row["C11_last_pass_change_pct"] = rel_change_last(c11)
        row["R11_final_ohm_per_m"] = r11[-1]
        row["R11_last_pass_change_pct"] = rel_change_last(r11)
        row["flags"] = " ".join(flags)
        summary.append(row)
    return summary, passes, typical_passes


def write_csv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else ("" if v is None else v))
                        for k, v in r.items()})


def plot_q1(arch: AedtArchive, summaries: list[dict], out: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    flags = {(r["design"], r["variation"]): r["flags"] for r in summaries}
    style = {
        "": ("0.7", 0.8, "converged"),
        "CG_EARLY_STOP": ("#e67e22", 1.4, "early stop (MinPasses below pass count)"),
        "CG_NOT_CONVERGED": ("#c0392b", 1.6, "did not reach target"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, (title, design) in zip(axes, (("Stripline", "Q1_Stripline"),
                                          ("Embedded microstrip", "Q1_Embedded_Microstrip"))):
        notes = []
        for v in sorted(arch.variations(design), key=lambda v: bool(flags[(design, label(v))])):
            t = v.convergence["CG"]
            pts = [(p.triangles, p.delta_energy_pct) for p in t.passes if p.delta_energy_pct is not None]
            x, y = zip(*pts)
            f = flags[(design, label(v))].split(" ")[0] if flags[(design, label(v))] else ""
            col, lw, lab = style.get(f, style[""])
            ax.plot(x, y, "-o" if f else "-", color=col, lw=lw, ms=3, label=lab, zorder=3 if f else 2)
            if f:
                notes.append(f"{'✗' if 'NOT' in f else '!'} {label(v).replace(' Wgnd=2000um', '')}"
                             f"  ({len(t.passes)} passes, ΔE {t.passes[-1].delta_energy_pct:.2f} %)")
        ax.axhline(t.target_pct, color="k", ls="--", lw=1, label=f"target {t.target_pct:g} %")
        ax.set_yscale("log")
        ax.set_xlabel("Mesh triangles")
        ax.set_title(f"{title}: {len(arch.variations(design))} variations (CG solve)")
        ax.grid(alpha=0.3, which="both")
        h, l = ax.get_legend_handles_labels()
        uniq = dict(zip(l, h))
        ax.legend(uniq.values(), uniq.keys(), fontsize=8, loc="upper right")
        if notes:
            ax.text(0.02, 0.02, "\n".join(notes), transform=ax.transAxes, fontsize=7.5,
                    family="monospace", va="bottom",
                    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.95))
    axes[0].set_ylabel("Delta energy between passes [%]")
    fig.suptitle("Q1 adaptive-mesh convergence: delta energy vs mesh size (pass 2 onward)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", type=Path, default=default_archive_path())
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)
    arch = AedtArchive(args.archive)

    by_q: dict[str, tuple[list, list]] = {}
    print(f"{'design':32s} {'var':>4s} {'typ':>4s}  {'max ΔC11 %':>10s} {'max ΔR11 %':>10s}  flags")
    for design in arch.designs():
        q = design[:2]
        if q not in QUESTION_DIR:
            continue  # scratch designs (e.g. 2DExtractorDesign1)
        summary, passes, typical = audit(arch, design)
        s, p = by_q.setdefault(q, ([], []))
        s += summary
        p += passes
        flagged = [r for r in summary if r["flags"]]
        print(f"{design:32s} {len(summary):4d} {typical:4d}  "
              f"{max(r['C11_last_pass_change_pct'] for r in summary):10.3f} "
              f"{max(r['R11_last_pass_change_pct'] or 0 for r in summary):10.3f}  "
              f"{len(flagged)} flagged")
        for r in flagged:
            print(f"    {r['variation']:40s} {r['flags']:22s} "
                  f"CG {r['CG_passes']}/{r['CG_max_passes']} passes, {r['CG_final_triangles']} tri, "
                  f"ΔE {r['CG_final_delta_energy_pct']:.3f} %")

    for q, (summary, passes) in by_q.items():
        out = ROOT / "results" / QUESTION_DIR[q] / "analysis"
        write_csv(summary, out / "convergence_summary.csv")
        write_csv(passes, out / "convergence_passes.csv")
    if not args.no_plot:
        plot_q1(arch, by_q["Q1"][0], ROOT / "results" / QUESTION_DIR["Q1"] / "analysis" / "q1_convergence.png")

    tri_max = max(int(r["CG_final_triangles"]) for s, _ in by_q.values() for r in s)
    print(f"\nLargest final CG mesh in the project: {tri_max} triangles")
    print("Wrote results/q*/analysis/convergence_summary.csv and convergence_passes.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
