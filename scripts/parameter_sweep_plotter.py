#!/usr/bin/env python3
"""Plot the Q1 geometry/dielectric sweeps and the Q3 soldermask sweep from solved data.

Values come from the extracted capacitance matrices stored in the Ansys archive,
not from readings of screenshot plots. For each sweep point the script also
records how well that point converged (delta energy of the last adaptive pass),
so unconverged points are marked on the plots.

Q1 (Organic/Thin Film, nominal W = 20 um, H = 12 um, er = 3.2, Wgnd = 2000 um)
    C11 vs W, H, er and ground width for stripline and embedded microstrip.
    A sweep is every variation that differs from the nominal point in exactly
    one variable, plus the nominal point itself.

    The microstrip design (Q1_Microstrip1) was saved without solution data, so
    its sweeps are not in the archive. Export them from AEDT (right-click the
    report -> Export -> .csv) into one folder as
    ``microstrip_W_sweep.csv``, ``microstrip_H_sweep.csv``, ``microstrip_er_sweep.csv``
    and ``microstrip_Wgnd_sweep.csv``, then pass ``--microstrip-csv-dir``.

Ground-plane width
    Plotted as the % deviation of C11 from the widest-ground value, with the
    +/-1 % adaptive-convergence target shaded: differences inside that band
    are not resolvable with this mesh accuracy.

Q3 (PWB, soldermask er = 3.5)
    C11 and |C12| vs soldermask thickness. 5-30 um uses the corrected
    side-mask geometry (Q3_PWB_Soldermask_Fixed), 35-60 um the covered-trace
    geometry (Q3_PWB_Soldermask), and 0 um the separate no-mask model. They are
    drawn as separate series because they are different constructions.

Outputs: CSV data in ``results/q*/analysis/`` and PNG figures next to them.

Example
-------
    python scripts/parameter_sweep_plotter.py
    python scripts/parameter_sweep_plotter.py --microstrip-csv-dir exports/
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

from aedt_results import AedtArchive, default_archive_path, repo_root

ROOT = repo_root()
Q1_OUT = ROOT / "results" / "q1_parametric_studies" / "analysis"
Q3_OUT = ROOT / "results" / "q3_soldermask" / "analysis"

NOMINAL = {"W": 20e-6, "H": 12e-6, "er": 3.2, "Wgnd": 2000e-6}
SWEEP_VARS = ["W", "H", "er", "Wgnd"]
AXIS = {
    "W": ("Trace width W [µm]", 1e6),
    "H": ("Dielectric height H [µm]", 1e6),
    "er": ("Relative permittivity εr", 1.0),
    "Wgnd": ("Ground-plane width [µm]", 1e6),
}
Q1_DESIGNS = {"Stripline": "Q1_Stripline", "Embedded microstrip": "Q1_Embedded_Microstrip"}
COLORS = {"Microstrip": "#2c7fb8", "Stripline": "#d95f0e", "Embedded microstrip": "#1a9850"}


def _close(a, b):
    return abs(a - b) <= 1e-6 * max(abs(a), abs(b), 1e-30)


# --------------------------------------------------------------------------- #
# Data collection
# --------------------------------------------------------------------------- #
def q1_points(arch: AedtArchive) -> list[dict]:
    rows = []
    for topo, design in Q1_DESIGNS.items():
        for v in arch.variations(design):
            vals = {k: v.value(k) for k in SWEEP_VARS}
            off = [k for k in SWEEP_VARS if not _close(vals[k], NOMINAL[k])]
            if len(off) > 1:
                continue
            cg = v.convergence["CG"]
            last = cg.passes[-1]
            sweeps = off if off else SWEEP_VARS  # nominal belongs to every sweep
            for sv in sweeps:
                rows.append({
                    "topology": topo,
                    "sweep": sv,
                    "W_um": vals["W"] * 1e6,
                    "H_um": vals["H"] * 1e6,
                    "er": vals["er"],
                    "Wgnd_um": vals["Wgnd"] * 1e6,
                    "C11_pF_per_m": v.matrix("C")[0, 0] * 1e12,
                    "cg_passes": len(cg.passes),
                    "cg_triangles": last.triangles,
                    "cg_delta_energy_pct": last.delta_energy_pct,
                    "cg_target_pct": cg.target_pct,
                    "converged": cg.converged,
                    "source": f"{design} (archive, {v.link})",
                })
    return rows


def read_ansys_report_csv(path: Path) -> tuple[np.ndarray, np.ndarray, str, str]:
    """First column = swept variable, second = C(SIG,SIG). Units read from headers."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f) if r]
    hx, hy = rows[0][0], rows[0][1]
    data = np.array([[float(r[0]), float(r[1])] for r in rows[1:]])
    return data[:, 0], data[:, 1], hx, hy


def microstrip_points(folder: Path) -> list[dict]:
    rows = []
    for sv in SWEEP_VARS:
        p = folder / f"microstrip_{sv}_sweep.csv"
        if not p.exists():
            print(f"  (no {p.name}; microstrip {sv} sweep skipped)")
            continue
        x, y, hx, hy = read_ansys_report_csv(p)
        xs = 1e3 if "[mm]" in hx else 1.0  # -> um for lengths
        ys = 1e3 if re.search(r"\[nF", hy) else 1.0  # -> pF/m
        for xi, yi in zip(x, y):
            vals = {k: NOMINAL[k] * (1e6 if k != "er" else 1) for k in SWEEP_VARS}
            vals[sv] = xi * (xs if sv != "er" else 1)
            rows.append({
                "topology": "Microstrip", "sweep": sv,
                "W_um": vals["W"], "H_um": vals["H"], "er": vals["er"], "Wgnd_um": vals["Wgnd"],
                "C11_pF_per_m": yi * ys, "cg_passes": "", "cg_triangles": "",
                "cg_delta_energy_pct": "", "cg_target_pct": "", "converged": "",
                "source": f"{p.name} (AEDT report export)",
            })
    return rows


def q3_points(arch: AedtArchive) -> list[dict]:
    rows = []
    groups = [
        ("Q3_PWB_Nomask", "no mask (0 µm)"),
        ("Q3_PWB_Soldermask_Fixed", "side-mask geometry (5–30 µm)"),
        ("Q3_PWB_Soldermask", "covered-trace geometry (35–60 µm)"),
    ]
    for design, label in groups:
        for v in arch.variations(design):
            C = v.matrix("C")
            cg = v.convergence["CG"]
            t = 0.0 if design == "Q3_PWB_Nomask" else v.value("Tmask") * 1e6
            rows.append({
                "geometry": label, "Tmask_um": t,
                "C11_pF_per_m": C[0, 0] * 1e12, "C12_abs_pF_per_m": abs(C[0, 1]) * 1e12,
                "kC_pct": abs(C[0, 1]) / C[0, 0] * 100,
                "cg_delta_energy_pct": cg.passes[-1].delta_energy_pct,
                "converged": cg.converged, "source": f"{design} (archive, {v.link})",
            })
    return sorted(rows, key=lambda r: r["Tmask_um"])


def write_csv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in r.items()})


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
    return plt


def plot_q1(rows, out: Path):
    plt = _plt()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    topologies = [t for t in COLORS if any(r["topology"] == t for r in rows)]
    for ax, sv in zip(axes.flat, SWEEP_VARS):
        label, scale = AXIS[sv]
        key = {"W": "W_um", "H": "H_um", "er": "er", "Wgnd": "Wgnd_um"}[sv]
        for topo in topologies:
            pts = sorted((r for r in rows if r["topology"] == topo and r["sweep"] == sv), key=lambda r: r[key])
            if not pts:
                continue
            x = np.array([r[key] for r in pts])
            y = np.array([r["C11_pF_per_m"] for r in pts])
            if sv == "Wgnd":
                y = (y / y[np.argmax(x)] - 1) * 100
            ax.plot(x, y, "-o", ms=4, color=COLORS[topo], label=topo)
            bad = [i for i, r in enumerate(pts) if r["converged"] is False]
            if bad:
                ax.plot(x[bad], y[bad], "x", ms=10, mew=2, color="k",
                        label="not converged (ΔE > 1 %)" if topo == "Stripline" else None)
        ax.set_xlabel(label)
        if sv == "Wgnd":
            ax.axhspan(-1, 1, color="0.85", zorder=0, label="±1 % convergence target")
            ax.set_ylabel("C11 change vs widest ground [%]")
            ax.set_title("C11 vs ground-plane width")
        else:
            ax.set_ylabel("C11 [pF/m]")
            ax.set_title(f"C11 vs {label.split(' [')[0].split(' (')[0]}")
        ax.legend(fontsize=8)
    fig.suptitle("Q1: self-capacitance sensitivity (Organic/Thin Film, 1 GHz, Q2D)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_q3(rows, out: Path):
    plt = _plt()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.3))
    styles = {"no mask (0 µm)": ("s", "#555555"),
              "side-mask geometry (5–30 µm)": ("o", "#2c7fb8"),
              "covered-trace geometry (35–60 µm)": ("D", "#d95f0e")}
    for geom, (mk, col) in styles.items():
        pts = [r for r in rows if r["geometry"] == geom]
        x = [r["Tmask_um"] for r in pts]
        a1.plot(x, [r["C11_pF_per_m"] for r in pts], "-" + mk, color=col, ms=5, label=geom)
        a2.plot(x, [r["C12_abs_pF_per_m"] for r in pts], "-" + mk, color=col, ms=5, label=geom)
    for ax, yl, t in ((a1, "C11 [pF/m]", "Self-capacitance"), (a2, "|C12| [pF/m]", "Mutual capacitance")):
        ax.axvline(35, color="0.6", ls=":", lw=1)
        ax.text(35.5, ax.get_ylim()[0], " T = 35 µm\n (trace top)", fontsize=8, color="0.4", va="bottom")
        ax.set_xlabel("Soldermask thickness [µm]")
        ax.set_ylabel(yl)
        ax.set_title(t)
    a1.legend(fontsize=8, loc="upper left")
    fig.suptitle("Q3: PWB coupled microstrip vs soldermask thickness (εr,mask = 3.5, 1 GHz)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", type=Path, default=default_archive_path())
    ap.add_argument("--microstrip-csv-dir", type=Path, default=None,
                    help="folder of AEDT-exported microstrip sweep CSVs")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)
    arch = AedtArchive(args.archive)

    q1 = q1_points(arch)
    if args.microstrip_csv_dir:
        q1 = microstrip_points(args.microstrip_csv_dir) + q1
    else:
        print("Microstrip sweeps are not stored in the archive (Q1_Microstrip1 has no saved solutions); "
              "pass --microstrip-csv-dir to include AEDT exports.")
    write_csv(q1, Q1_OUT / "q1_sweep_data.csv")
    q3 = q3_points(arch)
    write_csv(q3, Q3_OUT / "q3_tmask_sweep.csv")

    if not args.no_plot:
        plot_q1(q1, Q1_OUT / "q1_capacitance_sweeps.png")
        plot_q3(q3, Q3_OUT / "q3_capacitance_vs_tmask.png")

    # Console summary
    for topo in dict.fromkeys(r["topology"] for r in q1):
        print(f"\n{topo}")
        for sv in SWEEP_VARS:
            pts = sorted((r for r in q1 if r["topology"] == topo and r["sweep"] == sv),
                         key=lambda r: r[{"W": "W_um", "H": "H_um", "er": "er", "Wgnd": "Wgnd_um"}[sv]])
            if not pts:
                continue
            c = [r["C11_pF_per_m"] for r in pts]
            nc = sum(r["converged"] is False for r in pts)
            print(f"  {sv:5s} {len(pts):2d} pts  C11 {min(c):7.2f} – {max(c):7.2f} pF/m"
                  + (f"  ({nc} not converged)" if nc else ""))
    print("\nQ3 soldermask")
    for r in q3:
        print(f"  Tmask {r['Tmask_um']:4.0f} µm  C11 {r['C11_pF_per_m']:7.3f}  |C12| {r['C12_abs_pF_per_m']:6.3f} pF/m"
              f"  ({r['geometry']})")
    print(f"\nWrote {Q1_OUT.relative_to(ROOT)}/ and {Q3_OUT.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
