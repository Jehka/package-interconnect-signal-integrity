#!/usr/bin/env python3
"""Per-unit-length RLGC analysis of the Q2, Q3 and Q4 extractions.

All numbers are read from the solved matrices stored in the Ansys archive
(``models/ansys/PBL1.aedtz``), so nothing is transcribed by hand. For each
coupled pair the script computes

    k_C          = |C'12| / C'11 x 100 %                 normalized capacitive coupling
    k_L          = L'12 / L'11 x 100 %                   normalized inductive coupling
    Z_c          = (Z Y)^(-1/2) Z,  Z = R + jwL, Y = G + jwC   characteristic impedance matrix
    gamma_m      = sqrt(eig(Z Y))                        modal propagation constants
    alpha_m      = Re(gamma_m)  [Np/m]  (x 8.686 -> dB/m)
    eps_eff,m    = (c0 * Im(gamma_m) / w)^2

Recomputing Z_c and alpha from the raw R, L, G, C matrices checks AEDT's own
post-processed T-line and modal data independently.

Each mode is labelled **even** or **odd** from its eigenvector (both lines in
phase vs. in anti-phase). AEDT numbers modes in solver order ("Mode 1", "Mode 2"),
and that order is not guaranteed to be the same across designs, so comparing
"Mode 1" of one design with "Mode 1" of another can compare different modes.

For the Q4 circuit the matrices are converted to lumped elements over the line
length l (default 15 mm), using the Maxwell-to-circuit conversion:

    R_line = R'11 * l
    C_m    = |C'12| * l
    C_g,i  = (C'ii - |C'12|) * l

Outputs (CSV + Markdown) go to each question's ``results/.../analysis`` folder.

Example
-------
    python scripts/rlgc_analysis.py
    python scripts/rlgc_analysis.py --length-mm 10
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import inv, sqrtm

from aedt_results import AedtArchive, default_archive_path, repo_root

C0 = 299_792_458.0
NP_TO_DB = 20 / np.log(10)  # 8.6859
ROOT = repo_root()


# --------------------------------------------------------------------------- #
# Core calculations
# --------------------------------------------------------------------------- #
def modal_analysis(R, L, G, C, f):
    """Even/odd modal attenuation, eps_eff and the characteristic impedance matrix."""
    w = 2 * np.pi * f
    Z = R + 1j * w * L
    Y = G + 1j * w * C
    ZY = Z @ Y
    Zc = inv(sqrtm(ZY)) @ Z
    out = {"Zc": Zc, "modes": {}}
    if R.shape == (1, 1):
        gamma = np.sqrt(ZY[0, 0])
        out["modes"]["single"] = _mode(gamma, w)
        return out
    vals, vecs = np.linalg.eig(ZY)
    for i in range(len(vals)):
        v = vecs[:, i] / vecs[np.argmax(np.abs(vecs[:, i])), i]
        kind = "even" if np.real(v[0] * np.conj(v[1])) > 0 else "odd"
        out["modes"][kind] = _mode(np.sqrt(vals[i]), w)
    return out


def _mode(gamma, w):
    return {
        "alpha_Np_per_m": gamma.real,
        "alpha_dB_per_m": gamma.real * NP_TO_DB,
        "eps_eff": (C0 * gamma.imag / w) ** 2,
        "v_over_c": w / gamma.imag / C0,
    }


def summarize(var):
    """Per-unit-length quantities for one solved variation."""
    R, L, G, C = (var.matrix(k) for k in "RLGC")
    f = var._archive.solution_frequency(var.design, var.vid)
    m = modal_analysis(R, L, G, C, f)
    s = {
        "design": var.design,
        "frequency_GHz": f / 1e9,
        "C11_pF_per_m": C[0, 0] * 1e12,
        "R11_ohm_per_m": R[0, 0],
        "L11_nH_per_m": L[0, 0] * 1e9,
        "G11_mS_per_m": G[0, 0] * 1e3,
        "Re_Z0_11_ohm": m["Zc"][0, 0].real,
        "_R": R, "_L": L, "_G": G, "_C": C,
    }
    if C.shape[0] >= 2:
        s.update({
            "C22_pF_per_m": C[1, 1] * 1e12,
            "C12_pF_per_m": C[0, 1] * 1e12,
            "kC_pct": abs(C[0, 1]) / C[0, 0] * 100,
            "R22_ohm_per_m": R[1, 1],
            "R12_ohm_per_m": R[0, 1],
            "L12_nH_per_m": L[0, 1] * 1e9,
            "kL_pct": L[0, 1] / L[0, 0] * 100,
        })
    for kind, md in m["modes"].items():
        for k, v in md.items():
            s[f"{kind}_{k}"] = v
    return s


def lumped(s, length_m):
    C, R = s["_C"], s["_R"]
    cm = abs(C[0, 1])
    return {
        "length_mm": length_m * 1e3,
        "R1_ohm": R[0, 0] * length_m,
        "R2_ohm": R[1, 1] * length_m,
        "Cm_pF": cm * length_m * 1e12,
        "Cg1_pF": (C[0, 0] - cm) * length_m * 1e12,
        "Cg2_pF": (C[1, 1] - cm) * length_m * 1e12,
        "C11_total_pF": C[0, 0] * length_m * 1e12,
        "R12_mutual_ohm_not_modelled": R[0, 1] * length_m,
    }


# --------------------------------------------------------------------------- #
# Output helpers
# --------------------------------------------------------------------------- #
def write_table(rows, labels, path_stem: Path, title: str, fmt: dict, notes: list[str] = ()):
    """rows: list of (quantity_label, key) ; labels: column name -> summary dict."""
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    cols = list(labels)
    with open(path_stem.with_suffix(".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity"] + cols)
        for name, key in rows:
            w.writerow([name] + [f"{float(labels[c][key]):.10g}" if key in labels[c] else "" for c in cols])
    lines = [f"# {title}", "", "| Quantity | " + " | ".join(cols) + " |", "|---" * (len(cols) + 1) + "|"]
    for name, key in rows:
        cells = [(labels[c].get("_fmt") or fmt.get(key, "{:.4g}")).format(labels[c][key])
                 if key in labels[c] else "" for c in cols]
        md_name = name.replace("|", "\\|")
        lines.append(f"| {md_name} | " + " | ".join(cells) + " |")
    if notes:
        lines += [""] + [f"- {n}" for n in notes]
    path_stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[2:]), end="\n\n")


def pct_change(a, b):
    return (b - a) / a * 100


ROWS_COUPLED = [
    ("C11 (pF/m)", "C11_pF_per_m"),
    ("C12 (pF/m)", "C12_pF_per_m"),
    ("|C12|/C11 (%)", "kC_pct"),
    ("L11 (nH/m)", "L11_nH_per_m"),
    ("L12/L11 (%)", "kL_pct"),
    ("R11 (Ω/m)", "R11_ohm_per_m"),
    ("G11 (mS/m)", "G11_mS_per_m"),
    ("Re(Z0,11) (Ω)", "Re_Z0_11_ohm"),
    ("Even-mode attenuation (dB/m)", "even_alpha_dB_per_m"),
    ("Odd-mode attenuation (dB/m)", "odd_alpha_dB_per_m"),
    ("Even-mode attenuation (Np/m)", "even_alpha_Np_per_m"),
    ("Odd-mode attenuation (Np/m)", "odd_alpha_Np_per_m"),
    ("Even-mode ε_eff", "even_eps_eff"),
    ("Odd-mode ε_eff", "odd_eps_eff"),
]
FMT = {
    "C11_pF_per_m": "{:.3f}", "C12_pF_per_m": "{:.4f}", "kC_pct": "{:.2f}",
    "L11_nH_per_m": "{:.2f}", "kL_pct": "{:.2f}", "R11_ohm_per_m": "{:.3f}",
    "G11_mS_per_m": "{:.3f}", "Re_Z0_11_ohm": "{:.3f}",
    "even_alpha_dB_per_m": "{:.3f}", "odd_alpha_dB_per_m": "{:.3f}",
    "even_alpha_Np_per_m": "{:.5f}", "odd_alpha_Np_per_m": "{:.5f}",
    "even_eps_eff": "{:.3f}", "odd_eps_eff": "{:.3f}",
}


def pick(arch, design, **match):
    for v in arch.variations(design):
        if all(v.variables.get(k) == val for k, val in match.items()):
            return v
    raise LookupError(f"{design}: no variation with {match}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", type=Path, default=default_archive_path(),
                    help=".aedtz archive or restored .aedtresults folder")
    ap.add_argument("--length-mm", type=float, default=15.0, help="Q4 line length for lumped conversion")
    ap.add_argument("--q3-mask", default="50um", help="soldermask thickness used for the Q3 comparison")
    args = ap.parse_args(argv)
    arch = AedtArchive(args.archive)

    # ---- Q2: ceramic coupled microstrip vs stripline --------------------- #
    ms = summarize(pick(arch, "Q2_Coupled_Microstrip"))
    sl = summarize(pick(arch, "Q2_Coupled_Stripline"))
    write_table(
        ROWS_COUPLED, {"Coupled microstrip": ms, "Coupled stripline": sl},
        ROOT / "results/q2_coupled_lines/analysis/q2_rlgc_comparison",
        f"Q2: ceramic coupled lines at {ms['frequency_GHz']:g} GHz", FMT,
        notes=[
            "Modes are labelled even/odd from their eigenvectors. AEDT's 'Mode 1/Mode 2' order is the "
            "solver's and differs between these two designs.",
            f"Stripline lowers |C12|/C11 from {ms['kC_pct']:.2f} % to {sl['kC_pct']:.2f} % although its "
            f"absolute |C12| is {abs(sl['C12_pF_per_m'])/abs(ms['C12_pF_per_m']):.2f}x higher.",
        ],
    )

    # ---- Q3: PWB without / with soldermask ------------------------------ #
    nm = summarize(pick(arch, "Q3_PWB_Nomask"))
    sm = summarize(pick(arch, "Q3_PWB_Soldermask", Tmask=args.q3_mask))
    delta = {k: pct_change(nm[k], sm[k]) for k in nm
             if not k.startswith("_") and isinstance(nm[k], float) and k in sm and nm[k] != 0}
    delta["C12_pF_per_m"] = pct_change(abs(nm["C12_pF_per_m"]), abs(sm["C12_pF_per_m"]))
    delta["_fmt"] = "{:+.2f}"
    write_table(
        ROWS_COUPLED, {"No mask": nm, f"{args.q3_mask} mask": sm, "Change (%)": delta},
        ROOT / "results/q3_soldermask/analysis/q3_rlgc_comparison",
        f"Q3: PWB coupled microstrip, no mask vs {args.q3_mask} soldermask", FMT,
        notes=[
            "The 'Change' column for C12 is the change in |C12|.",
            "Compare attenuation mode-for-mode: even with even, odd with odd.",
        ],
    )

    # ---- Q4: embedded coupled microstrip -> lumped RC ------------------- #
    q4 = summarize(pick(arch, "Q4_Embedded_Coupled_Microstrip"))
    lp = lumped(q4, args.length_mm * 1e-3)
    write_table(
        ROWS_COUPLED, {"Q4 embedded coupled microstrip": q4},
        ROOT / "results/q4_rc_crosstalk/analysis/q4_rlgc",
        f"Q4: Organic/Thin Film embedded coupled microstrip at {q4['frequency_GHz']:g} GHz", FMT,
    )
    rows = [
        ("Line resistance R1 (Ω)", "R1_ohm"), ("Line resistance R2 (Ω)", "R2_ohm"),
        ("Coupling capacitance Cm (pF)", "Cm_pF"),
        ("Ground capacitance Cg1 (pF)", "Cg1_pF"), ("Ground capacitance Cg2 (pF)", "Cg2_pF"),
        ("Total C11·l (pF)", "C11_total_pF"),
        ("Mutual resistance R12·l (Ω), not in circuit", "R12_mutual_ohm_not_modelled"),
    ]
    write_table(
        rows, {f"{args.length_mm:g} mm lumped": lp},
        ROOT / "results/q4_rc_crosstalk/analysis/q4_lumped_rc",
        f"Q4: lumped RC elements for a {args.length_mm:g} mm line",
        {k: "{:.5g}" for _, k in rows},
        notes=["Cg,i = (C'ii − |C'12|)·l and Cm = |C'12|·l (Maxwell → circuit capacitances)."],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
