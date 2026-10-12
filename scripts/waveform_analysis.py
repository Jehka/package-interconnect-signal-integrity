#!/usr/bin/env python3
"""Timing and crosstalk metrics from an Ansys Circuit transient export.

Reads the CSV written by *Transient Voltage Plot -> Export* (columns
``Time [ns]``, ``V(Vactive) [V]``, ``V(Vquiet) [V]``) and computes

    t_pd       = t(V_out = 50 %) - t(V_in = 50 %)        propagation delay
    t_r        = t(V_out = 90 %) - t(V_out = 10 %)       receiver rise time
    V_xt,peak  = max |V_quiet|                           peak crosstalk

Crossing times are found by linear interpolation between adjacent samples, not
by snapping to the nearest sample or placing markers by hand.

Input reference
---------------
The Q4 export contains only the two receiver nodes, so the input 50 % time is
computed from the *specified ideal pulse*: t_in,50 = delay + rise/2
(0.5 ns + 0.15 ns / 2 = 0.575 ns). If an input-node waveform is exported as
well, pass ``--input-column`` and the measured crossing is used instead.

Levels
------
10/50/90 % thresholds are taken between the specified logic levels (0 V and
1.2 V by default). ``--levels measured`` uses the steady-state minimum and
maximum of the active waveform instead.

Example
-------
    python scripts/waveform_analysis.py
    python scripts/waveform_analysis.py path/to/export.csv --vhigh 1.2 --delay 0.5 --rise 0.15
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV = ROOT / "results" / "q4_rc_crosstalk" / "Q4_transient_voltage.csv"
DEFAULT_OUT = ROOT / "results" / "q4_rc_crosstalk" / "analysis"


# --------------------------------------------------------------------------- #
def load_csv(path: Path) -> dict[str, np.ndarray]:
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], np.array([[float(x) for x in r] for r in rows[1:] if r], dtype=float)
    return {h.strip(): data[:, i] for i, h in enumerate(header)}


def find_column(cols: dict[str, np.ndarray], needle: str) -> str:
    hits = [c for c in cols if needle.lower() in c.lower()]
    if len(hits) != 1:
        raise SystemExit(f"column matching {needle!r} not found uniquely in {list(cols)}")
    return hits[0]


def time_scale(column_name: str) -> float:
    """Seconds per unit of the time column, from its '[ns]'-style unit tag."""
    for unit, scale in (("[ps]", 1e-12), ("[ns]", 1e-9), ("[us]", 1e-6), ("[ms]", 1e-3), ("[s]", 1.0)):
        if unit in column_name:
            return scale
    raise SystemExit(f"cannot determine time unit of column {column_name!r}")


def crossing(t: np.ndarray, v: np.ndarray, level: float, rising: bool = True, after: float = -np.inf) -> float:
    """First time after `after` where v crosses `level`, linearly interpolated."""
    for i in range(len(t) - 1):
        if t[i + 1] < after:
            continue
        a, b = v[i], v[i + 1]
        hit = (a < level <= b) if rising else (a > level >= b)
        if hit:
            return t[i] + (level - a) * (t[i + 1] - t[i]) / (b - a)
    raise ValueError(f"no {'rising' if rising else 'falling'} crossing of {level} after t={after}")


# --------------------------------------------------------------------------- #
def analyse(args) -> dict:
    cols = load_csv(args.csv)
    tcol = find_column(cols, "time")
    acol = find_column(cols, args.active_column)
    qcol = find_column(cols, args.quiet_column)
    s = time_scale(tcol)
    t = cols[tcol] * s  # seconds
    va, vq = cols[acol], cols[qcol]

    if args.levels == "measured":
        v_lo, v_hi = float(va.min()), float(va.max())
    else:
        v_lo, v_hi = args.vlow, args.vhigh
    swing = v_hi - v_lo
    lvl = {p: v_lo + p / 100 * swing for p in (10, 50, 90)}

    if args.input_column:
        vin = cols[find_column(cols, args.input_column)]
        t_in50 = crossing(t, vin, lvl[50])
        in_ref = f"measured on column {args.input_column!r}"
    else:
        t_in50 = (args.delay + args.rise / 2) * 1e-9
        in_ref = (
            f"ideal specified pulse: delay {args.delay} ns + rise {args.rise} ns / 2 "
            "(no input-node waveform exported)"
        )

    edge_start = args.delay * 1e-9 - 1e-15  # look at the first edge of the pulse train
    t_out50 = crossing(t, va, lvl[50], after=edge_start)
    t10 = crossing(t, va, lvl[10], after=edge_start)
    t90 = crossing(t, va, lvl[90], after=edge_start)

    i_pos, i_neg = int(np.argmax(vq)), int(np.argmin(vq))
    i_abs = int(np.argmax(np.abs(vq)))

    return {
        "source_csv": str(args.csv.relative_to(ROOT) if args.csv.is_relative_to(ROOT) else args.csv),
        "samples": int(len(t)),
        "window_ns": [float(t[0] / 1e-9), float(t[-1] / 1e-9)],
        "levels_V": {"low": v_lo, "high": v_hi, "mode": args.levels},
        "input_reference": in_ref,
        "t_in_50_ns": t_in50 / 1e-9,
        "t_out_50_ns": t_out50 / 1e-9,
        "propagation_delay_ps": (t_out50 - t_in50) / 1e-12,
        "t_10_ns": t10 / 1e-9,
        "t_90_ns": t90 / 1e-9,
        "rise_time_10_90_ps": (t90 - t10) / 1e-12,
        "input_rise_10_90_ps": args.rise * 1e3 * 0.8 if not args.input_column else None,  # linear ramp
        "crosstalk_peak_abs_mV": abs(vq[i_abs]) * 1e3,
        "crosstalk_peak_abs_t_ns": t[i_abs] / 1e-9,
        "crosstalk_peak_pos_mV": vq[i_pos] * 1e3,
        "crosstalk_peak_pos_t_ns": t[i_pos] / 1e-9,
        "crosstalk_peak_neg_mV": vq[i_neg] * 1e3,
        "crosstalk_peak_neg_t_ns": t[i_neg] / 1e-9,
        "crosstalk_pct_of_swing": abs(vq[i_abs]) / swing * 100,
        "_arrays": (t, va, vq, lvl),
    }


def plot(res: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t, va, vq, lvl = res["_arrays"]
    tn = t / 1e-9
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2]})

    ax1.plot(tn, va, color="#c0392b", lw=1.6, label="V(Vactive)")
    for p in (10, 50, 90):
        ax1.axhline(lvl[p], color="0.75", lw=0.8, ls="--")
        ax1.text(tn[-1], lvl[p], f" {p}%", va="center", fontsize=8, color="0.4")
    ax1.axvline(res["t_in_50_ns"], color="0.3", lw=0.9, ls=":")
    ax1.plot([res["t_out_50_ns"]], [lvl[50]], "o", color="k", ms=4)
    ax1.annotate(f"t_pd = {res['propagation_delay_ps']:.2f} ps",
                 xy=(res["t_out_50_ns"], lvl[50]), xytext=(res["t_out_50_ns"] + 0.35, lvl[50] - 0.25),
                 arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
    ax1.annotate(f"t_r (10–90%) = {res['rise_time_10_90_ps']:.2f} ps",
                 xy=(res["t_90_ns"], lvl[90]), xytext=(res["t_90_ns"] + 0.35, lvl[90] - 0.1),
                 arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
    ax1.set_ylabel("Active receiver [V]")
    ax1.set_title("Q4 coupled embedded microstrip, 15 mm lumped RC: receiver waveforms")
    ax1.grid(alpha=0.3)

    ax2.plot(tn, vq * 1e3, color="#27ae60", lw=1.4, label="V(Vquiet)")
    ax2.axhline(0, color="0.6", lw=0.8)
    for key in ("pos", "neg"):
        ax2.plot([res[f"crosstalk_peak_{key}_t_ns"]], [res[f"crosstalk_peak_{key}_mV"]], "o", color="k", ms=4)
        ax2.annotate(f"{res[f'crosstalk_peak_{key}_mV']:+.2f} mV",
                     xy=(res[f"crosstalk_peak_{key}_t_ns"], res[f"crosstalk_peak_{key}_mV"]),
                     xytext=(8, 0), textcoords="offset points",
                     fontsize=9, va="center")
    ax2.margins(y=0.3)
    ax2.set_ylabel("Quiet line [mV]")
    ax2.set_xlabel("Time [ns]")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", type=Path, default=DEFAULT_CSV, help="Ansys transient CSV export")
    ap.add_argument("--active-column", default="Vactive")
    ap.add_argument("--quiet-column", default="Vquiet")
    ap.add_argument("--input-column", default=None, help="input-node column, if exported")
    ap.add_argument("--vlow", type=float, default=0.0, help="pulse low level [V]")
    ap.add_argument("--vhigh", type=float, default=1.2, help="pulse high level [V]")
    ap.add_argument("--delay", type=float, default=0.5, help="pulse delay [ns]")
    ap.add_argument("--rise", type=float, default=0.15, help="pulse rise time [ns], 0-100 %%")
    ap.add_argument("--levels", choices=["specified", "measured"], default="specified")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output folder")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)
    args.csv = args.csv.resolve()

    res = analyse(args)
    args.out.mkdir(parents=True, exist_ok=True)
    public = {k: v for k, v in res.items() if not k.startswith("_")}
    (args.out / "waveform_metrics.json").write_text(json.dumps(public, indent=2) + "\n")
    if not args.no_plot:
        plot(res, args.out / "waveform_metrics.png")

    print(f"Input 50 % reference : {res['t_in_50_ns']:.6f} ns  ({res['input_reference']})")
    print(f"Active output 50 %   : {res['t_out_50_ns']:.6f} ns")
    print(f"Propagation delay    : {res['propagation_delay_ps']:.2f} ps")
    print(f"10 % / 90 % crossing : {res['t_10_ns']:.6f} / {res['t_90_ns']:.6f} ns")
    print(f"10–90 % rise time    : {res['rise_time_10_90_ps']:.2f} ps")
    print(f"Peak crosstalk       : {res['crosstalk_peak_abs_mV']:.4f} mV at {res['crosstalk_peak_abs_t_ns']:.3f} ns"
          f"  ({res['crosstalk_pct_of_swing']:.4f} % of {res['levels_V']['high'] - res['levels_V']['low']:g} V)")
    print(f"  positive / negative: {res['crosstalk_peak_pos_mV']:+.4f} / {res['crosstalk_peak_neg_mV']:+.4f} mV")
    print(f"Wrote {args.out / 'waveform_metrics.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
