"""Read solved 2D Extractor (Q2D) results straight out of an Ansys project archive.

No Ansys installation is needed. The reader works on either

* the archived project ``models/ansys/PBL1.aedtz`` (read in place, nothing is unzipped), or
* a restored project's results folder, e.g. ``PBL1.aedtresults/``.

What it exposes, per 2D Extractor design:

* the final (LastAdaptive) per-unit-length matrices: capacitance ``C`` [F/m],
  conductance ``G`` [S/m], resistance ``R`` [ohm/m] and inductance ``L`` [H/m];
* the same matrices after every adaptive pass;
* the adaptive-convergence table AEDT shows in *Solution Data -> Convergence*
  (pass, triangles, energy error %, delta energy %);
* for parametric sweeps, which solved variation belongs to which set of design
  variables (W, H, er, Tmask, ...).

How sweep variations are linked
-------------------------------
The archive stores each variation's matrices under a numeric variation ID
(``DV<id>_SOL10_MCAP0_V<ver>.sd``), and stores the convergence table under the
variable values, but it does not store a table connecting the two. They are
connected here by mesh fingerprint: the solver profile of each variation ID
records the triangle count of every adaptive pass, and the convergence table of
each variable set records the same sequence. A sequence match is exact and in
practice unique.

When two variations differ only in permittivity the mesh refinement is
identical (the geometry is the same), so the fingerprints tie. Those ties are
resolved physically: capacitance increases monotonically with ``er`` and ``W``
and decreases with ``H``, so the IDs are assigned in order of C11. Every link
records how it was made (``Variation.link``), and anything that cannot be linked
is reported rather than guessed.
"""

from __future__ import annotations

import os
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

RESULTS_PREFIX = "$RESULTSDIR/"

#: AEDT matrix file tags -> short names used in this module.
MATRIX_KINDS = {"MCAP": "C", "MCOND": "G", "MRES": "R", "MIND": "L"}

#: Sign of dC11/d(variable) used to break fingerprint ties.
_MONOTONIC = {"er": +1, "W": +1, "H": -1}

_LENGTH_UNITS = {"um": 1e-6, "mm": 1e-3, "nm": 1e-9, "mil": 25.4e-6, "m": 1.0}


# --------------------------------------------------------------------------- #
# Data containers
# --------------------------------------------------------------------------- #
@dataclass
class ConvergencePass:
    pass_no: int
    triangles: int
    total_energy: float | None
    energy_error_pct: float | None
    delta_energy_pct: float | None  # None on pass 1


@dataclass
class ConvergenceTable:
    solver: str  # "CG" or "RL"
    target_pct: float | None
    min_passes: int | None
    max_passes: int | None
    passes: list[ConvergencePass]

    @property
    def converged(self) -> bool:
        """True if the last pass met the delta-energy target."""
        if not self.passes or self.target_pct is None:
            return False
        last = self.passes[-1].delta_energy_pct
        return last is not None and last <= self.target_pct


@dataclass
class Variation:
    """One solved point of a design (a single variable set)."""

    design: str
    variables: dict[str, str]  # raw AEDT strings, e.g. {"W": "20um"}
    convergence: dict[str, ConvergenceTable]  # keyed "CG" / "RL"
    vid: int | None = None
    link: str = "unlinked"  # how vid was assigned
    _archive: "AedtArchive | None" = field(default=None, repr=False)

    def value(self, name: str) -> float:
        """Variable in SI units (lengths in metres, er dimensionless)."""
        return parse_quantity(self.variables[name])

    def matrix(self, kind: str = "C") -> np.ndarray:
        """Final per-unit-length matrix (C, G, R or L) in SI units."""
        if self.vid is None or self._archive is None:
            raise LookupError(f"{self.design} {self.variables}: no solved data linked")
        return self._archive.final_matrix(self.design, self.vid, kind)

    def pass_matrices(self, kind: str = "C") -> list[np.ndarray]:
        """Matrix after each adaptive pass (pass 1 first)."""
        if self.vid is None or self._archive is None:
            raise LookupError(f"{self.design} {self.variables}: no solved data linked")
        return self._archive.pass_matrices(self.design, self.vid, kind)


# --------------------------------------------------------------------------- #
# Parsing helpers
# --------------------------------------------------------------------------- #
def parse_quantity(text: str) -> float:
    """'12um' -> 1.2e-05, '3.2000000000000002' -> 3.2, '2000.0000000000002um' -> 0.002."""
    m = re.fullmatch(r"\s*([-+0-9.eE]+)\s*([a-zA-Z]*)\s*", text)
    if not m:
        raise ValueError(f"cannot parse quantity {text!r}")
    number, unit = float(m.group(1)), m.group(2)
    if not unit:
        return number
    if unit not in _LENGTH_UNITS:
        raise ValueError(f"unknown unit {unit!r} in {text!r}")
    return number * _LENGTH_UNITS[unit]


def parse_matrix_file(text: str) -> np.ndarray:
    """Parse an AEDT ``MatrixData`` .sd file into a square numpy array."""
    m = int(re.search(r"\bm=(\d+)", text).group(1))
    n = int(re.search(r"\bn=(\d+)", text).group(1))
    rows = re.findall(r"\$begin 'r\d+'(.*?)\$end 'r\d+'", text, re.S)
    data = [[float(v) for v in re.findall(r"c\(([^)]*)\)", r)] for r in rows]
    arr = np.array(data, dtype=float)
    if arr.shape != (m, n):
        raise ValueError(f"matrix shape {arr.shape} != ({m}, {n})")
    return arr


def _parse_variables(instance: str) -> dict[str, str]:
    # Instance strings look like:  $er=\\\'3.2\\\' er=\\\'3.2\\\' H=\\\'12um\\\' ...
    clean = instance.replace("\\", "").replace("'", "")
    out = {}
    for tok in clean.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            if not k.startswith("$"):  # project variables duplicate design ones
                out[k] = v
    return out


def _parse_convergence(block: str) -> dict[str, ConvergenceTable]:
    tables = {}
    for body in re.findall(r"\$begin '\d+'(.*?)\$end '\d+'", block, re.S):
        name = re.search(r"ConvSetupName='(\w+)'", body)
        if not name:
            continue
        solver = "CG" if name.group(1).upper().startswith("CG") else "RL"

        def num(key, cast=float):
            mm = re.search(rf"{key}='([^']*)'", body)
            try:
                return cast(mm.group(1)) if mm else None
            except ValueError:
                return None

        passes = []
        for entry in re.findall(r"^\s*c\((.*)\)\s*$", body, re.M):
            kv = dict(re.findall(r"(\w+)=([-+0-9.eE]+)", re.sub(r"ce='[^']*'", "", entry)))
            passes.append(
                ConvergencePass(
                    pass_no=int(kv["p"]),
                    triangles=int(kv["tri"]),
                    total_energy=float(kv["te"]) if "te" in kv else None,
                    energy_error_pct=float(kv["ee"]) if "ee" in kv else None,
                    delta_energy_pct=float(kv["de"]) if "de" in kv else None,
                )
            )
        passes.sort(key=lambda p: p.pass_no)
        tables[solver] = ConvergenceTable(
            solver=solver,
            target_pct=num("ConvTarget"),
            min_passes=num("MinPasses", int),
            max_passes=num("MaxPasses", int),
            passes=passes,
        )
    return tables


def _contains(seq: tuple, sub: tuple) -> bool:
    if not sub:
        return True
    n = len(sub)
    return any(seq[i : i + n] == sub for i in range(len(seq) - n + 1))


# --------------------------------------------------------------------------- #
# Archive reader
# --------------------------------------------------------------------------- #
class AedtArchive:
    """Read-only view of the solved results inside an .aedtz or .aedtresults folder."""

    def __init__(self, path: str | os.PathLike):
        self.path = Path(path)
        if self.path.is_file() and zipfile.is_zipfile(self.path):
            self._zip = zipfile.ZipFile(self.path)
            self._names = {
                n[len(RESULTS_PREFIX) :]: n
                for n in self._zip.namelist()
                if n.startswith(RESULTS_PREFIX) and not n.endswith("/")
            }
        elif self.path.is_dir():
            self._zip = None
            self._names = {
                p.relative_to(self.path).as_posix(): p
                for p in self.path.rglob("*")
                if p.is_file()
            }
        else:
            raise FileNotFoundError(f"{self.path} is neither an .aedtz archive nor a results folder")
        self._asol_cache: dict[str, str] = {}
        self._var_cache: dict[str, list[Variation]] = {}

    # -- low level -------------------------------------------------------- #
    def _read(self, rel: str) -> str:
        src = self._names[rel]
        if self._zip is not None:
            return self._zip.read(src).decode("latin-1")
        return Path(src).read_text(encoding="latin-1")

    def _asol(self, design: str) -> str:
        if design not in self._asol_cache:
            self._asol_cache[design] = self._read(f"{design}.asol")
        return self._asol_cache[design]

    def _solution_block(self, design: str, solution_name: str) -> str:
        asol = self._asol(design)
        for name, body in re.findall(
            r"\$begin '(SolutionBlock\d+)'(.*?)\$end '\1'", asol, re.S
        ):
            if re.search(rf"SolutionName='{solution_name}'", body):
                return body
        raise LookupError(f"{design}: no {solution_name} solution stored")

    # -- public API ------------------------------------------------------- #
    def designs(self) -> list[str]:
        """2D Extractor designs that have stored matrix solutions."""
        out = []
        for rel in self._names:
            if rel.endswith(".asol") and "/" not in rel:
                d = rel[:-5]
                if re.search(r"MCAP\(VID=", self._asol(d)):
                    out.append(d)
        return sorted(out)

    def final_matrix(self, design: str, vid: int, kind: str = "C") -> np.ndarray:
        tag = {v: k for k, v in MATRIX_KINDS.items()}[kind]
        block = self._solution_block(design, "LastAdaptive")
        m = re.search(rf"{tag}\(VID={vid}, Ver=\d+, File='([^']+)'\)", block)
        if not m:
            raise LookupError(f"{design}: no final {kind} matrix for variation {vid}")
        return parse_matrix_file(self._read(f"{design}.results/{m.group(1)}"))

    def solution_frequency(self, design: str, vid: int) -> float:
        """Extraction frequency [Hz] of the final solution."""
        block = self._solution_block(design, "LastAdaptive")
        m = re.search(rf"MCAP\(VID={vid}, Ver=\d+, File='([^']+)'\)", block)
        text = self._read(f"{design}.results/{m.group(1)}")
        return float(re.search(r"Frequency=([-+0-9.eE]+)", text).group(1))

    def pass_matrices(self, design: str, vid: int, kind: str = "C") -> list[np.ndarray]:
        tag = {v: k for k, v in MATRIX_KINDS.items()}[kind]
        block = self._solution_block(design, "AdaptivePass")
        found = []
        for sub_id, body in re.findall(r"\$begin 'SubData'\s*ID=(\d+)(.*?)\$end 'SubData'", block, re.S):
            m = re.search(rf"{tag}\(VID={vid}, Ver=\d+, File='([^']+)'\)", body)
            if m:
                found.append((int(sub_id), m.group(1)))
        return [
            parse_matrix_file(self._read(f"{self.design_dir(design)}/{f}"))
            for _, f in sorted(found)
        ]

    def design_dir(self, design: str) -> str:
        return f"{design}.results"

    def _vid_triangles(self, design: str) -> dict[int, tuple[int, ...]]:
        """Triangle count of every adaptive pass, in solve order, per variation ID."""
        out = {}
        pat = re.compile(rf"^{re.escape(design)}\.results/DV(\d+)_S\d+_V\d+\.profile$")
        for rel in self._names:
            m = pat.match(rel)
            if m:
                text = self._read(rel)
                out[int(m.group(1))] = tuple(
                    int(t) for t in re.findall(r"Triangles\\', (\d+)", text)
                )
        return out

    def _final_vids(self, design: str) -> set[int]:
        try:
            block = self._solution_block(design, "LastAdaptive")
        except LookupError:
            return set()
        return {int(v) for v in re.findall(r"MCAP\(VID=(\d+)", block)}

    def variations(self, design: str) -> list[Variation]:
        """All solved variable sets of a design, each linked to its matrices."""
        if design in self._var_cache:
            return self._var_cache[design]
        asol = self._asol(design)
        keys = re.findall(
            r"'SolutionKey\(SimSetup=\d+, Instance=\\'(.*?)\\', Solution=-3, VersionID=\d+\)'="
            r"'(SolutionBlock\d+)'",
            asol,
        )
        variations = []
        for inst, blk in keys:
            body = re.search(rf"\$begin '{blk}'(.*?)\$end '{blk}'", asol, re.S).group(1)
            variations.append(
                Variation(design, _parse_variables(inst), _parse_convergence(body), _archive=self)
            )

        final_vids = self._final_vids(design)
        tri = {v: t for v, t in self._vid_triangles(design).items() if v in final_vids}

        if len(variations) == 1 and len(final_vids) == 1:
            variations[0].vid = next(iter(final_vids))
            variations[0].link = "single variation"
        else:
            self._link_by_fingerprint(design, variations, tri)

        self._var_cache[design] = variations
        return variations

    def _link_by_fingerprint(self, design, variations, tri):
        def fp(v):
            cg = tuple(p.triangles for p in v.convergence.get("CG", ConvergenceTable("CG", None, None, None, [])).passes)
            rl = tuple(p.triangles for p in v.convergence.get("RL", ConvergenceTable("RL", None, None, None, [])).passes)
            return cg, rl

        cands = {}
        for i, v in enumerate(variations):
            cg, rl = fp(v)
            if not cg:
                continue
            cands[i] = {vid for vid, seq in tri.items() if _contains(seq, cg) and _contains(seq, rl)}

        # Iteratively peel off unique matches.
        assigned: dict[int, int] = {}
        changed = True
        while changed:
            changed = False
            for i, c in cands.items():
                if i in assigned:
                    continue
                c = c - set(assigned.values())
                if len(c) == 1:
                    assigned[i] = next(iter(c))
                    variations[i].link = "mesh fingerprint"
                    changed = True

        # Resolve remaining ties: groups of variations sharing one candidate set.
        groups: dict[frozenset, list[int]] = {}
        for i, c in cands.items():
            if i not in assigned:
                rest = frozenset(c - set(assigned.values()))
                if rest:
                    groups.setdefault(rest, []).append(i)
        for vids, idxs in groups.items():
            if len(vids) != len(idxs):
                continue
            differing = {
                k
                for k in variations[idxs[0]].variables
                if len({variations[i].variables.get(k) for i in idxs}) > 1
            }
            if len(differing) != 1 or next(iter(differing)) not in _MONOTONIC:
                continue
            var = next(iter(differing))
            sign = _MONOTONIC[var]
            by_value = sorted(idxs, key=lambda i: sign * variations[i].value(var))
            by_c11 = sorted(vids, key=lambda vid: self.final_matrix(design, vid, "C")[0, 0])
            for i, vid in zip(by_value, by_c11):
                assigned[i] = vid
                variations[i].link = f"fingerprint tie resolved by C11 monotonic in {var}"

        for i, vid in assigned.items():
            variations[i].vid = vid


def default_archive_path() -> Path:
    """models/ansys/PBL1.aedtz relative to the repository root."""
    return Path(__file__).resolve().parent.parent / "models" / "ansys" / "PBL1.aedtz"


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent
