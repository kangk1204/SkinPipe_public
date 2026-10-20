'Uni-Dock: the Vina scoring function on a GPU, with a caveat attached.'
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from ..types import Receptor

log = logging.getLogger(__name__)

#: Uni-Dock's own search presets. `exhaustiveness` exists on the command line
#: but does not behave like Vina's -- see the module docstring.
SEARCH_MODES = ("fast", "balance", "detail")

DEFAULT_SEARCH_MODE = "balance"

#: Meeko's macrocycle pseudo-atoms. To dock a macrocycle, meeko breaks one ring
#: bond and adds a pair of these with a "glue" potential that pulls the opening
#: shut again during the search. The potential is a Vina 1.2 feature; Uni-Dock
#: does not implement it, so these atoms arrive unparameterised and the energy
#: is meaningless. This -- not ligand size -- is what the GPU path gets wrong.
GLUE_ATOM_TYPES = frozenset({"CG", "CG0", "CG1", "CG2", "CG3", "G", "G0", "G1",
                             "G2", "G3"})

#: Operational review/routing window, motivated by +110.03 and -66.69 in a
#: historical batch. Finite values outside it are retained as flagged raw
#: outputs; this is not a universal physical affinity bound.
PLAUSIBLE_AFFINITY = (-20.0, 0.0)

#: Ligands per GPU batch. Bounds VRAM (300 ligands used 5.9 GB of 12 GB on an
#: RTX 4070) and keeps the argument list well short of ARG_MAX.
BATCH_SIZE = 300


def is_plausible(affinity: float) -> bool:
    """Does this score fall inside the current operational review window?"""
    low, high = PLAUSIBLE_AFFINITY
    return low <= float(affinity) < high


def torsion_count(pdbqt: str) -> int:
    """Rotatable bonds in a prepared ligand, as the PDBQT records them.

    Kept for the log line, not as a gate: a torsion count of 9 was the first
    boundary measured here and it was a confound. See `has_glue_atoms`.
    """
    return sum(1 for line in pdbqt.splitlines() if line.startswith("BRANCH"))


def has_glue_atoms(pdbqt: str) -> bool:
    'Was this ligand prepared as an opened macrocycle?'
    for line in pdbqt.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        # The AutoDock type is the last field of a PDBQT atom line. Read it as a
        # token rather than by column: the type occupies columns 78-79 in the
        # spec but meeko writes three-character glue types like `CG0`, and a
        # fixed two-column slice truncates those to `CG`. That happens to still
        # match here, and would stop matching the day a type is added whose
        # first two characters are not themselves a glue type.
        fields = line.split()
        if fields and fields[-1] in GLUE_ATOM_TYPES:
            return True
    return False


def best_affinity(docked_pdbqt: str) -> Optional[float]:
    """The best `REMARK VINA RESULT` in a docked PDBQT, or None if there is none."""
    best: Optional[float] = None
    for line in docked_pdbqt.splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            try:
                value = float(line.split(":", 1)[1].split()[0])
            except (IndexError, ValueError):
                continue
            best = value if best is None else min(best, value)
    return best


class UnidockUnavailable(RuntimeError):
    """No `unidock` executable, or it will not run. Distinct from a bad dock."""


def executable(explicit: Optional[str] = None) -> str:
    """The `unidock` binary, or refuse with something actionable.

    It is a separate conda package with a CUDA build and is deliberately not a
    dependency of this one: installing it pulls a CUDA runtime, and the
    pipeline runs perfectly well without a GPU.
    """
    found = explicit or os.environ.get("UNIDOCK_BIN") or shutil.which("unidock")
    if not found:
        raise UnidockUnavailable(
            "no `unidock` on PATH. It ships separately:\n"
            "    mamba create -n unidock -c conda-forge unidock\n"
            "then point UNIDOCK_BIN at "
            "<env>/bin/unidock, or pass --unidock-bin.")
    return found


def dock_batch(
    receptor: Receptor,
    ligand_pdbqts: "list[Path]",
    *,
    search_mode: str = DEFAULT_SEARCH_MODE,
    n_poses: int = 9,
    seed: int = 42,
    binary: Optional[str] = None,
    timeout_s: int = 3600,
    batch_size: int = BATCH_SIZE,
) -> "dict[str, str]":
    """Dock many prepared ligands per GPU launch. Returns {stem: docked PDBQT}.

    This is the only mode worth using the GPU for. One ligand per invocation
    pays the grid build and the CUDA context every time: measured here, 300
    ligands went through one `--gpu_batch` in 90 seconds -- 200 per minute,
    against 4.9 per minute for CPU Vina across 31 worker processes on 32 cores.

    A ligand missing from the returned mapping is a ligand Uni-Dock wrote no
    output for. That is reported as absence by omission, so the caller can
    re-dock it rather than receive a zero or a guess. Nothing here judges the
    numbers -- `is_plausible` is the caller's to apply, because the caller is
    the one that can fall back.
    """
    exe = executable(binary)
    if search_mode not in SEARCH_MODES:
        raise ValueError(f"search_mode must be one of {SEARCH_MODES}, "
                         f"got {search_mode!r}")
    if receptor.pdbqt is None or not Path(receptor.pdbqt).is_file():
        raise FileNotFoundError(f"receptor PDBQT missing for {receptor.pdb_id}")
    paths = [Path(p).resolve() for p in ligand_pdbqts]
    if not paths:
        return {}

    cx, cy, cz = (float(v) for v in receptor.center)
    sx, sy, sz = (float(v) for v in receptor.size)
    out: dict[str, str] = {}
    for start in range(0, len(paths), max(1, int(batch_size))):
        chunk = paths[start:start + max(1, int(batch_size))]
        with tempfile.TemporaryDirectory(prefix="unidock_batch_") as tmp:
            out_dir = Path(tmp) / "out"
            out_dir.mkdir()
            cmd = [
                exe,
                "--receptor", str(Path(receptor.pdbqt).resolve()),
                "--gpu_batch", *[str(p) for p in chunk],
                "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}",
                "--center_z", f"{cz:.3f}",
                "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}",
                "--size_z", f"{sz:.3f}",
                "--search_mode", search_mode,
                "--num_modes", str(int(n_poses)),
                "--seed", str(int(seed)),
                "--dir", str(out_dir),
            ]
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True,
                                      timeout=timeout_s)
            except FileNotFoundError as exc:
                raise UnidockUnavailable(f"could not execute {exe}: {exc}") from exc
            except subprocess.TimeoutExpired as exc:
                raise UnidockUnavailable(
                    f"unidock did not finish {len(chunk)} ligand(s) within "
                    f"{timeout_s}s") from exc
            if proc.returncode != 0:
                raise UnidockUnavailable(
                    f"unidock exited {proc.returncode} on a batch of "
                    f"{len(chunk)}: {(proc.stderr or proc.stdout or '')[:300]}")
            # `<stem>_out.pdbqt`, and the stems are the ligand ids the caller
            # passed in. Matched by name rather than by position: a batch that
            # silently drops one would otherwise shift every later ligand onto
            # the wrong molecule.
            wanted = {p.stem for p in chunk}
            for produced in sorted(out_dir.glob("*.pdbqt")):
                stem = produced.stem
                if stem.endswith("_out"):
                    stem = stem[: -len("_out")]
                if stem in wanted:
                    out[stem] = produced.read_text(encoding="utf-8")
            missing = wanted - set(out)
            if missing:
                log.warning("unidock returned no pose for %d of %d ligand(s) in "
                            "this batch (%s%s); they are reported as undocked, "
                            "not as zero", len(missing), len(chunk),
                            ", ".join(sorted(missing)[:4]),
                            "..." if len(missing) > 4 else "")
    return out


def dock_one_text(
    receptor: Receptor,
    ligand_pdbqt: Path,
    *,
    search_mode: str = DEFAULT_SEARCH_MODE,
    n_poses: int = 9,
    seed: int = 42,
    binary: Optional[str] = None,
    timeout_s: int = 900,
) -> str:
    """Dock one prepared ligand and return the docked PDBQT text.

    The return type is the same as `engine.dock_one`'s intermediate -- docked
    PDBQT text with Vina REMARK lines -- so `_poses_from_docked_text` parses
    either engine's output with no branch. That is the whole reason this
    function exists at this seam rather than deeper.
    """
    exe = executable(binary)
    if search_mode not in SEARCH_MODES:
        raise ValueError(f"search_mode must be one of {SEARCH_MODES}, "
                         f"got {search_mode!r}")
    if receptor.pdbqt is None or not Path(receptor.pdbqt).is_file():
        raise FileNotFoundError(f"receptor PDBQT missing for {receptor.pdb_id}")

    with tempfile.TemporaryDirectory(prefix="unidock_") as tmp:
        tmp = Path(tmp)
        index = tmp / "ligands.txt"
        index.write_text(f"{Path(ligand_pdbqt).resolve()}\n", encoding="utf-8")
        out_dir = tmp / "out"
        out_dir.mkdir()
        cx, cy, cz = (float(v) for v in receptor.center)
        sx, sy, sz = (float(v) for v in receptor.size)
        cmd = [
            exe,
            "--receptor", str(Path(receptor.pdbqt).resolve()),
            "--ligand_index", str(index),
            "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}",
            "--center_z", f"{cz:.3f}",
            "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}",
            "--size_z", f"{sz:.3f}",
            "--search_mode", search_mode,
            "--num_modes", str(int(n_poses)),
            "--seed", str(int(seed)),
            "--dir", str(out_dir),
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=timeout_s)
        except FileNotFoundError as exc:
            raise UnidockUnavailable(f"could not execute {exe}: {exc}") from exc
        except subprocess.TimeoutExpired as exc:
            raise UnidockUnavailable(
                f"unidock did not finish within {timeout_s}s") from exc
        if proc.returncode != 0:
            raise UnidockUnavailable(
                f"unidock exited {proc.returncode}: "
                f"{(proc.stderr or proc.stdout or '')[:300]}")

        # Uni-Dock names its output <input stem>_out.pdbqt. Read whatever landed
        # rather than reconstructing the name, so a change in its convention
        # surfaces as "no output" instead of as a silent empty dock.
        produced = sorted(out_dir.glob("*.pdbqt"))
        if not produced:
            raise UnidockUnavailable(
                "unidock wrote no PDBQT; stdout: "
                f"{(proc.stdout or '')[-300:]}")
        return produced[0].read_text(encoding="utf-8")
