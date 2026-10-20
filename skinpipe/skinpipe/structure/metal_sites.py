'What happened to the catalytic metals, recorded where the next stage can see it.'
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Optional

from .rcsb import CATALYTIC_METALS

#: Appended to a receptor PDB's stem to name its metal sidecar:
#: ``2Y9X_A_clean.pdb`` -> ``2Y9X_A_clean.metals.pdb``.
SIDECAR_SUFFIX = ".metals.pdb"


def sidecar_path(receptor_pdb: Path) -> Path:
    """Where ``receptor_pdb``'s held-out metals are (or would be) recorded."""
    path = Path(receptor_pdb)
    return path.with_name(path.stem + SIDECAR_SUFFIX)


def metal_lines_in(pdb_text: str) -> list[str]:
    """Catalytic-metal coordinate records in ``pdb_text``, first model only.

    Alternate records follow the same occupancy selection as the receptor,
    including metals modelled only in B or another non-A location.
    """
    from .clean import atom_lines, select_altlocs

    out: list[str] = []
    for line in select_altlocs(atom_lines(pdb_text)):
        if line[17:20].strip().upper() in CATALYTIC_METALS:
            out.append(line)
    return out


def _dedup(lines: Iterable[str]) -> list[str]:
    """Drop repeats keyed on residue identity plus literal coordinates."""
    seen: set[str] = set()
    out: list[str] = []
    for line in lines:
        key = line[17:27] + line[30:54]
        if key in seen:
            continue
        seen.add(key)
        out.append(line)
    return out


def clear_sidecar(receptor_pdb: Path) -> bool:
    """Delete any metal record beside ``receptor_pdb``; True if one was there.

    The counterpart of :func:`write_sidecar`, and not merely a convenience: a
    sidecar outlives the file it describes unless something actively removes
    it, and the paths that MM/GBSA prepares into are reused. A caller that
    writes a metal-free receptor to a path a holo one used has to say so here,
    or :mod:`skinpipe.mmgbsa.metal_guard` goes on refusing a receptor that is
    fine, on the strength of a metal from a different structure.
    """
    path = sidecar_path(receptor_pdb)
    if not path.is_file():
        return False
    path.unlink(missing_ok=True)
    return True


def write_sidecar(receptor_pdb: Path, metal_lines: Iterable[str]) -> Optional[Path]:
    """Record the metals held out of ``receptor_pdb``. Returns the path, or None.

    Writing nothing when there are no metals is deliberate. An empty sidecar and
    an absent one would mean the same thing, and a file that is always present
    is a file nobody looks at. A *stale* sidecar is worse than no sidecar -- it
    would make :mod:`skinpipe.mmgbsa.metal_guard` refuse a receptor that has no
    metal -- so an existing one is deleted when the new record is empty.

    That deletion only helps if this function is reached on the empty path, so
    call it unconditionally, or call :func:`clear_sidecar` yourself. Guarding
    the call with ``if metals:`` is the bug it is written to prevent.
    """
    path = sidecar_path(receptor_pdb)
    lines = _dedup(line for line in metal_lines if line.strip())
    if not lines:
        clear_sidecar(receptor_pdb)
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"REMARK   1 catalytic metal atoms held out of {Path(receptor_pdb).name}\n"
        "REMARK   1 written by skinpipe.structure.metal_sites; read by\n"
        "REMARK   1 skinpipe.mmgbsa.metal_guard, which refuses to score a pose\n"
        "REMARK   1 that binds one of them.\n"
        + "\n".join(lines) + "\nEND\n"
    )
    return path


def read_sidecar(receptor_pdb: Path) -> list[str]:
    """Metal records held out of ``receptor_pdb``; empty when there is no sidecar."""
    path = sidecar_path(receptor_pdb)
    if not path.is_file():
        return []
    return metal_lines_in(path.read_text())


def all_metal_lines(receptor_pdb: Path) -> list[str]:
    """Metals still in the receptor file, plus any recorded beside it.

    Both halves are needed. A caller who hands MM/GBSA a raw holo PDB still has
    its metals in the file (and ``mmgbsa.system.prepare_receptor`` is about to
    delete them); the pipeline hands it a cleaned PDB whose metals live only in
    the sidecar. Either way the binding site has a metal in it.
    """
    path = Path(receptor_pdb)
    lines = metal_lines_in(path.read_text()) if path.is_file() else []
    return _dedup(lines + read_sidecar(path))


@dataclass(frozen=True)
class MetalAtom:
    """One catalytic metal atom and where it sits, in Angstrom."""

    resname: str
    chain: str
    resnum: str
    xyz: tuple[float, float, float]

    @property
    def label(self) -> str:
        return f"{self.resname} {self.chain}:{self.resnum}".strip()


def parse_metal_atoms(lines: Iterable[str]) -> list[MetalAtom]:
    """Turn metal coordinate records into :class:`MetalAtom` values.

    A record whose coordinates will not parse raises rather than being skipped.
    A metal we cannot locate is precisely the one that must not be quietly
    dropped: skipping it here would let the guard downstream conclude the
    binding site is metal-free, which is the failure this whole path exists to
    prevent.
    """
    out: list[MetalAtom] = []
    for line in lines:
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError("non-finite metal coordinates")
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"unparseable metal coordinate record: {line[:54]!r}") from exc
        out.append(MetalAtom(resname=line[17:20].strip().upper(), chain=line[21],
                             resnum=line[22:27].strip(), xyz=xyz))
    return out


def receptor_metals(receptor_pdb: Path) -> list[MetalAtom]:
    """Every catalytic metal belonging to ``receptor_pdb``, in file or sidecar."""
    return parse_metal_atoms(all_metal_lines(receptor_pdb))
