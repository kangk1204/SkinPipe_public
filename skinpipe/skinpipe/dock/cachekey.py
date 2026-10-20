'Everything that can change a docking answer, in one comparable record.'
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Optional

from ..types import Receptor







CONDITIONS_VERSION = 5

_HASH_CHUNK = 1 << 20


def file_sha256(path: Optional[Path]) -> Optional[str]:
    """First 16 hex digits of a file's SHA-256; ``None`` when it is not there.

    ``None`` is deliberately not a digest of anything, so a receptor that has
    been deleted never compares equal to one that is present: a cache entry
    written against a file that has since vanished is a mismatch, not a match.
    """
    if path is None:
        return None
    path = Path(path)
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()[:16]


#: PDB element symbols and AutoDock types that mean hydrogen. ``HG`` is absent
#: on purpose: in an element column that is mercury.
_HYDROGEN = frozenset({"H", "HD", "HS"})


def _hydrogen_key(name: str) -> str:
    """A hydrogen's name in one spelling, so two writers' files still match.

    PDB v2 puts a branch number in front (``1HB``) where v3 puts it at the end
    (``HB1``); they name the same atom. Rotating the leading digit to the back
    makes the two spellings compare equal. It cannot merge two *different*
    hydrogens: ``HD1`` and ``HE2`` are untouched by it.
    """
    name = name.strip().upper()
    if name[:1].isdigit():
        name = name[1:] + name[0]
    return name


#: Separates the heavy-atom rows from the hydrogen-topology rows in the hashed
#: payload, so a residue line can never be read as an atom line. Heavy rows
#: always contain "|"; this marker never does.
_HYDROGEN_SECTION = "--hydrogen topology--"


def structure_sha256(path: Optional[Path]) -> Optional[str]:
    "Hash of a structure's chemistry, deliberately not of its bytes."
    if path is None or not Path(path).is_file():
        return None
    heavy: list[str] = []
    hydrogens: dict[str, list[str]] = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        kind = line[76:79].strip().upper()
        name = line[12:16]
        # No element column (older hand-made PDBs): fall back to the atom-name
        # convention, where a hydrogen's name starts with H or with a digit.
        is_h = kind in _HYDROGEN if kind else name.strip()[:1] in ("H", "1", "2", "3")
        if is_h:
            # line[17:27] is residue name, chain, number and insertion code:
            # which residue owns this proton. The AutoDock type is left out --
            # HD versus H is a typing choice, not a change of chemistry.
            hydrogens.setdefault(line[17:27], []).append(_hydrogen_key(name))
            continue
        heavy.append(f"{name}|{line[17:27]}|{line[30:54]}|{kind}")
    if not heavy and not hydrogens:
        return file_sha256(path)
    rows = heavy + [_HYDROGEN_SECTION]
    rows += [f"{residue}|{' '.join(sorted(names))}"
             for residue, names in sorted(hydrogens.items())]
    return hashlib.sha256("\n".join(rows).encode()).hexdigest()[:16]


def dock_conditions(receptor: Receptor, exhaustiveness: int, n_poses: int,
                    seed: int, engine: str = "vina",
                    search_mode: Optional[str] = None,
                    unidock_bin: Optional[str] = None) -> dict[str, Any]:
    """The full set of inputs a cached pose list is only valid for.

    Vina receives the exact float box arguments. Uni-Dock receives formatted
    three-decimal arguments, so its effective values are recorded at that
    precision. Lists preserve JSON round-trip equality without rounding Vina.
    """
    from . import software
    cond: dict[str, Any] = {
        "version": CONDITIONS_VERSION,
        "pdb_id": receptor.pdb_id,
        "chain": receptor.chain,
        "receptor_atoms_sha256": structure_sha256(receptor.pdbqt),
        "receptor_pdbqt_sha256": software.file_identity(receptor.pdbqt),
        "center": [float(f"{float(v):.3f}") if engine == "unidock" else float(v)
                   for v in receptor.center],
        "size": [float(f"{float(v):.3f}") if engine == "unidock" else float(v)
                 for v in receptor.size],
        "exhaustiveness": int(exhaustiveness),
        "n_poses": int(n_poses),
        "seed": int(seed),







        "engine": str(engine),
        "software_identity": software.engine_identity(engine, unidock_bin),
    }
    # Uni-Dock's `fast`, `balance` and `detail` are different searches, and
    # nothing above distinguishes them: `exhaustiveness` is the Vina knob and
    # is not what `--search_mode` sets. Without this, a `detail` run in a
    # workdir that already holds `fast` results is served the `fast` numbers
    # under the `detail` heading -- the same failure the `engine` key exists
    # to prevent, one level down. Added only for the engine that has the
    # setting, so no Vina cache is invalidated by it.
    if engine == "unidock":
        cond["search_mode"] = str(search_mode or "balance")
        cond["execution_mode"] = "single_ligand"
    return cond


def unidock_batch_conditions(receptor: Receptor, ligands: list[Any],
                             exhaustiveness: int, n_poses: int, seed: int,
                             search_mode: str = "balance",
                             unidock_bin: Optional[str] = None) -> dict[str, Any]:
    """Bind GPU-batch reuse to its ordered input library and batching policy."""
    from . import unidock
    conditions = dock_conditions(receptor, exhaustiveness, n_poses, seed,
                                 engine="unidock", search_mode=search_mode,
                                 unidock_bin=unidock_bin)
    inputs = [(ligand.ligand_id, ligand.smiles) for ligand in ligands]
    conditions.update(execution_mode="gpu_batch", batch_size=unidock.BATCH_SIZE,
                      batch_policy="fixed_full_eligible_library_v1",
                      ordered_library_sha256=hashlib.sha256(json.dumps(
                          inputs, separators=(',', ':'), ensure_ascii=True).encode()).hexdigest())
    return conditions


def describe_mismatch(cached: Any, current: dict[str, Any]) -> list[str]:
    'Every condition that differs, phrased as ``field was -> now``.'
    if not isinstance(cached, dict):
        return ["written before docking conditions were recorded"]
    if cached.get("receptor_atoms_sha256") is None:
        return ["the receptor PDBQT it was docked against is no longer on disk"]
    changed = []
    for field in sorted(set(cached) | set(current)):
        was = cached.get(field, "<absent>")
        now = current.get(field, "<absent>")
        if was != now:
            changed.append(f"{field} {was!r} -> {now!r}")
    return changed
