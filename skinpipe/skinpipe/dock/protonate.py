"""pH-aware ligand protonation.

Why this matters more than it looks
-----------------------------------
Docking a ligand in the wrong ionisation state is a third, distinct failure
mode, and it looks like neither of the other two. Redocking 4-hydroxytamoxifen
into 3ERT lands the top pose 4.85 A from the crystal, and raising exhaustiveness
from 16 to 64 moves the best pose from 2.57 to 2.59 A -- no improvement at four
times the search. The search is not the problem and the box is not the problem.

The problem is that the tertiary amine on the tail has pKa ~8.9 and is
protonated at skin pH, where it forms a salt bridge to Asp351. Prepared from
neutral SMILES, that charge is absent, the salt bridge cannot be scored, and
nothing anchors the tail. No amount of sampling recovers an interaction the
force field was never told about.

Skin actives are full of ionisable groups: azelaic acid is a diacid, tranexamic
acid is a zwitterion, aminophenols and alkaloids carry basic nitrogens. Docking
them neutral is simply wrong.

Method
------
Open Babel's rule-based `-p` transform applies functional-group pKa rules,
atom by atom, to each ionisable centre. The output is one rule-transformed
structure; it is not an experimentally validated dominant protomer or a
population-weighted ensemble. See https://openbabel.org/docs/Command-line_tools/babel.html.
When
Open Babel is unavailable the SMILES passes through untouched and the caller is
told, so a silently neutral screen is never mistaken for a titrated one.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
from functools import lru_cache
from typing import Any, Optional

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

_TIMEOUT_S = 30


@lru_cache(maxsize=1)
def obabel_path() -> Optional[str]:
    """Find the Open Babel binary.

    Look beside the running interpreter before consulting PATH. A conda or
    micromamba environment invoked by absolute path -- which is how this
    pipeline is usually run -- has its own ``bin`` directory populated but not
    exported, so a bare PATH lookup misses a binary that is sitting right there.
    """
    local = os.path.join(os.path.dirname(sys.executable), "obabel")
    if os.path.isfile(local) and os.access(local, os.X_OK):
        return local
    return shutil.which("obabel")


def available() -> bool:
    return obabel_path() is not None


def protonate_smiles(smiles: str, ph: float = 7.4,
                     execution: Optional[dict[str, Any]] = None) -> tuple[str, Optional[str]]:
    """Return ``(smiles_at_ph, note)``.

    ``note`` is ``None`` when the molecule was titrated and left unchanged;
    otherwise it says what happened, including the case where Open Babel is
    missing and the input was passed through as written.
    """
    trace = execution if execution is not None else {}
    trace.clear()
    trace.update(requested_ph=float(ph), input_smiles=smiles, applied=False)
    def finish(output: str, note: Optional[str], state: str, applied: bool = False):
        trace.update(output_smiles=output, note=note, state=state, applied=applied)
        return output, note
    exe = obabel_path()
    trace['binary_requested'] = exe
    if exe is None:
        return finish(smiles, ("Open Babel not found; requested pH transform was not applied. "
                        "Input returned unchanged; this does not qualify it for requested-pH docking"), 'unavailable')
    try:
        proc = subprocess.run(
            [exe, "-ismi", "-osmi", "-p", str(ph)],
            input=smiles, capture_output=True, text=True, timeout=_TIMEOUT_S)
    except (subprocess.TimeoutExpired, OSError) as exc:
        return finish(smiles, f"protonation failed ({type(exc).__name__}); used the SMILES as written",
                      'tool_error')

    trace['returncode'] = proc.returncode
    if proc.returncode != 0:
        return finish(smiles, f"Open Babel exited with code {proc.returncode}; used the SMILES as written",
                      'tool_error')

    out = (proc.stdout or "").strip().split("\n")[0].split("\t")[0].strip()
    if not out:
        return finish(smiles, "Open Babel returned nothing; used the SMILES as written", 'empty_output')

    mol = Chem.MolFromSmiles(out)
    if mol is None:
        return finish(smiles, "Open Babel produced an unparseable SMILES; used the input as written",
                      'invalid_output')

    before = Chem.MolFromSmiles(smiles)
    q_before = Chem.GetFormalCharge(before) if before else 0
    q_after = Chem.GetFormalCharge(mol)
    canonical = Chem.MolToSmiles(mol)

    if before is not None and Chem.MolToSmiles(before) == canonical:
        return finish(canonical, None, 'applied_unchanged', True)

    if q_before != q_after:
        return finish(canonical, f"titrated at pH {ph}: net charge {q_before:+d} -> {q_after:+d}",
                      'applied_changed', True)

    # A zwitterion gains a cation and an anion, so the net charge is unchanged
    # while the electrostatics are transformed. Tranexamic acid is exactly this,
    # and reporting "no change" would be misleading.
    n_charged = sum(1 for a in mol.GetAtoms() if a.GetFormalCharge())
    return finish(canonical, (f"titrated at pH {ph}: net charge unchanged ({q_after:+d}) but "
                       f"{n_charged} atom(s) now carry a formal charge"), 'applied_changed', True)
