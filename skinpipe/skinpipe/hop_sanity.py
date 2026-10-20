"""Chemical sanity and synthetic accessibility for generated molecules.

Why this module exists
----------------------
Everything in `hopping.py` that *generates* structures produces junk at a high
rate. BRICS recombination in particular will happily bolt three fragments
together into something with eleven rotatable bonds, a molecular weight of 800
and no synthetic route at any price. Reporting a hop list without saying how
much was thrown away to get it is the single easiest way to make a generative
method look better than it is, so this module exists to make the discard
counts explicit and to put them in front of the reader.

Two independent gates
---------------------
1.  **Chemical sanity** -- is this a molecule at all? RDKit's sanitiser
    catches impossible valences, but it happily accepts things no chemist
    would draw: leftover BRICS attachment dummies, radicals, nine-membered
    carbocycles, tungsten. `is_sane` checks those explicitly and returns the
    reason, because "43 molecules were dropped" is much less useful than
    "43 molecules were dropped, 39 of them for retained attachment points".

2.  **Synthetic accessibility** -- could anyone make it? The SA_Score of Ertl
    & Schuffenhauer (2009, J Cheminform 1:8) scores 1 (trivial) to 10
    (hopeless) by combining a fragment-frequency term, learned from the ~1M
    most common ECFP4 fragments in PubChem, with a complexity penalty for
    ring fusion, stereocentres and macrocycles. RDKit ships the reference
    implementation in its Contrib tree; this module imports it from
    `RDConfig.RDContribDir/SA_Score` and falls back to a documented heuristic
    if the Contrib tree is missing from the install.

What SA_Score is not
--------------------
It is a similarity-to-known-chemistry score, not a retrosynthesis. It has no
model of reactions, protecting groups or commercial availability. A molecule
made entirely of common fragments in an uncommon arrangement scores well and
may still have no route; a genuinely easy molecule built from rare fragments
scores badly. Treat the default cutoff of 6.0 as a coarse filter that removes
the obviously hopeless, and get a chemist to look at anything you intend to
order. For reference, on this project's own seeds: kojic acid 2.71,
oxyresveratrol 2.40, imatinib 2.33 -- real drugs and actives all sit low, which
is the calibration that makes a cutoff of 6.0 permissive rather than strict.
"""
from __future__ import annotations

import logging
import math
import os
import sys
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, RDConfig

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

# Elements a topical small molecule may plausibly contain. Deliberately tight:
# BRICS recombination cannot invent elements, but a user-supplied fragment pool
# can, and a boron or a selenium in a leave-on cosmetic is a regulatory
# conversation nobody wants to have by accident.
_ALLOWED_ELEMENTS = frozenset({1, 6, 7, 8, 9, 15, 16, 17, 35, 53})

_MAX_RING_SIZE = 7          # tropolone is 7; anything larger is a macrocycle
_MIN_HEAVY, _MAX_HEAVY = 6, 50
_MAX_ROTATABLE = 12

SA_CUTOFF_DEFAULT = 6.0     # Ertl & Schuffenhauer scale, 1 (easy) .. 10 (hopeless)


def _load_sascorer():
    """Import RDKit's reference SA_Score implementation from the Contrib tree.

    Contrib is not a package, so it is not importable without a path hack.
    This is the documented way to reach it and is what every published use of
    SA_Score does.
    """
    path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if not os.path.isdir(path):                          # pragma: no cover
        log.warning("RDKit Contrib/SA_Score not found at %s; "
                    "falling back to the heuristic complexity score", path)
        return None
    if path not in sys.path:
        sys.path.append(path)
    try:
        import sascorer                                  # type: ignore
        return sascorer
    except Exception as exc:                             # pragma: no cover
        log.warning("could not import sascorer (%s); using the heuristic fallback", exc)
        return None


_SASCORER = _load_sascorer()
USING_REFERENCE_SA = _SASCORER is not None


def _heuristic_sa(mol: Chem.Mol) -> float:
    """Fallback complexity score on the same 1-10 scale, used only if the
    RDKit Contrib tree is absent.

    This reproduces the *shape* of SA_Score's complexity penalty -- ring
    fusion, stereocentres, spiro atoms, macrocycles and sheer size -- but has
    no fragment-frequency term, which is the half of SA_Score that actually
    encodes "chemists have made things like this before". It is therefore
    strictly worse, and `USING_REFERENCE_SA` is exported so a caller can say
    so in a report instead of silently presenting the two as equivalent.
    """
    ri = mol.GetRingInfo()
    n_heavy = max(mol.GetNumHeavyAtoms(), 1)
    n_rings = ri.NumRings()
    bridge = Descriptors.rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    spiro = Descriptors.rdMolDescriptors.CalcNumSpiroAtoms(mol)
    stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
    macro = sum(1 for r in ri.AtomRings() if len(r) > 7)
    penalty = (math.log10(n_heavy) * 0.8 + 0.5 * bridge + 0.6 * spiro
               + 0.3 * stereo + 1.2 * macro + 0.2 * max(0, n_rings - 2))
    return round(min(10.0, max(1.0, 1.0 + penalty)), 2)


def sa_score(smiles: str) -> Optional[float]:
    """Ertl & Schuffenhauer synthetic accessibility, 1 (easy) to 10 (hopeless).

    Returns None for an unparseable SMILES rather than a sentinel number, so a
    parse failure can never be mistaken for an easy molecule.

    >>> round(sa_score("OCC1=CC(=O)C(O)=CO1"), 1)      # kojic acid
    2.7
    """
    return sa_assessment(smiles)["score"]


def sa_assessment(smiles: str) -> dict:
    """Record the method actually used for this molecule, keeping its score.

    A reference module can import successfully while its fragment model is
    unavailable at calculation time. Module availability alone cannot label a
    per-molecule fallback as the reference SA_Score.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"score": None, "method": None, "status": "NOT_EVALUATED",
                "fallback_reason": None}
    reason = "reference implementation unavailable"
    try:
        if _SASCORER is not None:
            value = round(float(_SASCORER.calculateScore(mol)), 2)
            valid = math.isfinite(value) and 1 <= value <= 10
            return {"score": value if valid else None, "method": "Ertl SA_Score",
                    "status": "EVALUATED" if valid else "NOT_EVALUATED",
                    "fallback_reason": None}
    except Exception as exc:                             # pragma: no cover
        log.debug("sascorer failed on %s (%s); using the heuristic", smiles[:60], exc)
        reason = f"reference calculation failed: {type(exc).__name__}"
    return {"score": _heuristic_sa(mol), "method": "fallback heuristic",
            "status": "EVALUATED", "fallback_reason": reason}


def is_sane(smiles: str) -> tuple[bool, str]:
    """Would a chemist accept this as a drawable, orderable molecule?

    Returns (ok, reason). `reason` is "" when ok, otherwise a short tag that
    can be counted, so callers can report *why* a generation run lost most of
    its output rather than just how much.

    The checks, and why each one is here:

    *   ``unparseable`` / ``valence`` -- RDKit could not build or sanitise it.
        BRICS recombination does produce these.
    *   ``attachment_point`` -- a ``[*]`` dummy survived. This is the most
        common BRICS failure and the most misleading, because the SMILES looks
        fine and every descriptor computes.
    *   ``multi_fragment`` -- two disconnected pieces. `library.standardize`
        would silently keep only the larger one, which is not what the
        generator meant to make.
    *   ``radical`` -- unpaired electrons, usually from a botched bond split.
    *   ``element`` -- outside the organic set a topical product may contain.
    *   ``ring_size`` -- a ring larger than 7. Macrocycles are legitimate
        chemistry but are not what a recombination of drug fragments should be
        producing, and they will not cross the stratum corneum anyway.
    *   ``too_small`` / ``too_large`` / ``too_floppy`` -- outside the range
        where the rest of this pipeline has any validity. The 500 Da topical
        ceiling is applied later by `filters.topical`, not here; this is only
        a sanity bound.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "unparseable"
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return False, "valence"

    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        return False, "attachment_point"
    if len(Chem.GetMolFrags(mol)) > 1:
        return False, "multi_fragment"
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return False, "radical"
    if any(a.GetAtomicNum() not in _ALLOWED_ELEMENTS for a in mol.GetAtoms()):
        return False, "element"
    if any(len(r) > _MAX_RING_SIZE for r in mol.GetRingInfo().AtomRings()):
        return False, "ring_size"

    n_heavy = mol.GetNumHeavyAtoms()
    if n_heavy < _MIN_HEAVY:
        return False, "too_small"
    if n_heavy > _MAX_HEAVY:
        return False, "too_large"
    if Descriptors.NumRotatableBonds(mol) > _MAX_ROTATABLE:
        return False, "too_floppy"
    return True, ""


def filter_generated(smiles_list: list[str], max_sa: float = SA_CUTOFF_DEFAULT
                     ) -> tuple[list[str], dict[str, int]]:
    """Apply both gates and return the survivors plus a full discard census.

    The census is returned, never merely logged, because the yield of a
    generation run is a result in its own right. A method that produces 300
    structures and keeps 12 is not the same method as one that keeps 250, and
    the reader cannot tell them apart from the hit list alone.

    Counts in the returned dict: ``input``, ``kept``, ``duplicate``, one key
    per `is_sane` failure tag, and ``sa_reject`` for molecules that were
    chemically fine but scored above `max_sa`.
    """
    census: dict[str, int] = {"input": len(smiles_list), "kept": 0, "duplicate": 0}
    seen: set[str] = set()
    kept: list[str] = []

    for smi in smiles_list:
        ok, reason = is_sane(smi)
        if not ok:
            census[reason] = census.get(reason, 0) + 1
            continue
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        if canonical in seen:
            census["duplicate"] += 1
            continue
        seen.add(canonical)
        score = sa_score(canonical)
        if score is None or not math.isfinite(score) or not 1 <= score <= 10 or score > max_sa:
            census["sa_reject"] = census.get("sa_reject", 0) + 1
            continue
        kept.append(canonical)

    census["kept"] = len(kept)
    log.info("sanity/SA filter: %d in -> %d kept (%s)", census["input"], census["kept"],
             ", ".join(f"{k}={v}" for k, v in sorted(census.items())
                       if k not in ("input", "kept") and v))
    return kept, census
