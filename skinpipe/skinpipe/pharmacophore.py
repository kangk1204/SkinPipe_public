"""Anchor a catalog search on the part of the molecule that does the binding.

Why similarity search is not enough
-----------------------------------
A 2D fingerprint measures the whole molecule, and on a small one the bulk wins.
Screening PubChem by similarity to 4-butylresorcinol -- a clinically validated
human-tyrosinase inhibitor -- returned a run of simple alkylphenols: 2-butyl-,
4-pentyl-, 4-hexylphenol and their neighbours. They scored well because the
butylbenzene half of the query is most of its fingerprint, and dropping one
hydroxyl barely moves it.

Dropping that hydroxyl is the whole story. Tyrosinase is inhibited because the
resorcinol 1,3-diol chelates the binuclear copper. A phenol has nothing to
chelate with. The search retrieved the scaffold and discarded the pharmacophore,
and every hit was worthless in a way no downstream filter would catch, because
each one is a perfectly reasonable small molecule.

So the search is anchored instead. A substructure query on the core cannot
return a compound that lacks it.

Deriving the core
-----------------
Hand-declaring a core per target does not scale and encodes one person's guess.
`derive_core` takes several known actives and returns their maximum common
substructure, which is the part the chemistry has preserved across the series --
usually because it is the part that binds. On four validated tyrosinase
inhibitors (4-butylresorcinol, oxyresveratrol, phenylethyl resorcinol and a
resorcinyl ketone) it returns a nine-atom alkyl-substituted benzene-1,3-diol
that matches all four and rejects 4-hexylphenol.

The MCS is a hypothesis about the pharmacophore, not a measurement of it. It is
only as good as the series it is given: actives that share a scaffold by
historical accident will produce a core that encodes the accident. Feed it
chemically diverse actives, and read the core before trusting it.
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

#: What happened to the 3D overlay on a screened row, and where that puts it in
#: the ranking. This is the same policy `hop_core.hop_score` applies to
#: `hop_quality` -- a measurement that failed is not a measurement, and must
#: never be substituted for by whatever other number happens to be in the row.
#: See `_shape_sort_key` for why the two cases below "measured" are separate.
SHAPE_MEASURED = "measured"
SHAPE_NOT_REQUESTED = "not requested"     # no query_smiles, so nothing was tried
SHAPE_FAILED = "failed"                   # tried, and came back None or NaN



_SHAPE_TIER = {SHAPE_MEASURED: 0, SHAPE_NOT_REQUESTED: 1, SHAPE_FAILED: 2}

SHAPE_FAILURE_NOTE = (
    "3D shape could not be measured (embedding or MMFF failed), so this "
    "compound is ranked below every compound whose overlay succeeded. Its "
    "topical score is NOT a stand-in for the missing shape: the two are on "
    "different scales and substituting one for the other pays a compound for "
    "failing to be measured")

#: Cores for targets whose binding chemistry is settled, so a caller need not
#: supply a series. Each is the motif the literature attributes the binding to,
#: with the reason, because a core with no rationale is a guess someone will
#: later mistake for a fact.
KNOWN_CORES: dict[str, dict[str, str]] = {
    "TYR": {
        "smarts": "[#6]c1ccc([OX2H])cc1[OX2H]",
        "why": "resorcinol 1,3-diol chelates the binuclear copper; a single "
               "phenol cannot, which is why alkylphenols are false hits",
    },
    "MMP1": {
        "smarts": "[CX4,c][CX3](=[OX1])[NX3][OX2H]",
        "why": "hydroxamate chelates the catalytic zinc",
    },
    "MMP3": {
        "smarts": "[CX4,c][CX3](=[OX1])[NX3][OX2H]",
        "why": "hydroxamate chelates the catalytic zinc",
    },
}


@dataclass
class Core:
    smarts: str
    source: str                      # "declared" | "mcs" | "user"
    why: str = ""
    n_atoms: int = 0
    matched_actives: int = 0
    total_actives: int = 0
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return Chem.MolFromSmarts(self.smarts) is not None


def derive_core(actives: Iterable[str], *, ring_matches_ring_only: bool = True,
                complete_rings_only: bool = True, timeout_s: int = 30,
                min_atoms: int = 5) -> Optional[Core]:
    """Maximum common substructure of a series, as a pharmacophore hypothesis.

    Returns ``None`` when the series shares too little to be meaningful, which
    is itself informative: a core of three atoms is not a pharmacophore, it is
    a benzene ring, and anchoring a search on it would retrieve everything.
    """
    smis = [s for s in actives if s]
    mols = [m for m in (Chem.MolFromSmiles(s) for s in smis) if m is not None]
    if len(mols) < 2:
        log.warning("need at least two parseable actives to derive a core")
        return None
    res = rdFMCS.FindMCS(mols, ringMatchesRingOnly=ring_matches_ring_only,
                         completeRingsOnly=complete_rings_only, timeout=timeout_s)
    if not res.smartsString or res.numAtoms < min_atoms:
        log.warning("MCS of %d actives is only %d atoms; too small to anchor on",
                    len(mols), res.numAtoms)
        return None
    q = Chem.MolFromSmarts(res.smartsString)
    matched = sum(1 for m in mols if q is not None and m.HasSubstructMatch(q))
    core = Core(smarts=res.smartsString, source="mcs", n_atoms=res.numAtoms,
                matched_actives=matched, total_actives=len(mols),
                why=f"maximum common substructure of {len(mols)} actives")
    if res.canceled:
        core.warnings.append(
            f"the MCS search hit its {timeout_s}s timeout, so this core may be "
            "smaller than the true common substructure")
    if matched < len(mols):
        core.warnings.append(
            f"the derived core matches only {matched} of {len(mols)} actives, "
            "which should not happen; inspect the series for a parse failure")
    if res.numAtoms < 8:
        core.warnings.append(
            f"{res.numAtoms} atoms is a small core. Check it is a pharmacophore "
            "and not just a shared ring before screening on it")
    return core


def core_for(gene: Optional[str] = None, actives: Optional[Iterable[str]] = None,
             smarts: Optional[str] = None) -> Optional[Core]:
    """Resolve a core from an explicit SMARTS, a known target, or a series."""
    if smarts:
        c = Core(smarts=smarts, source="user", why="supplied by the caller")
        q = Chem.MolFromSmarts(smarts)
        c.n_atoms = q.GetNumAtoms() if q is not None else 0
        if q is None:
            c.warnings.append("the supplied SMARTS does not compile")
        return c
    if gene and gene.upper() in KNOWN_CORES:
        d = KNOWN_CORES[gene.upper()]
        q = Chem.MolFromSmarts(d["smarts"])
        return Core(smarts=d["smarts"], source="declared", why=d["why"],
                    n_atoms=q.GetNumAtoms() if q is not None else 0)
    if actives:
        return derive_core(actives)
    return None


def substructure_cids(core: Core, max_records: int = 500,
                      session: Optional[Any] = None) -> list[int]:
    """CIDs containing the core, queried as SMARTS rather than SMILES.

    PubChem exposes a separate `fastsubstructure/smarts` route, and the SMARTS
    must be form-encoded: put it in the URL path unescaped and the brackets and
    braces break the request with a 400 that reads like a malformed query. A
    POST carries it cleanly and has no length limit, which matters because a
    derived MCS core is often long.
    """
    import requests

    from . import catalog

    s = session or requests.Session()
    r = catalog._request(
        s, "POST",
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smarts/cids/JSON",
        data={"smarts": core.smarts, "MaxRecords": max_records})
    if r is None:
        return []
    # `catalog._json` raises PubChemBusy when the body is not JSON at all: a
    # 200 carrying a proxy error page is a transport failure, and returning []
    # for it would report "no compound contains this core", which is a claim
    # about chemistry made from a network fault. See catalog._json.
    body = catalog._json(r, "the anchored substructure search")
    return catalog.parse_cids(body, "the anchored substructure search")


def has_core(smiles: str, core: Core) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    q = Chem.MolFromSmarts(core.smarts)
    return bool(mol is not None and q is not None and mol.HasSubstructMatch(q))


def _finite_shape(value: Any) -> Optional[float]:
    'The overlay result, or None if it is not a usable number.'
    if value is None:
        return None
    try:
        shape = float(value)
    except (TypeError, ValueError):
        return None
    return shape if math.isfinite(shape) else None


def _shape_sort_key(row: dict[str, Any]) -> tuple[int, float, float]:
    """Rank a screened row: tier first, then shape, then topical fitness.

    The bug this replaces put the shape Tanimoto (0-1) and the topical score
    (0-10) into one key, choosing between them on whether the overlay had
    succeeded. Since the key was negated for a descending sort, a failed
    overlay with topical 10.00 scored -10.00 against a successful overlay's
    -0.90, and the failure sorted first. A compound was promoted to the top of
    the list *because* its 3D measurement had failed, and the better its skin
    score the higher it went.

    So the two scales never meet. The tier decides first and no other column
    can overturn it; the shape orders the measured rows against each other; the
    topical score only ever breaks a tie inside one tier. When no
    `query_smiles` was given, nothing was attempted and nothing failed, so
    every row shares the "not requested" tier and the ranking falls through to
    topical fitness for all of them -- the documented behaviour of that mode.
    """
    status = row["shape_status"]
    shape = row["shape_tanimoto"] if status == SHAPE_MEASURED else None
    return (_SHAPE_TIER[status],
            -(shape if shape is not None else 0.0),
            -(row["topical_score"] or 0.0))


def screen(core: Core, *, query_smiles: Optional[str] = None,
           max_records: int = 500, min_topical: float = 6.0,
           min_vendors: int = 1, vendor_check_top_n: int = 120,
           shape_top_n: int = 60, n_confs: int = 10,
           admet_gate: Optional[str] = "critical",
           target_gene: Optional[str] = None,
           progress: bool = True) -> list[dict[str, Any]]:
    'Catalog funnel anchored on a pharmacophore core.'
    from . import catalog
    from . import hopping
    from .filters import admet as admet_f
    from .filters import topical as topical_f
    from .types import LigandRecord

    if not core.valid:
        raise ValueError(f"core SMARTS does not compile: {core.smarts!r}")
    admet_gate = catalog.validate_gate(admet_gate)
    for w in core.warnings:
        log.warning("core: %s", w)

    cids = substructure_cids(core, max_records=max_records)
    if progress:
        log.info("substructure on core (%d atoms): %d compounds", core.n_atoms, len(cids))
    if not cids:
        return catalog.empty_funnel(cids)

    props, prop_meta = catalog.fetch_properties(cids)
    if prop_meta["unanswered"]:
        log.warning("%d of %d compounds received no property answer; unchecked, not rejected by the skin filter",
                    len(prop_meta["unanswered"]), prop_meta["requested"])
    if not props:
        raise catalog.PubChemBusy(
            "the substructure search returned compounds but no properties came "
            "back; that is throttling, not an empty result set")

    kept: dict[int, dict] = {}
    unscored = []
    for cid, p in props.items():
        if p["mw"] is None or p["logp"] is None:
            unscored.append(cid)
            continue
        if not has_core(p["smiles"], core):        # PubChem's match is looser than RDKit's
            continue
        a = topical_f.assess(LigandRecord(ligand_id=str(cid), smiles=p["smiles"], props=p))
        if a.topical_score >= min_topical:
            kept[cid] = p | {"topical_score": a.topical_score, "log_kp": a.log_kp,
                             "permeation": a.permeation_class}
    if progress:
        log.info("core confirmed + skin filter: %d/%d", len(kept), len(props) - len(unscored))
    if unscored:
        log.warning("%d compounds have no MW/XLogP and were not scored against the skin filter", len(unscored))
    if not kept:
        return catalog.empty_funnel(cids)









    shortlist = dict(sorted(kept.items(), key=lambda kv: -kv[1]["topical_score"]
                            )[:vendor_check_top_n])
    if progress and len(shortlist) < len(kept):
        log.info("vendor check limited to the top %d by topical fitness (of %d)",
                 len(shortlist), len(kept))
    vendors, vendor_meta = catalog.vendor_lookup(shortlist)
    # Same distinction `catalog.screen` makes: a compound the lookup never
    # reached is unchecked, not unsold, and must not be cut as though it were.
    unchecked = set(vendor_meta["unanswered"])
    if unchecked:
        log.warning("%d of %d vendor lookups did not answer; those compounds are "
                    "unchecked, not unpurchasable",
                    len(unchecked), vendor_meta["requested"])
    buyable = {c: d | {"vendors": vendors.get(c, []),
                       "n_vendors": len(vendors.get(c, [])),
                       "vendor_status": "unchecked" if c in unchecked else "checked"}
               for c, d in shortlist.items()
               if c not in unchecked and len(vendors.get(c, [])) >= min_vendors}
    if progress:
        log.info("vendors: %d/%d purchasable", len(buyable), len(kept))
    if not buyable:
        return catalog.empty_funnel(cids)

    ranked = sorted(buyable.items(), key=lambda kv: -kv[1]["topical_score"])[:shape_top_n]
    out: list[dict[str, Any]] = []
    for cid, d in ranked:
        shape = _finite_shape(
            hopping.best_shape_from_smiles(query_smiles, d["smiles"], n_confs=n_confs)
            if query_smiles else None)
        if not query_smiles:
            status = SHAPE_NOT_REQUESTED
        elif shape is None:
            status = SHAPE_FAILED
        else:
            status = SHAPE_MEASURED
        # The target decides whether a hazard the target itself activates is
        # read as a hazard at all; see `catalog.screen`.
        prof = admet_f.summarise(d["smiles"], ligand_id=str(cid),
                                 target_gene=target_gene)
        hazard = catalog._hazard(prof)
        if admet_gate and catalog._at_least(hazard, admet_gate):
            continue
        out.append({"cid": cid, "smiles": d["smiles"], "mw": d["mw"], "logp": d["logp"],
                    "log_kp": d.get("log_kp"), "permeation": d.get("permeation"),
                    "topical_score": d.get("topical_score"),
                    "shape_tanimoto": round(shape, 3) if shape is not None else None,
                    "shape_measured": status == SHAPE_MEASURED,
                    "shape_status": status,
                    "shape_note": SHAPE_FAILURE_NOTE if status == SHAPE_FAILED else "",
                    "n_vendors": d["n_vendors"], "vendors": d["vendors"][:4],
                    # Same three columns `catalog.screen` carries, and for the
                    # same reason: a hazard severity is unreadable without the
                    # target it was judged against, and the gate a row survived
                    # has to round-trip when the file is fed back in.
                    "target_gene": target_gene or "",
                    "admet_gate": admet_gate or "none",
                    "admet_hazard": hazard, "admet_overall": prof.overall,
                    "core": core.smarts})
    out.sort(key=_shape_sort_key)
    n_failed = sum(1 for r in out if r["shape_status"] == SHAPE_FAILED)
    if n_failed:
        log.warning("3D overlay failed for %d of %d survivors; they are ranked "
                    "last as unmeasured, not scored on topical fitness instead",
                    n_failed, len(out))
    if progress:
        log.info("shape + ADMET: %d survivors", len(out))
    return out
