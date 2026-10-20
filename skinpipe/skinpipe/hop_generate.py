"""Analogue generation: fragmentation, BRICS recombination, bioisostere edits.

Split out of `hopping.py`, which re-exports everything here and carries the
honesty contract for the whole feature -- read that module docstring first,
especially the part about how much of this output is junk.
"""
from __future__ import annotations

import contextlib
import logging
import random
from typing import Iterable, Iterator, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS, Recap

from .hop_bioisostere import bioisostere_replace
from .hop_sanity import SA_CUTOFF_DEFAULT, filter_generated, is_sane
from .library import standardize

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

BRICS_SEED = 0xC05


@contextlib.contextmanager
def _seeded(seed: int) -> Iterator[None]:
    """Seed the global RNG for a block, then put it back exactly as found.

    `BRICS.BRICSBuild` shuffles its reagents with the `random` module. Seeding
    it is the only way to make recombination reproducible, and restoring the
    previous state afterwards is the only way to do that without silently
    changing the behaviour of everything else in the process.
    """
    state = random.getstate()
    random.seed(seed)
    try:
        yield
    finally:
        random.setstate(state)


def fragment(smiles: str, method: str = "brics") -> list[str]:
    """Cut a molecule into the pieces a medicinal chemist would recognise.

    `method="brics"` (Degen et al. 2008, ChemMedChem 3:1503) breaks 16 bond
    types derived from a retrosynthetic analysis of known drugs, and labels
    each cut with the environment it came from (`[16*]`, `[8*]`...), so pieces
    can only be recombined where the chemistry would allow it. Those labels are
    what makes `enumerate_hops` produce assemblable molecules rather than
    arbitrary graphs, so they are kept in the returned SMILES.

    `method="recap"` (Lewell et al. 1998, JCICS 38:511) breaks a smaller set of
    genuinely common reaction bonds -- amide, ester, amine, ether, urea and a
    few more -- and returns the leaves of the cleavage tree. It cuts less, so
    it yields larger and more realistic fragments, but it also fails completely
    on molecules whose bonds are all outside its set: kojic acid returns zero
    RECAP leaves, while BRICS gives two fragments.

    Returns a sorted list so the output is deterministic.

    >>> fragment("Oc1ccc(/C=C/c2cc(O)cc(O)c2)c(O)c1")      # oxyresveratrol
    ['[16*]c1cc(O)cc(O)c1', '[16*]c1ccc(O)cc1O', '[7*]C[8*]']
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.warning("fragment: unparseable SMILES %r", smiles[:80])
        return []

    method = method.lower()
    if method == "brics":
        return sorted(BRICS.BRICSDecompose(mol))
    if method == "recap":
        leaves = sorted(Recap.RecapDecompose(mol).GetLeaves().keys())
        if not leaves:
            log.info("RECAP found no cleavable bond in %s; try method='brics'",
                     smiles[:60])
        return leaves
    raise ValueError(f"unknown fragmentation method {method!r}; use 'brics' or 'recap'")


# --------------------------------------------------------------- generation

def _brics_candidates(seed_smiles: str, n: int, fragment_pool: Optional[Iterable[str]],
                      max_depth: int, rng_seed: int) -> list[tuple[str, str]]:
    """BRICS recombination. Returns (smiles, derivation) pairs, unfiltered."""
    mol = Chem.MolFromSmiles(seed_smiles)
    if mol is None:
        return []
    pieces = sorted(BRICS.BRICSDecompose(mol))
    if fragment_pool:
        pieces = sorted(set(pieces) | set(fragment_pool))
    if len(pieces) < 2:
        log.info("BRICS gave %d fragment(s) for %s; nothing to recombine. This is "
                 "expected for small rigid molecules -- use the bioisostere route",
                 len(pieces), seed_smiles[:60])
        return []

    frag_mols = [m for m in (Chem.MolFromSmiles(p) for p in pieces) if m is not None]
    out: list[tuple[str, str]] = []
    with _seeded(rng_seed):
        builder = BRICS.BRICSBuild(frag_mols, scrambleReagents=True, maxDepth=max_depth)
        for i, product in enumerate(builder):
            if i >= n:
                break
            try:
                product.UpdatePropertyCache(strict=False)
                out.append((Chem.MolToSmiles(product), "brics_recombination"))
            except Exception:
                continue
    log.info("BRICS: %d fragments -> %d raw products (depth %d)",
             len(frag_mols), len(out), max_depth)
    return out


def _bioisostere_candidates(seed_smiles: str) -> list[tuple[str, str]]:
    """Single-edit bioisosteric replacement. Returns (smiles, derivation) pairs."""
    return [(h["smiles"], f"bioisostere:{h['rule']}")
            for h in bioisostere_replace(seed_smiles)]


def enumerate_hops(seed_smiles: str, n: int = 200,
                   methods: tuple[str, ...] = ("bioisostere", "brics"),
                   fragment_pool: Optional[Iterable[str]] = None,
                   max_depth: int = 2, max_sa: float = SA_CUTOFF_DEFAULT,
                   rng_seed: int = BRICS_SEED,
                   return_census: bool = False
                   ) -> list[str] | tuple[list[str], dict]:
    """Generate candidate analogues of `seed_smiles`, sanity- and SA-filtered.

    `n` caps the *raw* BRICS products pulled off the generator, before
    filtering, because the generator is effectively infinite and the cap is
    what bounds combinatorial explosion. The bioisostere route is not capped:
    it is bounded by the size of the rule table and terminates on its own.
    Expect far fewer than `n` survivors -- see `return_census`.

    `max_depth` bounds BRICS recombination depth. 2 means "at most two
    recombination steps from a fragment", which already covers most of the
    interesting rearrangements; 3 multiplies the count and the junk rate
    together.

    `fragment_pool` accepts extra BRICS-labelled fragments (from `fragment()`
    on other molecules) to recombine alongside the seed's own. This is how you
    graft the seed's recognition fragment onto a core from a different series.
    Fragments must carry BRICS dummy labels or they will not connect.

    The parent is always excluded from the result -- several bioisostere rules
    can regenerate it on symmetric substrates.

    With `return_census=True` returns `(smiles, census)` where census counts
    what was discarded and why. Use it. The yield of this function is part of
    its result, and a hop list quoted without it is not interpretable.

    The census also carries `census["derivations"]`, mapping each surviving
    SMILES to the route that made it (`"brics_recombination"`, or
    `"bioisostere:<rule name>"`). Pass it to `to_ligand_records` so every
    generated structure reaches the rest of the pipeline knowing where it came
    from.

    >>> hops = enumerate_hops("OCC1=CC(=O)C(O)=CO1", n=50)   # kojic acid
    >>> len(hops) > 0
    True
    """
    seed = standardize(seed_smiles)
    if seed is None:
        raise ValueError(f"could not parse seed SMILES {seed_smiles!r}")

    raw: list[tuple[str, str]] = []
    for method in methods:
        method = method.lower()
        if method == "bioisostere":
            raw += _bioisostere_candidates(seed)
        elif method == "brics":
            raw += _brics_candidates(seed, n, fragment_pool, max_depth, rng_seed)
        else:
            raise ValueError(f"unknown method {method!r}; use 'bioisostere' or 'brics'")

    by_origin: dict[str, str] = {}
    for smi, derivation in raw:
        by_origin.setdefault(smi, derivation)

    kept, census = filter_generated(list(by_origin), max_sa=max_sa)
    kept = [s for s in kept if s != seed]
    census["parent_regenerated"] = census["kept"] - len(kept)
    census["kept"] = len(kept)
    census["methods"] = ",".join(methods)
    # Which route and, for a bioisostere, which rule produced each survivor.
    # This is what lets `to_ligand_records` write a real `derivation` onto the
    # LigandRecord instead of a generic label, so a structure that reaches a
    # chemist can be traced back to the citation that justified it.
    canonical_origins: dict[str, str] = {}
    proposals = []
    for ordinal, (smi, derivation) in enumerate(raw, 1):
        mol = Chem.MolFromSmiles(smi)
        canonical = Chem.MolToSmiles(mol) if mol is not None else None
        if canonical is not None:
            canonical_origins.setdefault(canonical, derivation)
        if return_census:
            sane, reason = is_sane(smi)
            proposals.append({"ordinal": ordinal, "raw_smiles": smi,
                              "canonical_smiles": canonical,
                              "derivation": derivation,
                              "sanity_policy_pass": sane,
                              "sanity_policy_reason": reason or None,
                              "legacy_accepted": canonical in kept,
                              "regenerated_parent": canonical == seed})
    census["derivations"] = {s: canonical_origins[s] for s in kept}
    if return_census:
        # Keep the legacy survivor list for existing callers, while exposing
        # the proposal universe before textual deduplication and policy gates.
        # A policy rejection (e.g. ring size) is not an invalid structure.
        census["generation_parent_smiles"] = seed
        census["raw_proposal_count"] = len(raw)
        census["proposals"] = proposals

    log.info("enumerate_hops(%s): %d raw -> %d usable (%.0f%% yield)",
             seed[:50], census["input"], len(kept),
             100.0 * len(kept) / census["input"] if census["input"] else 0.0)
    return (kept, census) if return_census else kept
