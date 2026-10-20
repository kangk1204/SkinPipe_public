'PDB ID in, docking-ready :class:`~skinpipe.types.Receptor` out.'
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Optional

import numpy as np

from ..types import Receptor
from . import assembly, clean, cofactors, diagnostics, identity, metal_sites, metals, mmcif, modified, rcsb
from .box import MIN_LIGAND_ATOMS, MIN_LIGAND_MW, WIDE_BOX_EDGE, pick_ligand, place_box

log = logging.getLogger(__name__)


class TargetIdentityError(ValueError):
    """The chain that would be prepared is not the protein that was asked for.

    Raised rather than noted, because every downstream number -- the box, the
    docking score, the MM/GBSA energy, the whole ranking -- is about whatever
    chain this function hands back, and a screen that silently docks a
    different protein produces a confident answer to a question nobody asked.
    """


def fetch_pdb(pdb_id: str, workdir: Path) -> Path:
    """Download (and cache) ``pdb_id``, always returning a PDB-format file.

    Large structures have no PDB-format file at RCSB. The mmCIF is then
    downloaded and converted locally (see :mod:`skinpipe.structure.mmcif`),
    preserving author chain and residue numbering. A structure with more chains
    than the PDB format can address raises
    :class:`~skinpipe.structure.mmcif.StructureTooLargeError` rather than
    silently losing chains.
    """
    workdir = Path(workdir)
    raw_dir = workdir / "pdb"
    path, fmt = rcsb.download_coordinates(pdb_id, raw_dir)
    if fmt == "pdb":
        _refresh_comp_ids(raw_dir / f"{path.stem}.cif", path)
        return path
    converted = raw_dir / f"{path.stem}.pdb"
    if converted.is_file() and converted.stat().st_size > 0:
        _refresh_comp_ids(path, converted)
        return converted
    return mmcif.cif_to_pdb(path, converted)


def _refresh_comp_ids(cif_path: Path, pdb_path: Path) -> None:
    """Re-derive the truncation map for a *cached* PDB when it cannot be trusted.

    `rcsb.download_coordinates` returns a cached PDB before it looks at the
    mmCIF, and the PDB this package converts is written under the same name -- so
    on every run after the first, `cif_to_pdb` does not run and the sidecar it
    would have written is never revisited. A cache built by an earlier converter
    therefore keeps whatever map that converter produced: none at all (the
    truncated code stands, which is the defect the sidecar exists to fix), or a
    pre-guard map that named one component for a prefix two components share.
    Nothing in the file says which, so a stamp mismatch is treated as no
    evidence and the map is re-derived from the mmCIF. That touches only the
    sidecar; the coordinates are left exactly as they are.

    A missing mmCIF is not an error: an entry RCSB serves in PDB format has no
    truncation to record. It does mean a *converted* PDB whose mmCIF has been
    deleted cannot be re-derived, and `box._comp_id_map` is what refuses to
    trust an unverifiable map.
    """
    if not (cif_path.is_file() and cif_path.stat().st_size > 0):
        return
    try:
        if mmcif.sidecar_is_current(cif_path, pdb_path):
            return
        mmcif.write_comp_id_sidecar(cif_path, pdb_path)
        log.info("%s: re-derived the component-id map from %s (the cached PDB "
                 "was converted by an earlier version)",
                 pdb_path.name, cif_path.name)
    except Exception as exc:                      # noqa: BLE001 - advisory only
        log.warning("%s: could not re-derive the component-id map from %s (%s); "
                    "a full five-character component id may be unavailable and "
                    "the redock control may be skipped", pdb_path.name,
                    cif_path.name, exc)


def list_cocrystal_ligands(pdb_path: Path) -> list[dict[str, Any]]:
    """Every non-water HETATM group in the file, largest first.

    Each entry carries ``resname``, ``chain``, ``resnum``, ``n_atoms`` and
    ``centroid``, plus ``icode``, ``is_metal`` and ``drug_like`` so callers can
    filter without re-deriving the blacklist.
    """
    text = Path(pdb_path).read_text()
    groups: dict[tuple[str, str, str, str], list[list[float]]] = {}
    elements: dict[tuple[str, str, str, str], list[str]] = {}
    names: dict[tuple[str, str, str, str], set[str]] = {}
    protein: list[list[float]] = []
    for line in clean.atom_lines(text):
        element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
        if line.startswith("ATOM"):




            if element != "H":
                protein.append([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])])
            continue
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        if clean.is_water(resname):
            continue
        if line[16] not in (" ", "A"):
            continue
        key = (resname, line[21], line[22:26].strip(), line[26])
        groups.setdefault(key, []).append(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
        elements.setdefault(key, []).append(element)
        names.setdefault(key, set()).add(line[12:16].strip().upper())
    prot = np.asarray(protein, dtype=float) if protein else np.zeros((0, 3))

    out: list[dict[str, Any]] = []
    for (resname, chain_id, resnum, icode), coords in groups.items():
        arr = np.asarray(coords, dtype=float)
        metal = clean.is_metal(resname)
        els = elements.get((resname, chain_id, resnum, icode), [])
        heavy = np.asarray([c for c, e in zip(coords, els) if e != "H"],
                           dtype=float) if els else arr
        if heavy.size == 0:
            heavy = arr
        n_heavy = int(heavy.shape[0])
        carbons = sum(1 for e in els if e == "C")
        hetero = sum(1 for e in els if e not in ("C", "H"))


















        polymer_residue = clean._CORE_BACKBONE.issubset(
            names.get((resname, chain_id, resnum, icode), set()))
        out.append(
            {
                "resname": resname,
                "chain": chain_id,
                "resnum": int(resnum) if resnum.lstrip("-").isdigit() else resnum,
                "icode": icode.strip(),
                "n_atoms": int(arr.shape[0]),
                "centroid": tuple(float(v) for v in arr.mean(axis=0)),
                "extent": tuple(float(v) for v in (arr.max(axis=0) - arr.min(axis=0))),
                "is_metal": metal,









                "n_heavy_atoms": n_heavy,
                "n_contacts": _contacts(prot, heavy),
                "contacts_per_atom": (_contacts(prot, heavy) / n_heavy
                                      if n_heavy else 0.0),
                # A membrane lipid or detergent is nearly all carbon. Recorded
                # rather than used to exclude: all-trans retinoic acid is
                # C20H28O2 -- 91% carbon, two heteroatoms -- and it is the
                # ligand for RARG, a node in this registry. Burial decides; this
                # only breaks ties and tells the reader what was chosen.
                "carbon_fraction": (carbons / n_heavy) if n_heavy else 0.0,
                "n_heteroatoms": hetero,
                "lipid_like": bool(n_heavy and carbons / n_heavy >= 0.85
                                   and hetero <= 2),
                "polymer_residue": polymer_residue,
                "drug_like": (
                    not metal
                    and not polymer_residue
                    and resname not in rcsb.HET_BLACKLIST
                    and n_heavy >= MIN_LIGAND_ATOMS
                ),
            }
        )





    out.sort(key=lambda d: (d["drug_like"], d["contacts_per_atom"],
                            not d["lipid_like"], d["n_heavy_atoms"]), reverse=True)
    return out


def _contacts(protein: "np.ndarray", ligand: "np.ndarray",
              cutoff: float = 4.5) -> int:
    """Protein heavy atoms within `cutoff` of any ligand heavy atom.

    Box-prefiltered rather than a full distance matrix: a 100k-atom entry
    against a 30-atom ligand is three million distances per group, paid once
    per group, and the answer only involves the handful of protein atoms in the
    ligand's neighbourhood.
    """
    if protein.size == 0 or ligand.size == 0:
        return 0
    lo = ligand.min(axis=0) - cutoff
    hi = ligand.max(axis=0) + cutoff
    near = protein[np.all((protein >= lo) & (protein <= hi), axis=1)]
    if near.size == 0:
        return 0
    d = np.linalg.norm(near[:, None, :] - ligand[None, :, :], axis=2)





    return int((d <= cutoff).sum())


def select_target_chain(
    raw_text: str,
    ligands: list[dict[str, Any]],
    entities: list[dict[str, Any]],
    code: str,
    chain: Optional[str] = None,
    ligand_code: Optional[str] = None,
    uniprot: Optional[str] = None,
    gene: Optional[str] = None,
) -> tuple[str, Optional[dict[str, Any]], list[str]]:
    """Which chain to dock, checked against the protein that was asked for.

    Returns ``(chain, entity, notes)``.

    Without a wanted identity this is the historical rule: the chain carrying
    the largest drug-like cocrystal ligand, falling back to the first polymer
    chain. That rule is only safe on a single-protein entry. Given a ``uniprot``
    (or a ``gene``), the choice is made *within* the chains RCSB maps to that
    accession, so a bigger ligand sitting on a co-crystallised partner cannot
    pull the preparation onto the wrong protein -- and an explicit ``chain``
    that belongs to a different accession is refused outright.

    When RCSB cannot be asked at all the identity cannot be checked. That is
    recorded as ``NOT VERIFIED`` in the notes rather than assumed, because an
    entity lookup that failed is not evidence that the chain is right.
    """
    notes: list[str] = []
    want_acc = (uniprot or "").strip().upper() or None
    want_gene = (gene or "").strip().upper() or None
    target_chains = _chains_for_identity(entities, want_acc, want_gene)

    if want_acc or want_gene:
        wanted = want_acc or want_gene
        if not entities:
            # "RCSB says this entry has no polymer entity" and "RCSB did not
            # answer" are opposite claims about the same empty list, and only
            # the first is a fact about the structure. The second is fixed by
            # re-running; saying the first sends the reader to doubt an entry
            # that is fine.
            why = _ENTITY_FETCH_ERRORS.get(code)
            if why:
                notes.append(
                    f"target identity NOT VERIFIED **this run**: the polymer-entity "
                    f"lookup for {code} did not answer ({why[:120]}), so nothing "
                    f"confirms the prepared chain is {wanted}. Re-run to check it; "
                    "the entry itself may be perfectly well annotated.")
            else:
                notes.append(
                    f"target identity NOT VERIFIED: RCSB reported no polymer entity for "
                    f"{code}, so nothing confirms the prepared chain is {wanted}")
        elif not target_chains:
            raise TargetIdentityError(
                f"{code}: no chain maps to {wanted}"
                + (f" ({want_gene})" if want_acc and want_gene else "")
                + f"; entry carries {_describe_entities(entities)}. No verified "
                "target mapping was established; annotations may be incomplete.")

    available = clean.chain_ids(raw_text)
    if chain is None and (want_acc or want_gene) and not entities \
            and _ENTITY_FETCH_ERRORS.get(code):



















        raise TargetIdentityError(
            f"{code}: the polymer-entity lookup did not answer "
            f"({_ENTITY_FETCH_ERRORS[code][:120]}), so nothing constrains which "
            f"chain carries {want_acc or want_gene}. Refusing to pick one from "
            f"all {len(available)} chain(s) in the file: that is how a screen "
            "silently changes which protein it docked. Re-run, or pass `chain` "
            "explicitly if you already know which one it is.")

    if chain is None:
        pool = ligands if target_chains is None else [
            lig for lig in ligands if lig["chain"] in target_chains]
        chosen = pick_ligand(pool, None, ligand_code)
        pickable = available if target_chains is None else [
            c for c in available if c in target_chains]
        chain = chosen["chain"] if chosen else (pickable[0] if pickable else None)
        if chain is None and target_chains is not None:
            # RCSB maps the target to chains this coordinate file does not
            # contain. Falling back to whatever chain IS here would prepare a
            # different protein, which is the failure this guard exists for.
            raise TargetIdentityError(
                f"{code}: {want_acc or want_gene} is mapped to chain(s) "
                f"{', '.join(sorted(target_chains))}, none of which is in the "
                f"coordinate file (it has {', '.join(available) or 'none'})")
        if chain is None:
            raise ValueError(f"{code}: no polymer chain found")
        if target_chains is None:
            notes.append(f"chain {chain} selected automatically")
        else:
            notes.append(
                f"chain {chain} selected automatically from the "
                f"{len(target_chains)} chain(s) RCSB maps to "
                f"{want_acc or want_gene}: {', '.join(sorted(target_chains))}")

    entity = _entity_for_chain(entities, chain, strict=True)
    notes.extend(_verify_identity(code, chain, entity, entities, want_acc, want_gene))
    if entity is not None and want_acc:
        # The first reference of a fusion can describe an expression partner.
        # Site annotations must be queried for the requested component.
        entity = dict(entity, uniprot=want_acc)
    if entity is None:
        # The lenient lookup answers `entities[0]`, which its own docstring
        # calls "exactly wrong for deciding whether the chain is the protein
        # that was asked for" -- and this is not only labelling: the entity
        # sets `Receptor.gene` and feeds the tier-2 UniProt binding-site box,
        # so guessing it moves the box onto another polymer's annotations.
        # `_verify_identity` returns early and says nothing when the caller
        # named no accession or gene, so that substitution used to leave no
        # trace anywhere in the notes.
        entity = _entity_for_chain(entities, chain)
        if entity is not None and not (want_acc or want_gene):
            who = entity.get("uniprot") or entity.get("gene") or "an unnamed entity"
            notes.append(
                f"target identity NOT VERIFIED: chain {chain} of {code} is in no "
                f"listed polymer entity, so {who} was assumed for the receptor "
                "label and for the UniProt binding-site box. Pass an accession "
                "or gene to have this checked.")
    return chain, entity, notes


def _chains_for_identity(
    entities: list[dict[str, Any]],
    want_acc: Optional[str],
    want_gene: Optional[str],
) -> Optional[set[str]]:
    """Chains RCSB maps to the wanted protein; ``None`` when nothing was asked.

    An empty set means the lookup answered and the target is not in the entry,
    which is a very different thing from not having asked.
    """
    if not (want_acc or want_gene) or not entities:
        return None
    out: set[str] = set()
    for ent in entities:
        acc = (ent.get("uniprot") or "").strip().upper() or None
        accessions = {str(a).strip().upper() for a in ent.get("uniprots") or []}
        if acc:
            accessions.add(acc)
        sym = (ent.get("gene") or "").strip().upper() or None
        if (want_acc in accessions if want_acc else want_gene and sym == want_gene):
            out.update(ent.get("chains") or [])
    return out


def _verify_identity(
    code: str,
    chain: str,
    entity: Optional[dict[str, Any]],
    entities: list[dict[str, Any]],
    want_acc: Optional[str],
    want_gene: Optional[str],
) -> list[str]:
    """Confirm ``chain`` is the wanted protein, or refuse; always leave a note."""
    if not (want_acc or want_gene):
        return []
    wanted = want_acc or want_gene
    if not entities:
        return []                       # already noted as NOT VERIFIED
    if entity is None:
        return [f"target identity NOT VERIFIED: chain {chain} of {code} is in no "
                f"polymer entity RCSB reported, so nothing confirms it is {wanted}"]

    acc = (entity.get("uniprot") or "").strip().upper() or None
    accessions = {str(a).strip().upper() for a in entity.get("uniprots") or []}
    if acc:
        accessions.add(acc)
    sym = (entity.get("gene") or "").strip().upper() or None
    if (want_acc in accessions if want_acc else want_gene and sym == want_gene):
        return [f"chain {chain} confirmed as {want_acc if want_acc in accessions else acc or '(no accession)'}"
                f"{f' ({sym})' if sym else ''}, the requested {wanted}"]
    if accessions or sym:
        raise TargetIdentityError(
            f"{code}: chain {chain} is {acc or '(no accession)'}"
            f"{f' ({sym})' if sym else ''}, not the requested {wanted}"
            f"{f' ({want_gene})' if want_acc and want_gene else ''}. "
            "Preparing it would dock a different protein")
    return [f"target identity NOT VERIFIED: RCSB gives chain {chain} of {code} no "
            f"accession or gene, so nothing confirms it is {wanted}"]


def _describe_entities(entities: list[dict[str, Any]]) -> str:
    parts = []
    for ent in entities:
        label = ent.get("uniprot") or ent.get("gene") or ent.get("description") or "?"
        parts.append(f"{label}[{','.join(ent.get('chains') or [])}]")
    return ", ".join(parts) or "no polymer entity"


#: Bumped when a change makes an old prepared receptor unfit for reuse.
RECEPTOR_IDENTITY_VERSION: int = 4


def _identity(pdb_id: str, chain: Optional[str], ligand_code: Optional[str],
              ph: float, box_padding: float, keep_waters: bool, *,
              native_model: Optional[dict] = None,
              assembly_model: Optional[dict] = None) -> dict:
    """What a prepared receptor has to match to be reused.

    Everything that changes the file, and nothing that does not. The seed and
    the OpenMM platform are in here because they decide where PDBFixer's
    rebuilt side chains land, and the reason this reuse exists at all is that
    they are not *enough*: seeded and on the Reference platform, 8PQD chain B
    and 7BW1 chain A still prepare differently on an AMD 5950X and an Intel
    13700K -- 51 and 22 rebuilt heavy atoms whose 200-step Langevin settles
    into a different local minimum on each. 9GU4, with 146, happens to agree.
    Which basin wins is decided by numpy, BLAS and libm, not by anything this
    repository controls.

    So a receptor is shipped and reused rather than re-derived. A collaborator
    on another machine scores against the same atoms, and the pose cache -- keyed
    on `receptor_atoms_sha256` -- stays valid instead of re-docking the screen.
    """
    from .clean import PDBFIXER_PLATFORM, PDBFIXER_SEED
    from ..dock.software import file_identity, receptor_preparation_identity

    identity = {
        "version": RECEPTOR_IDENTITY_VERSION,
        "pdb_id": pdb_id, "chain": chain, "ligand_code": ligand_code,
        "ph": float(ph), "box_padding": float(box_padding),
        "keep_waters": bool(keep_waters),
        "pdbfixer_seed": PDBFIXER_SEED,
        "pdbfixer_platform": PDBFIXER_PLATFORM,
        "preparation_tools": receptor_preparation_identity(),
        "preparation_source_sha256": {
            path.name: file_identity(path)
            for path in sorted(Path(__file__).parent.glob("*.py"))
        },
    }
    if native_model is not None:
        identity["native_modified_repair_model"] = native_model
    if assembly_model is not None:
        identity["deposited_assembly_model"] = assembly_model
    return identity


def _new_receptor_directory(workdir: Path, code: str, chain: str, identity: dict) -> Path:
    """Allocate a fresh preparation directory without replacing prior artifacts."""
    import hashlib
    import json
    import tempfile

    parent = Path(workdir) / "receptors"
    parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
    return Path(tempfile.mkdtemp(prefix=f"{code}_{chain}_{digest}_", dir=parent))


def _box_ok(v: object) -> bool:
    """Three finite numbers, which is the whole contract a search box has."""
    if not isinstance(v, (list, tuple)) or len(v) != 3:
        return False
    try:
        return all(math.isfinite(float(x)) for x in v)
    except (TypeError, ValueError):
        return False


def _reuse_receptor(workdir: Path, want: dict) -> "Optional[Receptor]":
    """A previously prepared receptor for exactly this request, or None.

    The stored `receptor.sha256` is re-checked against the file on disk, so a
    truncated or hand-edited pdbqt is re-derived rather than trusted.
    """
    import json as _json

    from ..dock.cachekey import structure_sha256

    if (want.get("version", 0) >= 4 and
            want.get("preparation_tools", {}).get("package_inputs_state") != "observed"):
        log.warning("receptor preparation package identity is unverified; cache reuse refused")
        return None

    for card in sorted(Path(workdir).glob("receptors/*/receptor_identity.json")):
        try:
            stored = _json.loads(card.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if stored.get("identity") != want:
            continue
        if stored.get("preparation_interrupted"):
            continue
        run_id = stored.get("execution_run_id")
        if run_id is not None:
            # A supervisor killed by an outer timeout cannot close its receipt.
            # A late-written identity card is admitted only after its own
            # supervisor recorded a completed, successful bounded attempt.
            if not isinstance(run_id, str) or len(run_id) != 32 or not all(c in "0123456789abcdef" for c in run_id):
                continue
            try:
                deadline = _json.loads((Path(workdir) / "bounded_preparations" / run_id / "deadline.json").read_text())
            except (OSError, ValueError):
                continue
            if deadline.get("state") != "PREPARED_REQUIRING_DOCKING_VALIDATION":
                continue
        schema = stored.get("schema_version")
        if schema is not None and (type(schema) is not int or schema not in (1, 2)):
            log.warning("%s: unsupported identity-card schema %r", card.parent.name, schema)
            continue
        if (stored.get("receptor", {}).get("assembly")
                or stored.get("identity", {}).get("deposited_assembly_model")) and schema != 2:
            log.warning("%s: assembly receptor has no portable artifact manifest", card.parent.name)
            continue







        if stored["identity"].get("chain") is None:
            log.warning("%s: the stored receptor card was written before a "
                        "chain was resolved; re-deriving", card.parent.name)
            continue
        pdbqt = card.parent / str(stored.get("pdbqt") or "")
        clean = card.parent / str(stored.get("clean_pdb") or "")
        if not (pdbqt.is_file() and clean.is_file()):
            continue
        if structure_sha256(pdbqt) != stored.get("receptor_sha256"):
            log.warning("%s: the stored receptor no longer matches its recorded "
                        "hash; re-deriving", card.parent.name)
            continue
        kw = dict(stored.get("receptor") or {})
        if stored.get("schema_version") == 2:
            try:
                _verify_identity_artifacts(card, stored)
            except (OSError, TypeError, ValueError) as exc:
                log.warning("%s: stored receptor artifacts failed verification (%s)",
                            card.parent.name, exc)
                continue
        reference = kw.get("cocrystal_sdf")
        if reference:
            reference = Path(reference)
            if not reference.is_absolute():
                reference = card.parent / reference
            if not reference.is_file():
                log.warning("%s: a recorded reference is missing; refusing stale reuse",
                            card.parent.name)
                continue
            kw["cocrystal_sdf"] = reference
        # A card without a usable box is not a usable receptor. `Receptor`
        # defaults `center` to the origin and `size` to 22 A, so an absent box
        # does not announce itself -- it docks the whole screen in a box that
        # was never chosen, somewhere near (0, 0, 0). An empty list instead
        # reached the unpack in the engine as `ValueError: not enough values`,
        # hundreds of ligands into a run. Both are re-derived here, where the
        # cost is minutes.
        if not (_box_ok(kw.get("center")) and _box_ok(kw.get("size"))):
            log.warning("%s: the stored receptor card has no usable search box "
                        "(center=%r, size=%r); re-deriving",
                        card.parent.name, kw.get("center"), kw.get("size"))
            continue
        kw["path"] = clean
        kw["pdbqt"] = pdbqt
        kw["center"] = tuple(float(x) for x in kw["center"])
        kw["size"] = tuple(float(x) for x in kw["size"])
        try:
            rec = Receptor(**kw)
        except TypeError as exc:                          # an older card
            log.info("stored receptor card is from an older schema (%s); "
                     "re-deriving", exc)
            continue
        rec.notes = list(rec.notes or []) + [
            "receptor reused from this workdir rather than re-derived: "
            "preparation is reproducible on one machine but not across CPUs "
            "for entries with rebuilt side chains"]
        log.info("%s_%s: reusing the prepared receptor already in this workdir",
                 stored["identity"]["pdb_id"], stored["identity"]["chain"])
        return rec
    return None


def _artifact_sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _card_artifact(card: Path, raw: object) -> Path:
    path = Path(str(raw))
    if path.is_absolute():
        raise ValueError("identity artifacts must use relative paths")
    target = (card.parent / path).resolve()
    if not target.is_relative_to(card.parent.resolve()):
        raise ValueError("identity artifact escapes its package")
    return target


def _verify_identity_artifacts(card: Path, stored: dict) -> None:
    artifacts = stored.get("artifacts") or {}
    required = {"pdbqt", "clean_pdb", "metal_sidecar", "cofactor_sidecar"}
    if not required.issubset(artifacts):
        raise ValueError("incomplete receptor artifact manifest")
    if stored.get("identity", {}).get("target_model") is not None:
        if "target_assessment" not in artifacts or not artifacts["target_assessment"].get("present"):
            raise ValueError("missing source-bound target assessment")
        import json
        target = json.loads(_card_artifact(card, artifacts["target_assessment"]["path"]).read_text())
        if target != stored["identity"]["target_model"]:
            raise ValueError("target assessment differs from requested identity")
    fields = stored.get("receptor") or {}
    selection_model = stored.get("identity", {}).get("receptor_selection_model")
    if selection_model is not None:
        import json
        for key in ("receptor_selection", "altloc_selected_input"):
            if not artifacts.get(key, {}).get("present"):
                raise ValueError("missing source-bound receptor selection evidence")
        selection = json.loads(_card_artifact(card, artifacts["receptor_selection"]["path"]).read_text())
        if (selection != fields.get("preparation_evidence", {}).get("receptor_selection")
                or selection.get("altloc", {}).get("policy") != selection_model.get("altloc")
                or selection.get("cofactors", {}).get("policy") != selection_model.get("cofactors")):
            raise ValueError("receptor selection evidence differs from requested identity")
        decisions = stored.get("identity", {}).get("structural_component_decisions")
        if decisions is not None:
            component_selection = selection.get("cofactors", {})
            if (component_selection.get("structural_component_decisions") != decisions
                    or component_selection.get("decision_policy") != cofactors.COMPONENT_DECISION_POLICY
                    or stored["identity"].get("structural_component_decision_policy") != cofactors.COMPONENT_DECISION_POLICY
                    or component_selection.get("manual_model_status") != "EXPLICIT_MANUAL_MODEL_NOT_CHEMICALLY_QUALIFIED"
                    or component_selection.get("native_eligibility") != "NOT_ESTABLISHED"):
                raise ValueError("structural component decisions differ from cache selection receipt")
    if fields.get("assembly") or stored.get("identity", {}).get("deposited_assembly_model"):
        model = stored.get("identity", {}).get("deposited_assembly_model")
        if not isinstance(model, dict) or not model or fields.get("assembly") != model:
            raise ValueError("assembly identity differs from receptor context")
        required_assembly = {"assembly_context", "assembly_source", "assembly_input", "selected_pdb"}
        if not required_assembly.issubset(artifacts):
            raise ValueError("assembly receptor has no complete context/source artifact binding")
        if not all(artifacts[key].get("present") for key in required_assembly):
            raise ValueError("assembly receptor context/source artifact is absent")
        import json
        context = json.loads(_card_artifact(card, artifacts["assembly_context"]["path"]).read_text())
        if context != model or model.get("primary_chain") != stored["identity"].get("chain"):
            raise ValueError("assembly metadata differs from requested identity")
        if (artifacts["assembly_source"]["sha256"] != model.get("source_sha256")
                or artifacts["assembly_input"]["sha256"] != model.get("assembled_pdb_sha256")):
            raise ValueError("assembly source/input differs from resolved context")
    if fields.get("cocrystal_sdf") and "cocrystal_sdf" not in artifacts:
        raise ValueError("reference has no artifact digest")
    for key, entry in artifacts.items():
        path = _card_artifact(card, entry["path"])
        if path.is_file() != entry["present"]:
            raise ValueError(f"artifact presence changed: {key}")
        if entry["present"] and _artifact_sha256(path) != entry["sha256"]:
            raise ValueError(f"artifact digest changed: {key}")
    for key in ("pdbqt", "clean_pdb"):
        if artifacts[key]["path"] != stored[key] or not artifacts[key]["present"]:
            raise ValueError(f"artifact identity changed: {key}")
    if fields.get("cocrystal_sdf"):
        if (artifacts["cocrystal_sdf"]["path"] != fields["cocrystal_sdf"]
                or not artifacts["cocrystal_sdf"]["present"]):
            raise ValueError("reference identity changed")
    decisions = stored.get("identity", {}).get("structural_component_decisions")
    if decisions is not None:
        if selection_model is None:
            raise ValueError("manual component cache lacks receptor selection identity")
        _verify_component_decision_scope(card, stored, selection, decisions)


def _verify_component_decision_scope(card: Path, stored: dict, selection: dict, decisions: list[dict]) -> None:
    """Recheck protected groups before a manual receptor can be reused/admitted.

    Read the saved source and exact resolved box-ligand instance without
    selecting a new ligand, preparing coordinates or changing a model.
    """
    artifacts = stored["artifacts"]
    model = stored.get("identity", {}).get("deposited_assembly_model")
    if not isinstance(model, dict) or not all(artifacts.get(key, {}).get("present") for key in
                                              ("assembly_source", "altloc_selected_input")):
        raise ValueError("manual component cache lacks its deposited source")
    scope = selection.get("component_decision_scope")
    if not isinstance(scope, dict) or set(scope) != {"selected_ligand_instance"}:
        raise ValueError("manual component cache lacks exact protected-group scope")
    excluded = scope["selected_ligand_instance"]
    ligand = stored["receptor"].get("cocrystal_ligand")
    if (excluded is not None and (not isinstance(excluded, list) or len(excluded) != 4
                                 or not all(isinstance(x, str) for x in excluded))):
        raise ValueError("invalid cached box-ligand instance")
    if (excluded is None) != (ligand is None) or excluded is not None and excluded[0] != ligand:
        raise ValueError("cached box-ligand scope differs from receptor")
    raw = _card_artifact(card, artifacts["altloc_selected_input"]["path"]).read_text()
    records = clean.atom_lines(raw)
    if excluded is not None and not any(line.startswith("HETATM") and cofactors._group_key(line) == tuple(excluded)
                                        for line in records):
        raise ValueError("cached box-ligand instance is absent from source")
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    source = _card_artifact(card, artifacts["assembly_source"]["path"])
    data = MMCIF2Dict(str(source))
    full, ambiguous = mmcif.comp_id_translation(mmcif._column(data, "_atom_site.auth_comp_id", "_atom_site.label_comp_id"))
    full.update({key: None for key in ambiguous})
    actual = {}
    cofactors.select_structural_cofactors(
        records, exclude=tuple(excluded) if excluded is not None else None,
        polymer_residues={(name, cid, number[:-1] if number[-1:].isalpha() else number,
                           number[-1] if number[-1:].isalpha() else "")
                          for cid, number, name in modified.declarations(raw)},
        structural_component_decisions=decisions,
        assembly_sha256=model.get("assembled_pdb_sha256"), component_ids=full, evidence=actual)
    if actual != selection["cofactors"]:
        raise ValueError("cached structural component classification differs from source")


def _write_identity(workdir: Path, rec: "Receptor", want: dict) -> None:
    import json as _json
    from dataclasses import asdict, is_dataclass

    from ..dock.cachekey import structure_sha256

    pdbqt = Path(rec.pdbqt)
    card = pdbqt.parent / "receptor_identity.json"
    fields = asdict(rec) if is_dataclass(rec) else dict(vars(rec))
    for k in ("path", "pdbqt"):
        fields.pop(k, None)
    fields["center"] = [float(x) for x in rec.center]
    fields["size"] = [float(x) for x in rec.size]
    from . import cofactors, metal_sites
    paths = {"pdbqt": pdbqt, "clean_pdb": Path(rec.path),
             "metal_sidecar": metal_sites.sidecar_path(Path(rec.path)),
             "cofactor_sidecar": cofactors.sidecar_path(Path(rec.path))}
    if getattr(rec, "assembly", None):
        paths.update(assembly_context=card.parent / "assembly_context.json",
                     assembly_source=card.parent / "assembly_source.cif",
                     assembly_input=card.parent / "assembly_input.pdb",
                     selected_pdb=card.parent / "selected.pdb")
    if want.get("target_model") is not None:
        paths["target_assessment"] = card.parent / "target_assessment.json"
    if want.get("receptor_selection_model") is not None:
        paths["receptor_selection"] = card.parent / "receptor_selection.json"
        paths["altloc_selected_input"] = card.parent / "altloc_selected_input.pdb"
    for index, path in enumerate(sorted(card.parent.glob("*.preparation.*.json"))):
        paths[f"preparation_diagnostic_{index}"] = path
    if rec.cocrystal_sdf is not None:
        paths["cocrystal_sdf"] = Path(rec.cocrystal_sdf)
        fields["cocrystal_sdf"] = paths["cocrystal_sdf"].resolve().relative_to(
            card.parent.resolve()).as_posix()
    artifacts = {}
    for key, path in paths.items():
        relative = path.resolve().relative_to(card.parent.resolve()).as_posix()
        artifacts[key] = {"path": relative, "present": path.is_file(),
                          "sha256": _artifact_sha256(path) if path.is_file() else None}
    payload = {
        "schema_version": 2,
        "identity": want,
        "pdbqt": pdbqt.name,
        "clean_pdb": Path(rec.path).name,
        "receptor_sha256": structure_sha256(pdbqt),
        "receptor": fields,
        "artifacts": artifacts,
    }
    import os
    payload["execution_run_id"] = os.environ.get("SKINPIPE_PREPARATION_RUN_ID")
    try:
        card.write_text(_json.dumps(payload, indent=1, ensure_ascii=False,
                                    default=str), encoding="utf-8")
    except OSError as exc:                                # noqa: BLE001
        log.warning("could not record the receptor identity (%s); the next run "
                    "will re-derive it", exc)


def _in_docking_box(xyz, center, size):
    return all(abs(float(xyz[i]) - float(center[i])) <= float(size[i]) / 2 for i in range(3))


def retained_metal_evidence(lines, center, size):
    """Inventory observed ions without asserting catalytic role or removing any."""
    atoms = [{'serial': line[6:11].strip(), 'chain_id': line[21],
              'residue_id': line[22:27].strip(), 'residue_name': line[17:20].strip(),
              'atom_name': line[12:16].strip(), 'element': line[76:78].strip(),
              'xyz_angstrom': [float(line[a:a + 8]) for a in (30, 38, 46)]}
             for line in lines]
    for atom in atoms:
        atom['inside_docking_box'] = _in_docking_box(atom['xyz_angstrom'], center, size)
    return {'policy': 'retain_observed_catalytic_metal_name_list_v1', 'atoms': atoms,
            'retained_atom_count': len(atoms), 'inside_box_atom_count': sum(a['inside_docking_box'] for a in atoms),
            'biological_role': 'NOT_ESTABLISHED', 'automatic_metal_removal': False,
            'per_instance_metal_removal': 'NOT_SUPPORTED_BY_STRUCTURAL_COMPONENT_DECISION_CONTRACT'}


def removed_glycan_evidence(text, center, size):
    """Report existing glycan removal; do not add/restore an unqualified model."""
    groups = {}
    for line in clean.atom_lines(text):
        if not line.startswith('HETATM') or line[17:20].strip() not in rcsb.GLYCANS:
            continue
        key = (line[21], line[22:27].strip(), line[17:20].strip())
        atom = {'atom_name': line[12:16].strip(),
                'xyz_angstrom': [float(line[a:a + 8]) for a in (30, 38, 46)]}
        atom['inside_docking_box'] = _in_docking_box(atom['xyz_angstrom'], center, size)
        groups.setdefault(key, []).append(atom)
    return {'policy': 'existing_glycan_exclusion_no_pocket_role_validation',
            'groups': [{'chain_id': k[0], 'residue_id': k[1], 'residue_name': k[2],
                        'observed_atom_count': len(atoms),
                        'inside_box_atom_count': sum(a['inside_docking_box'] for a in atoms)}
                       for k, atoms in groups.items()],
            'pocket_chemistry_equivalence': 'NOT_ESTABLISHED'}


@diagnostics.record_attempt
def prepare_receptor(
    pdb_id: str,
    workdir: Path,
    chain: Optional[str] = None,
    ligand_code: Optional[str] = None,
    ph: float = 7.4,
    box_padding: float = 5.0,
    keep_waters: bool = False,
    uniprot: Optional[str] = None,
    gene: Optional[str] = None,
    assembly_id: Optional[str] = None,
    target_role: str = "protein",
    preparation_platform: Optional[str] = None,
    preparation_cpu_threads: int = 1,
    preparation_timeout: Optional[float] = None,
    explicit_pocket: Optional[dict] = None,
    structural_component_decisions: Optional[list[dict]] = None,
) -> Receptor:
    """Fetch, clean, protonate and box a structure for Vina.

    ``uniprot``/``gene`` are the protein the caller means to dock. Pass them
    whenever the entry was found *from* an accession: the chain is then chosen
    among that protein's chains and verified to be it, and a mismatch raises
    :class:`TargetIdentityError`. Leave them out only when the caller picked
    the entry itself and knows it is not the registry's accession.
    """
    # Reject unusable geometry/protonation settings before fetching structures
    # or starting the bounded worker. These values otherwise reach native code.
    ph = clean.validate_ph(ph)
    from .box import validate_box_padding, validate_explicit_pocket
    box_padding = validate_box_padding(box_padding)
    structural_component_decisions = cofactors.normalize_component_decisions(structural_component_decisions)
    if explicit_pocket is not None:
        explicit_pocket = validate_explicit_pocket(explicit_pocket)
        if ligand_code is not None:
            raise ValueError("Choose either a cocrystal ligand or an explicit residue pocket")
    if preparation_timeout is not None:
        from .bounded import run_preparation
        settings = dict(chain=chain, ligand_code=ligand_code, ph=ph, box_padding=box_padding,
                        keep_waters=keep_waters, uniprot=uniprot, gene=gene, assembly_id=assembly_id,
                        target_role=target_role, preparation_platform=preparation_platform,
                        preparation_cpu_threads=preparation_cpu_threads)
        if explicit_pocket is not None:
            settings["explicit_pocket"] = explicit_pocket
        if structural_component_decisions:
            settings["structural_component_decisions"] = structural_component_decisions
        return run_preparation(pdb_id, workdir, settings, preparation_timeout)
    workdir = Path(workdir)
    preparation_model = clean.preparation_model(preparation_platform, preparation_cpu_threads)
    if target_role not in ("protein", "peptide"):
        raise ValueError("target_role must be protein or peptide")
    code = pdb_id.strip().upper()
    raw_pdb = fetch_pdb(code, workdir)
    raw_text = raw_pdb.read_text()
    notes: list[str] = []

    ligands = list_cocrystal_ligands(raw_pdb)
    # Resolved once, before anything selects on it. The PDB resName field holds
    # three characters, so a five-character component arrives truncated and both
    # `pick_ligand` and the template fetch would otherwise be working from a code
    # that can name a different molecule. `None` marks a prefix shared by more
    # than one component, which cannot be resolved from this file at all.
    from .box import _comp_id_map
    _comp_ids = _comp_id_map(code, workdir, pdb_path=raw_pdb)
    for _lig in ligands:
        if _lig["resname"] in _comp_ids:
            full = _comp_ids[_lig["resname"]]
            _lig["component_id"] = full
            _lig["component_ambiguous"] = full is None
    cif_path = rcsb.download_mmcif(code, workdir / "pdb")
    # Resolve full component IDs before the first chain decision, including a
    # cold workdir where the three-character PDB names were just downloaded.
    _refresh_comp_ids(cif_path, raw_pdb)
    _comp_ids = _comp_id_map(code, workdir, pdb_path=raw_pdb)
    for ligand in ligands:
        if ligand["resname"] in _comp_ids:
            ligand["component_id"] = _comp_ids[ligand["resname"]]
            ligand["component_ambiguous"] = ligand["component_id"] is None
    entities = identity.source_entities(cif_path, _safe_entities(code))
    entities, binding_notes = identity.bind_gene(entities, gene, uniprot)
    notes.extend(binding_notes)

    if chain is None and not clean.chain_ids(raw_text):
        raise ValueError(f"{code}: no polymer chain found in {raw_pdb}")
    chain, entity, identity_notes = select_target_chain(
        raw_text, ligands, entities, code, chain=chain, ligand_code=ligand_code,
        uniprot=uniprot, gene=gene)
    notes.extend(identity_notes)
    target_model = identity.assess(entity, chain, uniprot, target_role)
    context = assembly.build_context(cif_path, raw_text, chain, assembly_id=assembly_id)
    unsupported = sorted({row["polymer_type"] for row in context.metadata["instances"]
                          if row["polymer_type"] and row["polymer_type"] != "polypeptide(L)"})
    if unsupported:
        raise assembly.AssemblyError(f"{code}: receptor preparation does not parameterize assembly polymer types {unsupported}")
    raw_text = context.pdb_text
    _refresh_comp_ids(cif_path, raw_pdb)
    _comp_ids = _comp_id_map(code, workdir, pdb_path=raw_pdb)
    partner_chains = sorted({row["output_chain"] for row in context.metadata["instances"]
                            if row["polymer_type"] and row["output_chain"] != chain})
    notes.append(f"deposited assembly {context.metadata['assembly_id']} resolved in target chain {chain}'s "
                 f"original frame; retained partner chains {','.join(partner_chains) or '(none)'}; "
                 f"selection policy {context.metadata['selection_policy']}; physiological relevance and "
                 "native chemical support require separate checks")

























    resolved_text = clean.resolved_altloc_text(raw_text)
    structural_component_decisions = cofactors.bind_component_decisions(
        structural_component_decisions, clean.atom_lines(resolved_text),
        assembly_sha256=context.metadata.get("assembled_pdb_sha256"), component_ids=_comp_ids)
    native_records = modified.groups(resolved_text)
    native_text = "\n".join(modified.declaration_lines(resolved_text)
                            + [line for lines in native_records.values() for line in lines])
    native_model = modified.parameterization_descriptor(native_text) if native_text else None
    want = _identity(code, chain, ligand_code, ph, box_padding, keep_waters,
                     native_model=native_model, assembly_model=context.metadata)
    want["target_model"] = target_model
    want["preparation_model"] = preparation_model
    if explicit_pocket is not None:
        want["explicit_pocket"] = explicit_pocket
    want["pdbfixer_platform"] = preparation_model["platform"]
    want["receptor_selection_model"] = {"altloc": clean.ALTLOC_POLICY,
                                         "cofactors": cofactors.COFACTOR_POLICY}
    if structural_component_decisions:
        want["structural_component_decisions"] = structural_component_decisions
        want["structural_component_decision_policy"] = cofactors.COMPONENT_DECISION_POLICY
    reused = _reuse_receptor(workdir, want)
    if reused is not None:
        return reused

    # Gene symbol where RCSB has one; otherwise the entity description, which is
    # what non-human targets (fungal tyrosinase, bacterial homologues) carry.
    # This is what the chain turned out to be; `gene` above is what was asked
    # for, and `select_target_chain` has already refused any disagreement.
    receptor_gene = (entity or {}).get("gene") or (entity or {}).get("description") or code

    out_dir = _new_receptor_directory(workdir, code, chain, want)
    import json
    import shutil
    context_pdb = out_dir / "assembly_input.pdb"
    context_pdb.write_text(raw_text)
    (out_dir / "assembly_context.json").write_text(json.dumps(context.metadata, indent=2) + "\n")
    shutil.copy2(cif_path, out_dir / "assembly_source.cif")
    (out_dir / "target_assessment.json").write_text(json.dumps(target_model, indent=2) + "\n")
    resolved_pdb = out_dir / "altloc_selected_input.pdb"
    resolved_pdb.write_text(resolved_text)
    ligands = list_cocrystal_ligands(resolved_pdb)
    for ligand in ligands:
        if ligand["resname"] in _comp_ids:
            ligand["component_id"] = _comp_ids[ligand["resname"]]
            ligand["component_ambiguous"] = ligand["component_id"] is None

    selection_evidence: dict = {}
    protein_lines, metal_lines = clean.select_receptor_lines(
        raw_text, None, keep_waters=keep_waters, keep_metals=True,
        selection_evidence=selection_evidence,
    )
    if not protein_lines:
        raise ValueError(f"{code}: chain {chain} has no protein atoms")
    protein_lines, trim_notes = clean.trim_ragged_termini(protein_lines)
    notes.extend(trim_notes)
    if not protein_lines:
        raise ValueError(f"{code}: chain {chain} has no fully resolved residues")
    # Resolve the coordinate-aware box ligand once, before writing the rigid
    # receptor. Cross-chain reranking or fallback must control cofactor removal
    # too; a second coordinate-free pick could leave the boxed ligand in it.
    center, size, source, box_notes, cocrystal = place_box(
        resolved_text, ligands, chain, ligand_code, box_padding, entity, out_dir, code,
        pdb_path=resolved_pdb,
        retained_chains={line[21] for line in protein_lines},
        retained_serials={int(line[6:11]) for line in protein_lines},
        explicit_pocket=explicit_pocket,
    )
    notes.extend(box_notes)
    metal_evidence = retained_metal_evidence(metal_lines, center, size)
    selection_evidence['metals'] = metal_evidence
    if metal_lines:
        descriptions = ', '.join(f"{a['residue_name']} {a['chain_id']}:{a['residue_id']} "
                                 f"({'inside' if a['inside_docking_box'] else 'outside'} box)" for a in metal_evidence['atoms'])
        notes.append(f"retained observed metal-name-list ions: {descriptions}; catalytic role is not established; "
                     "purification/crystallisation ions may remain. Per-instance metal removal is not supported by "
                     "the structural-component decision contract; no metal was removed automatically")
    glycan_evidence = removed_glycan_evidence(resolved_text, center, size)
    selection_evidence['removed_glycans'] = glycan_evidence
    if glycan_evidence['groups']:
        descriptions = ', '.join(f"{g['residue_name']} {g['chain_id']}:{g['residue_id']} "
                                 f"({g['inside_box_atom_count']} atom(s) in box)" for g in glycan_evidence['groups'])
        notes.append(f"excluded deposited glycan groups under the existing preparation policy: {descriptions}; "
                     "removal was not conditioned on pocket distance and pocket chemistry equivalence is not established")
    notes.append("ligand tautomer enumeration and conformer coverage are not assessed by receptor preparation; "
                 "consult the ligand preparation receipt before interpreting a single rigid receptor model")
    cofactor_evidence: dict = {}
    cofactor_exclude = ((cocrystal["resname"], cocrystal["chain"], str(cocrystal["resnum"]),
                         str(cocrystal.get("icode") or "").strip()) if cocrystal else None)
    if structural_component_decisions:
        selection_evidence["component_decision_scope"] = {
            "selected_ligand_instance": list(cofactor_exclude) if cofactor_exclude is not None else None}
    try:
        cofactor_lines = cofactors.select_structural_cofactors(
            clean.atom_lines(resolved_text), None,
            exclude=cofactor_exclude,
            polymer_residues={(name, cid, number[:-1] if number[-1:].isalpha() else number,
                               number[-1] if number[-1:].isalpha() else "")
                              for cid, number, name in modified.declarations(resolved_text)},
            evidence=cofactor_evidence,
            structural_component_decisions=structural_component_decisions,
            assembly_sha256=context.metadata.get("assembled_pdb_sha256"), component_ids=_comp_ids,
        )
    finally:
        selection_evidence["cofactors"] = cofactor_evidence
        (out_dir / "receptor_selection.json").write_text(json.dumps(selection_evidence, indent=2) + "\n")
    selected = clean.write_pdb(
        clean.seqres_lines(raw_text, None) + modified.declaration_lines(raw_text)
        + protein_lines + modified.retained_conect(raw_text, protein_lines), out_dir / "selected.pdb"
    )

    clean_pdb = out_dir / f"{code}_{chain}_clean.pdb"
    if preparation_platform is None:
        notes.extend(clean.run_pdbfixer(selected, clean_pdb, ph=ph, keep_waters=keep_waters))
    else:
        notes.extend(clean.run_pdbfixer(selected, clean_pdb, ph=ph, keep_waters=keep_waters,
                                      platform=preparation_platform, cpu_threads=preparation_cpu_threads))
    preparation_evidence = {"schema_version": 1, "unobserved_gaps": [],
                            "discarded_rebuilt_residue_count": 0,
                            "observed_heavy_atoms_discarded": False,
                            "observed_alternate_heavy_records_discarded": selection_evidence["altloc"]["discarded_heavy_atom_record_count"],
                            "receptor_selection": selection_evidence}
    if structural_component_decisions:
        preparation_evidence.update(manual_model_status="EXPLICIT_MANUAL_MODEL_NOT_CHEMICALLY_QUALIFIED",
                                    native_eligibility="NOT_ESTABLISHED")
        notes.append("source-bound per-instance structural-component decisions applied; explicit manual model "
                     "is not chemical/partial-structure qualification, native eligibility or biological efficacy")
    import hashlib
    journals = []
    input_sha = hashlib.sha256(selected.read_bytes()).hexdigest()
    output_sha = hashlib.sha256(clean_pdb.read_bytes()).hexdigest()
    for receipt in sorted(out_dir.glob(clean_pdb.name + ".preparation.*.json")):
        try:
            raw = receipt.read_bytes()
            value = json.loads(raw)
        except (OSError, ValueError) as exc:
            log.warning("unreadable prior preparation journal %s: %s", receipt, exc)
            continue
        if isinstance(value, dict) and value.get("state") == "COMPLETED" and value.get("input_sha256") == input_sha and value.get("output_sha256") == output_sha:
            journals.append((value.get("started_at", 0), raw, value))
    if journals:
        _, journal_raw, value = max(journals, key=lambda row: row[0])
        preparation_evidence["journal_sha256"] = hashlib.sha256(journal_raw).hexdigest()
        preparation_evidence["unobserved_gaps"] = [gap for stage in value.get("stages", [])
                                                   for gap in stage.get("gaps", [])]
        reconstruction_stages = [stage for stage in value.get('stages', [])
                                 if stage.get('name') == 'record_reconstructed_heavy_atoms']
        if reconstruction_stages:
            rebuilt = [dict(atom, inside_docking_box=_in_docking_box(atom['xyz_angstrom'], center, size))
                       for stage in reconstruction_stages for atom in stage.get('atoms', [])]
            preparation_evidence.update(reconstruction_evidence_status='RECORDED', reconstructed_heavy_atoms=rebuilt,
                                        reconstructed_inside_box_count=sum(atom['inside_docking_box'] for atom in rebuilt))
            if rebuilt:
                notes.append(f"{preparation_evidence['reconstructed_inside_box_count']} reconstructed heavy atom(s) lie "
                             "inside the docking box; their positions/rotamers are model assumptions, not observed pocket geometry")
        else:
            preparation_evidence['reconstruction_evidence_status'] = 'NOT_RECORDED'
    # PDBFixer picked histidine tautomers without seeing the metals, because they
    # were held out of it.  Put any proton it parked on a coordinating nitrogen
    # back where it belongs before meeko types the receptor.
    notes.extend(metals.repair_metal_sites(clean_pdb, metal_lines))
    if not keep_waters:
        notes.append("waters removed (rigid receptor cannot displace them)")

    pdbqt = out_dir / f"{code}_{chain}_receptor.pdbqt"
    notes.extend(clean.write_receptor_pdbqt(clean_pdb, metal_lines, pdbqt,
                                            cofactor_lines, observed_pdb=selected,
                                            preparation_evidence=preparation_evidence))

    # `clean_pdb` is protein-only: the metals are in the PDBQT, which only Vina
    # reads. MM/GBSA is handed `Receptor.path`, i.e. this file, so record the
    # metals beside it -- otherwise a dicopper tyrosinase and a protein that
    # never had a metal are the same file to the next stage. See
    # skinpipe.structure.metal_sites for the 2Y9X measurement behind this.
    # Same treatment for the structural cofactors, and for the same reason:
    # `write_receptor_pdbqt` puts them in the Vina receptor, `clean_pdb` has
    # none of them, and `Receptor.path` -- what MM/GBSA is handed -- is
    # `clean_pdb`. Unconditional on purpose: an empty record deletes a stale
    # one, and a stale one would make the guard refuse a clean receptor.
    cofactor_sidecar = cofactors.write_sidecar(clean_pdb, cofactor_lines)
    if cofactor_sidecar is not None:
        groups = sorted({f"{record_item[17:20].strip()} {record_item[21]}:{record_item[22:27].strip()}"
                         for record_item in cofactor_lines
                         if record_item.startswith(("ATOM", "HETATM"))})
        notes.append(
            f"{len(groups)} structural cofactor group(s) ({', '.join(groups)}) "
            f"are in the Vina receptor and absent from {clean_pdb.name}; "
            f"recorded in {cofactor_sidecar.name}, and MM/GBSA refuses to "
            "rescore a pose whose pocket wall one of them forms (see "
            "skinpipe.mmgbsa.cofactor_guard)")

    sidecar = metal_sites.write_sidecar(clean_pdb, metal_lines)
    if sidecar is not None:
        notes.append(
            f"{len(metal_lines)} retained metal-name-list atom(s) are absent from "
            f"{clean_pdb.name} and recorded in {sidecar.name}; MM/GBSA refuses to "
            "score a pose that binds one (see skinpipe.mmgbsa.metal_guard)")

    if max(size) > WIDE_BOX_EDGE:
        notes.append(
            f"box edge {max(size):.1f} A is wide - raise exhaustiveness or narrow the box"
        )

    meta = _safe_entry_meta(code)
    receptor = Receptor(
        gene=receptor_gene,
        pdb_id=code,
        path=clean_pdb,
        pdbqt=pdbqt,
        center=center,
        size=size,
        chain=chain,
        cocrystal_ligand=(cocrystal["resname"] if cocrystal else None),
        cocrystal_sdf=(cocrystal["sdf_path"] if cocrystal else None),
        cocrystal_completeness=(cocrystal.get("completeness") if cocrystal else None),
        cocrystal_is_fragment=bool(cocrystal.get("is_fragment") if cocrystal else False),
        cocrystal_atoms=((cocrystal.get("matched_heavy_atoms"),
                          cocrystal.get("ccd_heavy_atoms"))
                         if cocrystal and cocrystal.get("ccd_heavy_atoms") else None),
        resolution=meta.get("resolution"),
        method=meta.get("method"),
        pocket_source=source,
        notes=notes,
        assembly=context.metadata,
        preparation_evidence=preparation_evidence,
    )
    _write_identity(workdir, receptor, want)
    return receptor


def _entity_for_chain(
    entities: list[dict[str, Any]], chain: str, strict: bool = False
) -> Optional[dict[str, Any]]:
    """The entity ``chain`` belongs to.

    ``strict=True`` returns ``None`` when no entity actually lists the chain.
    The lenient fallback to ``entities[0]`` is fine for labelling a
    single-entity file and is exactly wrong for deciding whether the chain is
    the protein that was asked for, so identity checks pass ``strict``.
    """
    for ent in entities:
        if chain in ent.get("chains", []):
            return ent
    if strict:
        return None
    return entities[0] if entities else None


def prepare_target_receptor(pdb_id: str, workdir: Path, *, uniprot: str,
                            gene: Optional[str] = None, max_attempts: int = 3,
                            initial_failure: Optional[Exception] = None,
                            **settings) -> Receptor:
    """Try bounded alternative structures for an automatically selected target.

    Every attempt uses the same exact accession and preparation settings.
    This is not appropriate for an explicit chain/assembly/reference ligand:
    those requests must not be substituted. A changed entry is recorded in a
    separate selection receipt; the failed original and its atoms stay intact.
    """
    import hashlib
    import json
    import time
    from uuid import uuid4

    if type(max_attempts) is not int or not 1 <= max_attempts <= 8:
        raise ValueError("max_attempts must be an integer from 1 to 8")
    if not uniprot or any(settings.get(k) is not None for k in ("chain", "assembly_id", "ligand_code")):
        raise ValueError("alternative selection requires an exact accession and no explicit chain, assembly or ligand")
    pdb_id, uniprot = pdb_id.strip().upper(), uniprot.strip().upper()
    path = Path(workdir) / "structure_selections" / (uuid4().hex + ".json")
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt = {"policy": "bounded_exact_accession_alternatives_v1", "requested_accession": uniprot,
               "initial_pdb_id": pdb_id, "max_attempts": max_attempts, "attempts": [],
               "state": "RUNNING", "docking_validated": False, "started_at": time.time(),
               "preparation_settings": settings}

    def write():
        temp = path.with_suffix(".tmp")
        temp.write_text(json.dumps(receipt, indent=2) + "\n")
        temp.replace(path)

    candidates = [pdb_id]
    last = None
    def retryable(exc):
        from .bounded import PreparationTimeoutError
        return isinstance(exc, (TargetIdentityError, identity.TargetStructureError, rcsb.RCSBError,
                                PreparationTimeoutError)) \
            or any(term in str(exc) for term in ("native modified", "observed heavy atom", "connectivity", "assembly"))

    if initial_failure is not None:
        last = initial_failure
        receipt["attempts"].append({"pdb_id": pdb_id, "state": "FAILED",
            "error_type": type(initial_failure).__name__, "error": str(initial_failure),
            "source": "caller_preserved_failure"})
        if not retryable(initial_failure):
            receipt["state"] = "FAILED"
            write()
            raise initial_failure
        write()
    for number in range(1 if initial_failure is not None else 0, max_attempts):
        if number == len(candidates):
            try:
                ids = rcsb.entries_for_uniprot(uniprot)
                meta = rcsb.entry_metadata(ids)
            except Exception as exc:
                receipt.update(state="ALTERNATIVE_LOOKUP_FAILED", error_type=type(exc).__name__, error=str(exc))
                write()
                raise
            receipt["structure_query"] = getattr(ids, "source_scope", None)
            alternatives = []
            for code in ids:
                if code in candidates:
                    continue
                info = meta.get(code) or {}
                components = [lig for lig in info.get("ligands") or []
                    if lig.get("het_code") not in rcsb.HET_BLACKLIST
                    and type(lig.get("n_heavy_atoms")) is int
                    and lig["n_heavy_atoms"] >= MIN_LIGAND_ATOMS
                    and (lig.get("mw") or 0.) >= MIN_LIGAND_MW]
                method = str(info.get("method") or "").lower()
                key = (not bool(components), not ("diffraction" in method or "x-ray" in method),
                       info.get("resolution") if info.get("resolution") is not None else 99., code)
                alternatives.append((key, dict(info, pdb_id=code,
                    has_druglike_ligand=bool(components), druglike_ligands=components)))
            alternatives.sort(key=lambda item: item[0])
            receipt["alternative_candidates"] = [info["pdb_id"] for _key, info in alternatives]
            preferred = _first_with_ligand_on_target(alternatives, uniprot)
            receipt["preparation_skipped_structures"] = [
                {"pdb_id": info["pdb_id"], "reason": info["target_structure_rejection"],
                 "evidence_source": "rcsb_entity_annotations"}
                for _key, info in alternatives if info.get("target_structure_rejection")]
            admissible = [info["pdb_id"] for _key, info in alternatives
                          if not info.get("target_structure_rejection")]
            if preferred and preferred["pdb_id"] in admissible:
                admissible.remove(preferred["pdb_id"])
                admissible.insert(0, preferred["pdb_id"])
            candidates.extend(admissible[:max_attempts - len(candidates)])
        if number >= len(candidates):
            break
        code = candidates[number]
        attempt = {"pdb_id": code, "state": "RUNNING", "started_at": time.time()}
        receipt["attempts"].append(attempt)
        write()
        try:
            receptor = prepare_receptor(code, workdir, uniprot=uniprot, gene=gene, **settings)
        except Exception as exc:
            attempt.update(state="FAILED", error_type=type(exc).__name__, error=str(exc), finished_at=time.time())
            last = exc
            write()
            if not retryable(exc):
                receipt["state"] = "FAILED"
                write()
                raise
            continue
        attempt.update(state="PREPARED", finished_at=time.time())
        receipt.update(state="PREPARED_REQUIRING_DOCKING_VALIDATION", selected_pdb_id=code,
                       receptor_pdbqt_sha256=hashlib.sha256(Path(receptor.pdbqt).read_bytes()).hexdigest())
        write()
        if code != pdb_id:
            receptor.notes.append(f"alternative receptor {code} selected after failed preparation of {pdb_id}; "
                                  f"exact accession {uniprot} preserved; {len(receipt['attempts'])} attempts "
                                  f"recorded in {path.name}; docking and functional validation remain required")
        return receptor
    receipt["state"] = "ALL_ATTEMPTED_STRUCTURES_FAILED"
    write()
    raise RuntimeError(f"no prepared receptor after {len(receipt['attempts'])} accession-bound attempt(s); "
                       f"selection receipt {path.name}; last failure: {last}") from last


#: Why `_safe_entities` came back empty, keyed by PDB code. An outage and a
#: genuinely entity-less entry both arrive as `[]`, and the note written about
#: that emptiness says "RCSB reported no polymer entity" -- an authoritative
#: negative, when what happened was that nobody answered. Observed on 1HFC,
#: which has one: the same entry verified as P03956 on a run minutes earlier
#: and came back NOT VERIFIED on the next when the search endpoint timed out.
#: Module-level because the emptiness travels several frames before anything
#: writes prose about it, and threading a reason through every signature to
#: reach one sentence is not worth the churn.
_ENTITY_FETCH_ERRORS: dict[str, str] = {}


def _safe_entities(code: str) -> list[dict[str, Any]]:
    _ENTITY_FETCH_ERRORS.pop(code, None)
    try:
        return rcsb.polymer_entities(code)
    except rcsb.RCSBError as exc:
        log.warning("polymer entity lookup for %s failed: %s", code, exc)
        _ENTITY_FETCH_ERRORS[code] = str(exc)
        return []


def _safe_entry_meta(code: str) -> dict[str, Any]:
    try:
        return rcsb.entry_metadata([code]).get(code, {})
    except rcsb.RCSBError as exc:
        log.warning("entry metadata for %s failed: %s", code, exc)
        return {}


def best_structure_for_uniprot(uniprot: str, prefer_holo: bool = True) -> Optional[dict[str, Any]]:
    """Pick the most dockable PDB entry mapped to ``uniprot``.

    Ranking, in order: a candidate cocrystal component, experimental method
    (diffraction, then other experimental, then predicted), then resolution.
    Components must also satisfy the same six-heavy-atom floor as coordinate
    pocket selection; molecular weight alone misclassified ER3, a single ion.
    ``prefer_holo=False`` drops the first criterion. The returned metadata
    records query coverage and does not establish pharmacological relevance.
    """
    ids = rcsb.entries_for_uniprot(uniprot)
    if not ids:
        log.warning("no PDB entries mapped to UniProt %s", uniprot)
        return None
    meta = rcsb.entry_metadata(ids)

    ranked: list[tuple[tuple, dict[str, Any]]] = []
    for pdb_id in ids:
        info = dict(meta.get(pdb_id) or {"pdb_id": pdb_id, "method": "unknown"})
        druglike = [
            lig
            for lig in info.get("ligands") or []
            if lig["het_code"] not in rcsb.HET_BLACKLIST
            and (lig.get("mw") or 0.0) >= MIN_LIGAND_MW
            and type(lig.get("n_heavy_atoms")) is int
            and lig["n_heavy_atoms"] >= MIN_LIGAND_ATOMS
        ]
        druglike.sort(key=lambda d: d.get("mw") or 0.0, reverse=True)
        method = (info.get("method") or "unknown").lower()
        is_xray = "x-ray" in method or "diffraction" in method
        is_predicted = "predicted" in method or pdb_id.startswith(("AF_", "MA_"))
        info["has_druglike_ligand"] = bool(druglike)
        info["ligand_metadata_unknown"] = [lig["het_code"]
            for lig in info.get("ligands") or []
            if lig["het_code"] not in rcsb.HET_BLACKLIST
            and (lig.get("mw") or 0.0) >= MIN_LIGAND_MW
            and (type(lig.get("n_heavy_atoms")) is not int or lig["n_heavy_atoms"] < 0)]
        info["best_ligand"] = druglike[0]["het_code"] if druglike else None
        info["druglike_ligands"] = druglike
        key = (
            0 if (druglike and prefer_holo) else 1,
            0 if is_xray else (2 if is_predicted else 1),
            info.get("resolution") if info.get("resolution") is not None else 99.0,
        )
        ranked.append((key, info))

    ranked.sort(key=lambda kv: kv[0])












    if prefer_holo:
        best = _first_with_ligand_on_target(ranked, uniprot)
    else:
        best = ranked[0][1]
    if not best:
        return None
    best["n_candidates"] = len(ids)
    best["structure_rejections"] = [{"pdb_id": info["pdb_id"], "reason": info["target_structure_rejection"]}
                                    for _key, info in ranked if info.get("target_structure_rejection")]
    best["functional_structure_admissibility"] = "NOT_ESTABLISHED"
    best["structure_query"] = getattr(ids, "source_scope", None)
    on_target = best.get("ligand_on_target_chain")
    best["rationale"] = (
        f"{'holo' if best['has_druglike_ligand'] else 'component eligibility not established' if best['ligand_metadata_unknown'] else 'apo'}, "
        f"{best.get('method')}, "
        f"resolution {best.get('resolution')}"
        # Three-valued on purpose: verified on target, verified *off* target,
        # and not established because the lookup did not answer. Collapsing the
        # last two is the conflation this module refuses everywhere else --
        # "RCSB says no" and "RCSB did not say" are different facts and only
        # the first is about the structure.
        + (f", cocrystal ligand on the {uniprot} chain" if on_target is True
           else f", NO cocrystal ligand on a {uniprot} chain"
           if on_target is False
           else f", whether the cocrystal ligand is on a {uniprot} chain could "
                "not be established"
           if "ligand_on_target_chain" in best else "")
    )
    log.info("best structure for %s: %s (%s)", uniprot, best["pdb_id"], best["rationale"])
    return best


#: How many entries to verify the ligand chain for before giving up. Each costs
#: one `polymer_entities` call; the first candidate usually passes.
_LIGAND_CHAIN_PROBES: int = 8


def _first_with_ligand_on_target(ranked: list, uniprot: str) -> dict[str, Any]:
    """The best-ranked entry whose drug-like ligand is on a `uniprot` chain.

    Falls back to the top of the ranking when none of the entries examined
    qualifies, recording `ligand_on_target_chain=False` so the caller and the
    report can say the box will not be placed from this protein's own ligand.
    """
    want = (uniprot or "").strip().upper()
    fallback = ranked[0][1] if ranked else {}
    probed = 0
    for _key, info in ranked:
        if not info.get("has_druglike_ligand"):
            continue
        if probed >= _LIGAND_CHAIN_PROBES:
            break
        probed += 1
        entities = _safe_entities(info["pdb_id"])
        if not entities:
            # The lookup did not answer. Claiming the ligand is on the wrong
            # chain would be as wrong as claiming it is on the right one.
            continue
        chains = _chains_for_identity(entities, want, None) or set()
        rejected_chains = set()
        rejected_reasons = []
        for entity in entities:
            for cid in set(entity.get("chains") or []) & chains:
                try:
                    identity.assess(entity, cid, want)
                except identity.TargetStructureError as exc:
                    rejected_chains.add(cid)
                    rejected_reasons.append(str(exc))
        if chains and rejected_chains == chains:
            info["target_structure_rejection"] = "; ".join(rejected_reasons)
            continue
        chains -= rejected_chains
        for lig in info.get("druglike_ligands") or []:
            if set(lig.get("chains") or []) & chains:
                info["ligand_on_target_chain"] = True
                info["best_ligand"] = lig["het_code"]
                return info
        info["ligand_on_target_chain"] = False
        log.info("%s: its drug-like ligand(s) %s sit on chain(s) outside %s "
                 "(%s maps to %s); looking further down the ranking",
                 info["pdb_id"],
                 [record_item["het_code"] for record_item in info.get("druglike_ligands") or []],
                 want, want, sorted(chains) or "no chain in this entry")
    fallback = next((info for _key, info in ranked if not info.get("target_structure_rejection")), {})
    if fallback:
        fallback.setdefault("ligand_on_target_chain", None)
    return fallback
