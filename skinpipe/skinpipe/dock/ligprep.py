"""SMILES to Vina-ready ligand PDBQT.

Split out of :mod:`skinpipe.dock.vina` purely for file size; the public entry
point ``prepare_ligand_pdbqt`` is re-exported there and should be imported from
``skinpipe.dock.vina``.

ETKDGv3 is used with a fixed seed so a screen is reproducible, followed by
MMFF94s minimisation (UFF where MMFF cannot type an atom). Failures are logged
and reported as ``None`` rather than raised: a library screen that dies on its
first unparsable natural product is useless.
"""
from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import tempfile
import time
import uuid
from pathlib import Path
from typing import Any, Optional

from rdkit import Chem
from rdkit.Chem import AllChem

from . import protonate

log = logging.getLogger(__name__)


def prepare_ligand_pdbqt(
    smiles: str,
    ligand_id: str,
    workdir: Path,
    n_conformers: int = 1,
    ph: Optional[float] = 7.4,
    *,
    execution: Optional[dict[str, Any]] = None,
    verified_protonation_trace=None,
    expected_trace_spec_sha256: Optional[str] = None,
    expected_admission_sha256: Optional[str] = None,
) -> Optional[Path]:
    """SMILES -> 3D -> Meeko PDBQT. ``None`` (never an exception) on failure.

    ``ph`` titrates the ligand through Open Babel's rule-based transform before
    embedding, because docking the wrong ionisation state is its own failure
    mode: it is invisible in the pose count, unfixable by more sampling, and it
    silently deletes salt bridges. See :mod:`skinpipe.dock.protonate`. Pass
    ``ph=None`` to dock the SMILES exactly as written.

    The file and its generation receipt bind input SMILES, exact pH,
    conformer budget, executing preparation software and output bytes.
    The file on disk is keyed on those inputs, not on
    ``ligand_id`` alone. Keying on the id alone means that correcting a wrong
    structure in the input library re-docks the ligand against the **old**
    molecule: ``dock.engine`` correctly invalidates its own result cache on a
    SMILES change, then this function hands it back the stale PDBQT, and the
    output table reports the new compound's name and id against the old
    compound's poses. The same collision happens between a titrated and an
    untitrated run of the same library.

    An optional original-byte ``VerifiedProtonationTrace`` may reuse a closed
    pH substep. Its specification and effective-graph admission require two
    separately reviewed external SHA256 pins. This path fails closed on
    unassigned stereochemistry before embedding; it never chooses an isomer
    or certifies a physiological protomer. Source and effective graphs both
    reject unassigned persistent atom/bond stereo before embedding.
    """
    workdir = Path(workdir)
    if execution is not None:
        execution.clear()
    out_dir = workdir / "ligands"
    if not valid_ligand_id(ligand_id):
        log.warning("unsafe ligand id %r; preparation rejected", ligand_id)
        return None
    source_mol = Chem.MolFromSmiles(smiles)
    if source_mol is None:
        log.warning("%s: RDKit could not parse input SMILES %r", ligand_id, smiles)
        return None
    if source_mol.GetStereoGroups():
        log.warning("%s: enhanced stereochemistry requires an explicit isomer; preparation held", ligand_id)
        return None
    from ..platform.identity import stereochemistry_assessment
    source_stereo = stereochemistry_assessment(source_mol)
    if source_stereo["unassigned_stereo"] or source_stereo["unassigned_bond_stereo"]:
        log.warning("%s: input graph has unspecified atom/bond stereochemistry; preparation held", ligand_id)
        return None
    request = preparation_request(smiles, ph, n_conformers)
    verified = None
    if verified_protonation_trace is not None:
        try:
            from . import verified_ph
            if type(verified_protonation_trace) is not verified_ph.VerifiedProtonationTrace:
                raise ValueError('VerifiedProtonationTrace original bytes required')
            verified = verified_protonation_trace.verify(
                smiles=smiles, ligand_id=ligand_id, ph=ph,
                expected_spec_sha256=expected_trace_spec_sha256,
                expected_admission_sha256=expected_admission_sha256,
                current_software=request['software_identity'], current_binary_path=protonate.obabel_path(),
                current_stereo_consumer=verified_ph.stereo_consumer_identity())
            effective = Chem.MolFromSmiles(verified['effective_smiles'])
            if effective is None or effective.GetNumHeavyAtoms() == 0:
                raise ValueError('invalid verified effective graph')
            Chem.SanitizeMol(effective)
            verified_ph.require_effective_stereo(effective)
            # Native strict dispatch can take time or load another backend.
            # Recheck before a warm cache return as well as before publication.
            from . import software
            checked = verified_protonation_trace.verify(
                smiles=smiles, ligand_id=ligand_id, ph=ph,
                expected_spec_sha256=expected_trace_spec_sha256,
                expected_admission_sha256=expected_admission_sha256,
                current_software=software.preparation_identity(), current_binary_path=protonate.obabel_path(),
                current_stereo_consumer=verified_ph.stereo_consumer_identity())
            if checked != verified:
                raise ValueError('verified producer/consumer changed during strict stereo dispatch')
            request['verified_protonation_trace'] = verified
        except (ImportError, ValueError, RuntimeError, KeyError, TypeError, AttributeError) as exc:
            log.warning('%s: verified pH/admission rejected before embedding: %s: %s',
                        ligand_id, type(exc).__name__, exc)
            return None
    elif expected_trace_spec_sha256 is not None or expected_admission_sha256 is not None:
        log.warning('%s: verified pH expected pins supplied without original-byte trace', ligand_id)
        return None
    out_pdbqt = out_dir / f"{ligand_id}-{_request_digest(request)}.pdbqt"
    if read_preparation_receipt(out_pdbqt, request, ligand_id, for_reuse=True) is not None:
        return out_pdbqt
    if out_pdbqt.exists() or _receipt_path(out_pdbqt).exists():
        kind = 'incomplete' if not _pdbqt_usable(out_pdbqt) else 'unbound or changed'
        log.warning("%s: the cached ligand PDBQT is %s; preserving and re-preparing it", ligand_id, kind)
        _preserve_rejected(out_pdbqt)

    input_smiles = smiles
    protonation = {'requested_ph': None, 'input_smiles': smiles, 'output_smiles': smiles,
                   'state': 'not_requested', 'applied': False, 'note': None}
    if verified is not None:
        smiles = verified['effective_smiles']
        protonation = verified['protonation']
    elif ph is not None:
        smiles, note = protonate.protonate_smiles(smiles, ph, execution=protonation)
        if note:
            report = log.info if protonation.get('applied') else log.warning
            report("%s: %s", ligand_id, note)
        if protonation.get('applied') is not True:
            log.warning("%s: requested pH transform was not applied; preparation held. "
                        "Use ph=None only for an explicit as-written calculation.", ligand_id)
            # A failed transform is an attempt, never a usable ligand. Retain
            # its exact request and diagnostic without publishing a PDBQT or
            # overwriting a previous attempt when a later explicit retry works.
            held_dir = out_dir / 'held_preparations'
            held_dir.mkdir(parents=True, exist_ok=True)
            attempt_id = uuid.uuid4().hex
            hold_path = held_dir / f'{ligand_id}-{_request_digest(request)}-{attempt_id}.json'
            _atomic_write(hold_path,
                          json.dumps({'schema_version': 1, 'state': 'held',
                                      'reason': 'requested_ph_protonation_not_applied',
                                      'ligand_id': ligand_id, 'request': request,
                                      'request_sha256': _request_digest(request),
                                      'protonation': protonation, 'reusable': False,
                                      'pdbqt_file': None, 'attempt_id': attempt_id,
                                      'generated_at': time.time()}, indent=2, allow_nan=False) + '\n')
            if execution is not None:
                execution.update(state='held', reason='requested_ph_protonation_not_applied',
                                 protonation=protonation, hold_receipt=str(hold_path))
            return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.warning("%s: RDKit could not parse SMILES %r", ligand_id, smiles)
        return None
    try:
        Chem.SanitizeMol(mol)
    except (ValueError, Chem.AtomValenceException, Chem.KekulizeException) as exc:
        log.warning("%s: sanitisation failed: %s", ligand_id, exc)
        return None
    if mol.GetNumHeavyAtoms() == 0:
        # An empty SMILES parses to a valid, empty molecule; ETKDG then raises.
        log.warning("%s: SMILES %r has no heavy atoms", ligand_id, smiles)
        return None
    if mol.GetStereoGroups():
        log.warning("%s: enhanced stereochemistry requires an explicit isomer; preparation held", ligand_id)
        return None

    effective_stereo = stereochemistry_assessment(mol)
    undefined_stereo = effective_stereo["unassigned_stereo"]
    undefined_bonds = effective_stereo["unassigned_bond_stereo"]
    if undefined_stereo or undefined_bonds:
        log.warning("%s: effective graph has unspecified atom/bond stereochemistry; preparation held", ligand_id)
        return None
    effective_smiles = Chem.MolToSmiles(mol)

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useSmallRingTorsions = True
    n_conf = max(1, int(n_conformers))
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not conf_ids:
        params.useRandomCoords = True
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not conf_ids:
        log.warning("%s: ETKDGv3 could not embed a conformer (%s)", ligand_id, smiles)
        return None

    optimisation: dict[str, Any] = {}
    best_conf = _optimise_and_pick(mol, conf_ids, ligand_id, execution=optimisation)
    if best_conf is None:
        return None

    single = _mol_with_single_conf(mol, best_conf)

    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        setups = MoleculePreparation().prepare(single)
        if not setups:
            log.warning("%s: meeko produced no molecule setup", ligand_id)
            return None
        written = PDBQTWriterLegacy.write_string(setups[0])
    except (ImportError, ValueError, RuntimeError, KeyError) as exc:
        log.warning("%s: meeko ligand preparation failed: %s", ligand_id, exc)
        return None

    pdbqt_text = written[0] if isinstance(written, tuple) else written
    if not pdbqt_text or not pdbqt_text.strip():
        log.warning("%s: meeko wrote an empty PDBQT", ligand_id)
        return None

    if verified is not None:
        # Recheck live producer tables/modules and the independent consumer
        # before publishing either output file. Original trace bytes stay read-only.
        from . import software
        try:
            after = verified_protonation_trace.verify(
                smiles=input_smiles, ligand_id=ligand_id, ph=ph,
                expected_spec_sha256=expected_trace_spec_sha256,
                expected_admission_sha256=expected_admission_sha256,
                current_software=software.preparation_identity(), current_binary_path=protonate.obabel_path(),
                current_stereo_consumer=verified_ph.stereo_consumer_identity())
            if after != verified:
                raise ValueError('verified producer/consumer changed during preparation')
        except (ImportError, ValueError, RuntimeError, KeyError, TypeError, AttributeError) as exc:
            log.warning('%s: verified pH/admission changed before publication: %s: %s',
                        ligand_id, type(exc).__name__, exc)
            return None

    out_dir.mkdir(parents=True, exist_ok=True)
    # Each write has a unique name, including concurrent calls in one process.
    # The two-file publication is not a transaction; readers require the
    # matching full digest, so an interrupted or racing pair cannot be reused.
    _atomic_write(out_pdbqt, pdbqt_text)
    from . import software
    stable = ((True if verified is not None else software.preparation_identity() == request['software_identity'])
              and _stereo_policy_sha256() == request['stereochemistry_policy_sha256'])
    receipt = {'schema_version': 1, 'state': 'completed', 'ligand_id': ligand_id,
               'request': request, 'request_sha256': _request_digest(request),
               'input_smiles': input_smiles, 'effective_smiles': effective_smiles,
               'protonation': protonation, 'undefined_stereocentres': undefined_stereo,
               'undefined_bond_stereochemistry': undefined_bonds,
               'source_stereochemistry': source_stereo, 'effective_stereochemistry': effective_stereo,
               'requested_ph_qualification': ('AS_WRITTEN' if ph is None else
                   'RULE_TRANSFORM_APPLIED' if protonation.get('applied') is True else 'HELD_UNAPPLIED_REQUESTED_PH'),
               'reuse_scope': 'generation byte identity only; requested pH qualification checked separately',
               'model_scope': {
                   'tautomer_enumeration': 'NOT_PERFORMED',
                   'conformer_ensemble': 'one selected minimum from the recorded embedding budget',
                   'native_default_conformer_budget': 1,
                   'physiological_protomer_validation': 'NOT_ESTABLISHED'},
               'embedding': {'method': 'ETKDGv3', 'seed': 42,
                             'use_small_ring_torsions': True,
                             'used_random_coords_retry': bool(params.useRandomCoords),
                             'requested_conformers': n_conf, 'embedded_conformers': len(conf_ids),
                             'selected_conformer': int(best_conf)},
               'optimisation': optimisation, 'generated_at': time.time(),
               'generation_id': uuid.uuid4().hex,
               'pdbqt_file': out_pdbqt.name,
               'pdbqt_sha256': hashlib.sha256(pdbqt_text.encode()).hexdigest(),
               'software_snapshot_stable': stable,
               'reusable': _reusable_preparation(protonation, request, stable),
               'scope': 'generation/input byte identity; not proof of correct physiological protomer, stereochemistry or binding'}
    if verified is not None:
        receipt['protonation_origin'] = verified
    _atomic_write(_receipt_path(out_pdbqt), json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    log.debug("%s: ligand PDBQT written; protonation %s, requested pH %s",
              ligand_id, protonation['state'], ph)
    return out_pdbqt


def valid_ligand_id(value: Any) -> bool:
    """IDs are filename components, never paths inside a working directory."""
    return (isinstance(value, str) and bool(value) and value not in ('.', '..')
            and not any(c in value for c in '/\\\0') and Path(value).name == value)


def _atomic_write(path: Path, text: str) -> None:
    fd, name = tempfile.mkstemp(dir=path.parent, prefix='.'+path.name+'.', suffix='.part')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            handle.write(text)
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def _receipt_path(path: Path) -> Path:
    return path.with_suffix('.preparation.json')


def preparation_request(smiles: str, ph: Optional[float] = 7.4, n_conformers: int = 1,
                        software_identity: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    from . import software
    from ..platform.identity import STEREO_POLICY
    return {'smiles': smiles, 'ph': float(ph) if ph is not None else None,
            'ph_policy': 'requested-ph-requires-applied-v2',
            'stereochemistry_policy': STEREO_POLICY,
            'stereochemistry_policy_sha256': _stereo_policy_sha256(),
            'n_conformers': max(1, int(n_conformers)),
            'software_identity': software.preparation_identity() if software_identity is None else software_identity}


def _stereo_policy_sha256():
    """The shared checker lives outside the historical dock code inventory."""
    return hashlib.sha256((Path(__file__).resolve().parents[1] / "platform" / "identity.py").read_bytes()).hexdigest()


def _request_digest(request: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(request, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def _reusable_preparation(protonation: dict[str, Any], request: dict[str, Any], stable: Any) -> bool:
    if stable is not True or protonation['state'] in ('tool_error', 'empty_output', 'invalid_output'):
        return False
    if request.get('software_identity', {}).get('package_inputs_state') == 'unverified':
        return False
    if (request.get('ph_policy') == 'requested-ph-requires-applied-v2'
            and request.get('ph') is not None and protonation.get('applied') is not True):
        return False
    if protonation['state'].startswith('applied_'):
        babel = request.get('software_identity', {}).get('openbabel', {})
        return babel.get('available') is True and babel.get('runtime_inputs', {}).get('state') == 'observed'
    return True


def read_preparation_receipt(path: Path, request: Optional[dict[str, Any]] = None,
                             ligand_id: Optional[str] = None, for_reuse: bool = False) -> Optional[dict[str, Any]]:
    try:
        receipt = json.loads(_receipt_path(path).read_text())
        requested = receipt['request']
        protonation = receipt['protonation']
        states = ('not_requested', 'unavailable', 'applied_unchanged', 'applied_changed',
                  'tool_error', 'empty_output', 'invalid_output')
        if (protonation['state'] not in states
                or protonation['input_smiles'] != requested['smiles']
                or receipt['input_smiles'] != requested['smiles']
                or protonation['requested_ph'] != requested['ph']
                or type(protonation['applied']) is not bool
                or protonation['applied'] != protonation['state'].startswith('applied_')
                or type(receipt['reusable']) is not bool
                or not isinstance(receipt['effective_smiles'], str)):
            return None
        origin = requested.get('verified_protonation_trace')
        if origin is not None:
            from .verified_ph import canonical, sha
            if (receipt.get('protonation_origin') != origin
                    or origin['execution_origin'] != 'EXACT_VERIFIED_PRIOR_PH_SUBSTEP'
                    or origin['new_protonation_execution'] is not False
                    or origin['original_trace_bytes_modified'] is not False
                    or origin['protonation'] != protonation
                    or origin['effective_smiles'] != protonation['output_smiles']
                    or receipt['effective_smiles'] != Chem.MolToSmiles(Chem.MolFromSmiles(protonation['output_smiles']))
                    or origin['original_protonation_sha256'] != sha(canonical(protonation))):
                return None
        elif 'protonation_origin' in receipt:
            return None
        qualification = ('AS_WRITTEN' if requested['ph'] is None else
                         'RULE_TRANSFORM_APPLIED' if protonation.get('applied') is True
                         else 'HELD_UNAPPLIED_REQUESTED_PH')
        if ('requested_ph_qualification' in receipt
                and receipt['requested_ph_qualification'] != qualification):
            return None
        # Historical schema-1 receipts remain inspectable without being
        # retroactively certified for an unobserved applied pH transform.
        modern = 'software_snapshot_stable' in receipt
        reusable = (_reusable_preparation(protonation, requested, receipt['software_snapshot_stable'])
                    if modern else protonation['state'] not in ('tool_error', 'empty_output', 'invalid_output'))
        if receipt['reusable'] != reusable or modern and type(receipt['software_snapshot_stable']) is not bool:
            return None
        if (receipt.get('schema_version') != 1 or receipt.get('state') != 'completed'
                or not _pdbqt_usable(path) or receipt.get('pdbqt_file') != path.name
                or receipt.get('request_sha256') != _request_digest(receipt['request'])
                or hashlib.sha256(path.read_bytes()).hexdigest() != receipt.get('pdbqt_sha256')
                or request is not None and receipt['request'] != request
                or ligand_id is not None and receipt.get('ligand_id') != ligand_id
                or for_reuse and receipt.get('reusable') is not True):
            return None
        if for_reuse and not _reusable_preparation(protonation, requested, receipt.get('software_snapshot_stable', True)):
            return None
        # Old receipts remain readable as historical observations, but an
        # unapplied requested pH cannot enter any current docking cache.
        if for_reuse and requested.get('ph') is not None and protonation.get('applied') is not True:
            return None
        if for_reuse:
            from ..platform.identity import stereochemistry_assessment
            if (requested.get('stereochemistry_policy_sha256') is not None
                    and requested['stereochemistry_policy_sha256'] != _stereo_policy_sha256()):
                return None
            for smiles in (receipt['input_smiles'], receipt['effective_smiles']):
                molecule = Chem.MolFromSmiles(smiles)
                if molecule is None or molecule.GetStereoGroups():
                    return None
                stereo = stereochemistry_assessment(molecule)
                if stereo['unassigned_stereo'] or stereo['unassigned_bond_stereo']:
                    return None
        return receipt
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        return None


def _preserve_rejected(path: Path) -> None:
    directory = path.parent/'rejected_cache'/uuid.uuid4().hex
    directory.mkdir(parents=True)
    entries = {}
    for original in (path, _receipt_path(path)):
        if original.is_file():
            data = original.read_bytes()
            (directory/original.name).write_bytes(data)
            entries[original.name] = hashlib.sha256(data).hexdigest()
    (directory/'manifest.json').write_text(json.dumps(
        {'state': 'REJECTED_PREPARATION_CACHE_PRESERVED', 'files': entries}, indent=2)+'\n')


#: The two records every Vina-ready ligand PDBQT has: a torsion tree root and
#: at least one atom. Anything shorter is a truncated write, not a molecule.
def _pdbqt_usable(path: Path) -> bool:
    """Whether a cached ligand PDBQT is complete enough to dock.

    Deliberately not a parse: this runs once per ligand per run and only has to
    separate "meeko finished" from "the file stops in the middle". A ROOT line
    with no ATOM/HETATM is the shape a killed write leaves.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    has_root = any(ln.startswith("ROOT") for ln in text.splitlines())
    has_atom = any(ln.startswith(("ATOM", "HETATM")) for ln in text.splitlines())
    return has_root and has_atom


def _warn_about_undefined_stereo(mol: Chem.Mol, ligand_id: str) -> list[int]:
    """Compatibility diagnostic; actual preparation holds atom *or bond* gaps."""
    from ..platform.identity import stereochemistry_assessment
    assessment = stereochemistry_assessment(mol)
    unassigned = assessment['unassigned_stereo']
    if unassigned:
        log.warning("%s: %d persistent stereocentre(s) undefined; preparation held",
                    ligand_id, len(unassigned))
    return unassigned


def _structure_key(smiles: str, ph: Optional[float], n_conformers: int = 1) -> str:
    """Full request digest; the receipt separately binds exact output bytes."""
    return _request_digest(preparation_request(smiles, ph, n_conformers))


def _mol_with_single_conf(mol: Chem.Mol, conf_id: int) -> Chem.Mol:
    out = Chem.Mol(mol)
    out.RemoveAllConformers()
    out.AddConformer(Chem.Conformer(mol.GetConformer(conf_id)), assignId=True)
    return out


def _optimise_and_pick(mol: Chem.Mol, conf_ids: list[int], ligand_id: str,
                       execution: Optional[dict[str, Any]] = None) -> Optional[int]:
    """MMFF94s where possible, UFF otherwise. Returns the lowest-energy conf id."""
    results: list[tuple[int, float]] = []
    trace = execution if execution is not None else {}
    trace.update(max_iterations=800, conformers=[])
    def minimise(field, conf_id, method):
        code = field.Minimize(maxIts=800)
        energy = float(field.CalcEnergy())
        trace['conformers'].append({'conformer': int(conf_id), 'forcefield': method,
                                   'minimise_returncode': int(code),
                                   'energy': energy if math.isfinite(energy) else None})
        results.append((conf_id, energy))
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    except (ValueError, RuntimeError):
        props = None
    if props is not None:
        for conf_id in conf_ids:
            field = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
            if field is None:
                continue
            minimise(field, conf_id, 'MMFF94s')
    if not results:
        log.debug("%s: MMFF unavailable, falling back to UFF", ligand_id)
        for conf_id in conf_ids:
            field = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            if field is None:
                continue
            minimise(field, conf_id, 'UFF')
    if not results:
        log.warning("%s: neither MMFF nor UFF could type the molecule", ligand_id)
        trace['selection'] = 'first_embedded_unoptimised'
        return conf_ids[0]
    trace['selection'] = 'lowest_returned_energy'
    return min(results, key=lambda kv: kv[1])[0]
