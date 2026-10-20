"""Authenticate a prior pH substep without running its producer again.

External expected digests are trusted caller inputs, obtained by independent
review, not values calculated from the supplied packet. Authentication proves
execution provenance only. Preparation additionally requires a separately
reviewed effective-graph admission and the live strict stereo checks below.
"""
from __future__ import annotations

import copy
import hashlib
import importlib
import inspect
import json
import math
import re
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any


class VerifiedPhTraceRejected(ValueError):
    """Prior-execution authentication/admission did not pass."""


class EffectiveStereoRejected(ValueError):
    """The opt-in effective graph cannot be embedded without an isomer choice."""


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode()


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise VerifiedPhTraceRejected('verified pH trace rejected: ' + message)


def _json(raw: bytes) -> Any:
    _require(type(raw) is bytes and len(raw) <= 16 * 1024 * 1024, 'bounded original JSON bytes required')

    def pairs(items):
        result = {}
        for key, value in items:
            _require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result

    def constant(value):
        raise VerifiedPhTraceRejected('nonfinite JSON number: ' + value)

    return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)


def _pinned(raw: bytes, expected: str) -> Any:
    _require(isinstance(expected, str) and re.fullmatch('[0-9a-f]{64}', expected) is not None,
             'external expected SHA256 required')
    _require(sha(raw) == expected, 'original bytes differ from external expected SHA256')
    return _json(raw)


def ph_producer_identity(snapshot: dict[str, Any]) -> dict[str, Any]:
    """pH producer closure; embedding/Meeko are separate consumer inputs."""
    babel, rdkit = snapshot['openbabel'], snapshot['rdkit_package']
    _require(snapshot['package_inputs_state'] == 'observed', 'unobserved preparation package')
    _require(babel['available'] is True and babel['runtime_inputs']['state'] == 'observed',
             'unobserved Open Babel executable/rule tables')
    count = babel['runtime_inputs']['files_count']
    _require(type(count) is int and count > 0 and len(babel['runtime_inputs']['files']) == count,
             'incomplete Open Babel rule-table inventory')
    _require(rdkit['available'] is True and rdkit['snapshot_state'] == 'observed'
             and not rdkit['missing_manifest_files'] and rdkit['installed_files_hashed'] > 0,
             'incomplete RDKit package/native inventory')
    digests = [babel['sha256'], *snapshot['rdkit_selected_native_sha256'].values(),
               *(snapshot['code_sha256'][name] for name in (
                   'protonate.py', 'openbabel_inputs.py', 'software.py'))]
    _require(all(isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) for value in digests),
             'missing producer module/native/binary digest')
    return {'pH_modules_sha256': {name: snapshot['code_sha256'][name]
                                 for name in ('protonate.py', 'openbabel_inputs.py', 'software.py')},
            'openbabel': copy.deepcopy(babel), 'rdkit_package': copy.deepcopy(rdkit),
            'rdkit_version': snapshot['rdkit_version'],
            'rdkit_native_sha256': copy.deepcopy(snapshot['rdkit_selected_native_sha256']),
            'native_command_contract': ['-ismi', '-osmi', '-p', 'requested_ph'],
            'canonicalizer': 'actual RDKit MolFromSmiles/MolToSmiles',
            'excluded_from_pH_substep_identity': ['ligprep.py embedding consumer', 'Meeko typing consumer']}


def stereo_consumer_identity() -> dict[str, Any]:
    """Observe actual opt-in OpenFF dispatch bytes; never import on default path."""
    from openff.toolkit.topology import Molecule
    from openff.toolkit.utils import GLOBAL_TOOLKIT_REGISTRY
    from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
    from . import software

    capable = [toolkit for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits
               if callable(getattr(toolkit, 'from_smiles', None))]
    # An OpenEye/custom native backend is not covered by the observed RDKit
    # package. Reject it rather than pretend OpenFF wrapper bytes attest it.
    _require(capable and all(type(toolkit) is RDKitToolkitWrapper for toolkit in capable),
             'only the observed RDKit strict from_smiles backend is supported')
    _require(all(getattr(toolkit.from_smiles, '__func__', None) is RDKitToolkitWrapper.from_smiles
                 for toolkit in capable), 'instance from_smiles override is not supported')

    def callable_identity(call):
        function = getattr(call, '__func__', call)
        chain = [function, inspect.unwrap(function)]
        result = []
        for defining in chain:
            module = importlib.import_module(defining.__module__)
            path = inspect.getfile(defining)
            _require(Path(path).resolve() == Path(inspect.getfile(module)).resolve(),
                     'callable defining bytes are outside their observed module')
            digest = software.file_identity(path)
            _require(digest is not None, 'unobserved strict stereo callable')
            result.append({'module': defining.__module__, 'qualname': defining.__qualname__,
                           'path': path, 'sha256': digest})
        return result

    callables = {'Molecule.from_smiles': callable_identity(Molecule.from_smiles),
                 'registry.call': callable_identity(GLOBAL_TOOLKIT_REGISTRY.call),
                 'capable_from_smiles_in_registry_order': [callable_identity(toolkit.from_smiles)
                                                           for toolkit in capable]}
    objects = [Molecule, type(GLOBAL_TOOLKIT_REGISTRY),
               *(type(toolkit) for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits)]
    defining = {}
    for obj in objects:
        module = importlib.import_module(obj.__module__)
        path = inspect.getfile(module)
        digest = software.file_identity(path)
        _require(digest is not None, 'unobserved OpenFF defining source')
        defining[obj.__module__] = {'path': path, 'sha256': digest}
    package = software._package_identity('openff-toolkit', 'openff')
    _require(package['snapshot_state'] == 'observed' and not package['missing_manifest_files'],
             'unobserved OpenFF package namespace')
    helper = software.file_identity(__file__)
    _require(helper is not None, 'unobserved verified trace consumer source')
    from ..platform import identity as stereo_identity
    shared = software.file_identity(stereo_identity.__file__)
    _require(shared is not None, 'unobserved shared stereochemistry policy source')
    return {'package': package, 'defining_sources': defining, 'actual_callables': callables,
            'registered_toolkits': [
                {'class': type(toolkit).__module__ + '.' + type(toolkit).__qualname__,
                 'name': toolkit.toolkit_name, 'version': toolkit.toolkit_version}
                for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits],
            'from_smiles_backend': 'openff.toolkit.utils.rdkit_wrapper.RDKitToolkitWrapper',
            'verified_ph_module_sha256': helper, 'shared_stereo_policy_sha256': shared}


def require_effective_stereo(mol) -> None:
    """Opt-in only: no unspecified/enhanced configuration is chosen by ETKDG."""
    from rdkit import Chem
    def require(condition, message):
        if not condition:
            raise EffectiveStereoRejected(message)

    require(not mol.GetStereoGroups(), 'enhanced stereochemistry requires separate evaluation')
    from ..platform.identity import stereochemistry_assessment
    assessment = stereochemistry_assessment(mol)
    require(not assessment['unassigned_stereo'] and not assessment['unassigned_bond_stereo'],
             'effective graph has unspecified atom/bond stereochemistry')
    from openff.toolkit.topology import Molecule
    # Preserve native dispatch and strict from_smiles semantics. A missing
    # optional dependency or an UndefinedStereochemistryError is a rejection.
    try:
        # RDKit has independently rejected every persistent atom/bond gap.
        # OpenFF may call a reversible N+-H centre undefined; permit that
        # specific graph feature while retaining its actual formal charge.
        Molecule.from_smiles(Chem.MolToSmiles(mol),
                             allow_undefined_stereo=bool(assessment['reversible_protonated_nitrogen']))
    except Exception as error:
        raise EffectiveStereoRejected(f'{type(error).__name__}: {error}') from error


@dataclass(frozen=True)
class VerifiedProtonationTrace:
    'Original UTF-8 JSON bytes, never rewritten or labelled as a new call.'
    trace_raw: bytes
    producer_receipt_raw: bytes
    producer_manifest_raw: bytes
    producer_software_raw: bytes
    source_row_raw: bytes
    spec_raw: bytes
    admission_raw: bytes

    def to_packet(self) -> dict[str, str]:
        """Lossless JSON transport for a bounded native preparation request."""
        return {field.name: getattr(self, field.name).decode('utf-8') for field in fields(self)}

    @classmethod
    def from_packet(cls, packet: dict[str, str]) -> VerifiedProtonationTrace:
        _require(type(packet) is dict and set(packet) == {field.name for field in fields(cls)},
                 'exact original-byte packet fields required')
        _require(all(type(value) is str for value in packet.values()), 'UTF-8 packet strings required')
        return cls(**{name: value.encode('utf-8') for name, value in packet.items()})

    def verify(self, *, smiles: str, ligand_id: str, ph: float | None,
               expected_spec_sha256: str, expected_admission_sha256: str,
               current_software: dict[str, Any], current_binary_path: str | None,
               current_stereo_consumer: dict[str, Any]) -> dict[str, Any]:
        """Authenticate originals and independently pinned admission, fail closed."""
        spec = _pinned(self.spec_raw, expected_spec_sha256)
        admission = _pinned(self.admission_raw, expected_admission_sha256)
        _require(spec['schema'] == 'skinpipe-verified-ph-trace-v1', 'unsupported specification schema')
        originals = {'trace_sha256': self.trace_raw,
                     'producer_receipt_sha256': self.producer_receipt_raw,
                     'producer_manifest_sha256': self.producer_manifest_raw,
                     'producer_software_sha256': self.producer_software_raw,
                     'source_row_sha256': self.source_row_raw}
        for key, raw in originals.items():
            _pinned(raw, spec[key])
        row, receipt = _json(self.trace_raw), _json(self.producer_receipt_raw)
        manifest, previous = _json(self.producer_manifest_raw), _json(self.producer_software_raw)
        identity = ph_producer_identity(previous)
        _require(ph_producer_identity(current_software) == identity, 'fresh pH producer closure changed')
        _require(sha(canonical(identity)) == spec['pH_producer_identity_sha256'], 'pH producer pin mismatch')
        _require(receipt['status'] == 'INPUT_PROTONATION_PREFLIGHT_COMPLETE_REVIEW_REQUIRED'
                 and receipt['prepared_input_count'] is None
                 and receipt['original108_holds'] == 'UNCHANGED', 'producer is not a closed diagnostic receipt')
        _require(receipt['producer_sha256'] == row['producer_sha256'] == sha(self.producer_software_raw)
                 and receipt['diagnostic_rows_manifest_sha256'] == sha(self.producer_manifest_raw),
                 'producer receipt/trace/manifest lineage mismatch')
        population = spec['producer_population']
        _require(type(population) is int and 0 < population <= 100000
                 and receipt['diagnostic_population'] == len(manifest) == population,
                 'producer population mismatch')
        ids = [item['candidate_id'] for item in manifest]
        _require(len(set(ids)) == population
                 and [item['ordered_index'] for item in manifest] == list(range(population)),
                 'producer order/identity is not unique')
        _require(receipt['admitted_order_sha256'] == spec['admitted_order_sha256'] == sha(canonical(ids))
                 and receipt['control_order_sha256'] == spec['control_order_sha256'], 'population/control order changed')
        index = spec['ordered_index']
        _require(type(index) is int and 0 <= index < population, 'invalid source order index')
        selected = manifest[index]
        _require(selected['candidate_id'] == row['candidate_id'] == spec['candidate_id'] == ligand_id
                 and row['ordered_index'] == index
                 and selected['diagnostic_sha256'] == sha(self.trace_raw)
                 and row['source_row_sha256'] == spec['source_row_sha256'], 'candidate/source/order mismatch')
        _require(type(ph) in (int, float) and math.isfinite(ph)
                 and ph == spec['requested_ph'], 'exact finite requested pH required')
        trace = row['protonation']
        _require(trace['requested_ph'] == ph and trace['applied'] is True
                 and type(trace['returncode']) is int and trace['returncode'] == 0
                 and trace['state'] in ('applied_changed', 'applied_unchanged'), 'prior pH execution did not succeed')
        _require(trace['binary_requested'] == spec['binary_requested'] == current_binary_path
                 and previous['openbabel']['sha256'] == spec['binary_sha256'], 'exact owned binary binding changed')
        _require(trace['input_smiles'] == row['source_graph']['smiles'] == spec['input_smiles'] == smiles
                 and trace['output_smiles'] == row['effective_graph']['smiles'] == spec['effective_smiles'],
                 'exact source/effective graph strings changed')
        _require(sha(smiles.encode()) == row['source_graph']['smiles_sha256']
                 and sha(trace['output_smiles'].encode()) == row['effective_graph']['smiles_sha256'],
                 'source/effective graph digest mismatch')
        ph_request = {'candidate_id': ligand_id, 'source_row_sha256': spec['source_row_sha256'],
                      'input_smiles': smiles, 'requested_ph': float(ph)}
        _require(sha(canonical(ph_request)) == spec['pH_request_sha256'], 'pH request binding changed')
        consumer = {'preparation_software': current_software, 'strict_stereo': current_stereo_consumer}
        admission_binding = {'candidate_id': ligand_id, 'source_row_sha256': spec['source_row_sha256'],
                             'input_smiles': smiles, 'effective_smiles': trace['output_smiles'],
                             'requested_ph': float(ph), 'trace_sha256': sha(self.trace_raw),
                             'trace_spec_sha256': expected_spec_sha256,
                             'consumer_identity_sha256': sha(canonical(consumer))}
        _require(admission['schema'] == 'skinpipe-effective-graph-admission-v1'
                 and admission['decision'] == 'APPROVED_FOR_PREPARATION'
                 and admission['source_stereo_constraints'] == 'INDEPENDENTLY_REVIEWED'
                 and all(admission.get(key) == value for key, value in admission_binding.items()),
                 'separate effective-graph admission does not match this consumer')
        return {'execution_origin': 'EXACT_VERIFIED_PRIOR_PH_SUBSTEP',
                'new_protonation_execution': False, 'original_trace_bytes_modified': False,
                **{key: sha(raw) for key, raw in originals.items()},
                'expected_spec_sha256': expected_spec_sha256,
                'expected_admission_sha256': expected_admission_sha256,
                'pH_request_sha256': spec['pH_request_sha256'],
                'pH_producer_identity_sha256': spec['pH_producer_identity_sha256'],
                'consumer_identity_sha256': sha(canonical(consumer)),
                'strict_stereo_consumer': copy.deepcopy(current_stereo_consumer),
                'original_protonation_sha256': sha(canonical(trace)),
                'protonation': copy.deepcopy(trace), 'effective_smiles': trace['output_smiles'],
                'scope': 'verified prior execution and reviewed preparation admission; not physiological protomer proof'}
