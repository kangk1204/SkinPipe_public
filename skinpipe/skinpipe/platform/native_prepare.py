"""Bounded native preparation in an administrator-configured chemistry runtime.

The API runtime need not install Meeko/OpenMM. Requests contain only validated
identifiers/SMILES; executable and owned work directory come from the server.
Cancellation kills only this preparation process group and leaves receipts.
"""

from pathlib import Path
import json
import logging
import os
import subprocess
import sys
import time
import uuid

from .artifacts import atomic_bytes, canonical_json, sha256_file

PREPARATION_REASONS = {
    "unsupported_multi_fragment_chemistry", "ligand_preparation_empty_pdbqt",
    "forcefield_unsupported_chemistry", "invalid_structure", "conformer_generation_failed",
    "ligand_preparation_failed", "native_preparation_runtime_failed",
    # Admission may stop before embedding: an unspecified stereoisomer is not
    # an engine failure, and no automatic stereoisomer may replace its source.
    "undefined_stereochemistry",
    "enhanced_stereochemistry_requires_explicit_isomer",
    "requested_ph_protonation_not_applied",
}


def preparation_failure_code(messages):
    """Safe diagnostic categories; never change the exact input chemistry."""
    text = "\n".join(messages)
    for marker, reason in (
        ("enhanced stereochemistry requires an explicit isomer", "enhanced_stereochemistry_requires_explicit_isomer"),
        ("requested pH transform was not applied", "requested_ph_protonation_not_applied"),
        ("stereocentre(s) undefined in the input SMILES", "undefined_stereochemistry"),
        ("effective graph has undefined stereochemistry", "undefined_stereochemistry"),
        ("effective graph has unspecified atom/bond stereochemistry", "undefined_stereochemistry"),
        ("input graph has unspecified atom/bond stereochemistry", "undefined_stereochemistry"),
        ("UndefinedStereochemistryError", "undefined_stereochemistry"),
        ("fragments. Must have 1", "unsupported_multi_fragment_chemistry"),
        ("meeko wrote an empty PDBQT", "ligand_preparation_empty_pdbqt"),
        ("neither MMFF nor UFF could type", "forcefield_unsupported_chemistry"),
        ("RDKit could not parse", "invalid_structure"),
        ("sanitisation failed", "invalid_structure"),
        ("ETKDGv3 could not embed", "conformer_generation_failed"),
    ):
        if marker in text:
            return reason
    return "ligand_preparation_failed"


class PreparationUnavailable(RuntimeError):
    pass


def run_native(request, workdir, *, python=None, timeout=360, stop=lambda: None):
    stop()
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    attempt = workdir / ("native-" + uuid.uuid4().hex)
    attempt.mkdir()
    input_path, output_path = attempt / "request.json", attempt / "result.json"
    atomic_bytes(input_path, canonical_json({**request, "workdir": str(workdir)}))
    executable = str(Path(python or sys.executable).expanduser().absolute())
    environment = {
        **{key: os.environ[key] for key in (
            "PATH", "HOME", "TMPDIR", "LANG", "LC_ALL", "SYSTEMROOT", "WINDIR"
        ) if key in os.environ},
        "PYTHONPATH": str(Path(__file__).resolve().parents[2]),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    started = time.monotonic()
    with (attempt / "execution.log").open("wb") as log:
        process = subprocess.Popen(
            [
                executable,
                "-m",
                "skinpipe.platform.native_prepare",
                str(input_path),
                str(output_path),
            ],
            env=environment,
            stdout=log,
            stderr=log,
            start_new_session=True,
        )
        try:
            while process.poll() is None:
                stop()
                if time.monotonic() - started >= timeout:
                    raise PreparationUnavailable("native_preparation_timeout")
                time.sleep(0.1)
            stop()
            if (
                process.returncode
                or not output_path.is_file()
                or output_path.stat().st_size > 1024 * 1024
            ):
                raise PreparationUnavailable("native_preparation_runtime_failed")
            result = json.loads(output_path.read_bytes())
            if result.get("status") != "complete":
                raise PreparationUnavailable(
                    result.get("error_type", "native_preparation_failed")
                )
            for row in result.get("ligands", []) + (
                [result] if result.get("card_path") else []
            ):
                for key in ("path", "receipt_path", "card_path"):
                    if row.get(key) and not Path(row[key]).resolve().is_relative_to(
                        workdir
                    ):
                        raise PreparationUnavailable(
                            "native_artifact_outside_owned_workdir"
                        )
            result["diagnostic_log"] = {"path": str(attempt / "execution.log"),
                "sha256": sha256_file(attempt / "execution.log")}
            return result
        finally:
            from .engines import terminate_owned
            terminate_owned(process)


def execute(request):
    workdir = Path(request["workdir"])
    if request["kind"] == "receptor":
        from ..structure.prepare import prepare_receptor

        settings = dict(request["settings"])
        if "structural_component_decisions" in settings:
            from ..structure.cofactors import normalize_component_decisions
            settings["structural_component_decisions"] = normalize_component_decisions(settings["structural_component_decisions"])
        rec = prepare_receptor(
            settings.pop("pdb_id"),
            workdir,
            uniprot=request["accession"],
            preparation_cpu_threads=1,
            **settings,
        )
        card = Path(rec.pdbqt).parent / "receptor_identity.json"
        return {"status": "complete", "card_path": str(card)}
    if request["kind"] != "ligands" or not 1 <= len(request["ligands"]) <= 32:
        raise ValueError("Bounded native ligand batch required")
    from ..structure.clean import validate_ph
    ph = validate_ph(request.get("ph", 7.4))
    from ..dock import ligprep
    from .identity import digest

    class Messages(logging.Handler):
        def __init__(self):
            super().__init__()
            self.messages = []

        def emit(self, record):
            if len(self.messages) < 16:
                self.messages.append(record.getMessage()[:2048])

    rows = []
    for item in request["ligands"]:
        row = {"candidate_id": item["candidate_id"]}
        messages = Messages()
        ligprep.log.addHandler(messages)
        try:
            from .identity import chemical_identity
            source_identity = chemical_identity(item["smiles"])
            if source_identity.get("enhanced_stereo"):
                row.update(reason="enhanced_stereochemistry_requires_explicit_isomer", error_type="EnhancedStereochemistryHeld",
                           original_structure_id="raw:" + digest(item["smiles"]))
                rows.append(row)
                continue
            if source_identity.get("unassigned_stereo") or source_identity.get("unassigned_bond_stereo"):
                row.update(reason="undefined_stereochemistry", error_type="UndefinedStereochemistryHeld",
                           original_structure_id="raw:" + digest(item["smiles"]))
                rows.append(row)
                continue
            optional = {}
            if 'verified_protonation_trace' in item:
                from ..dock.verified_ph import VerifiedProtonationTrace
                optional = {
                    'verified_protonation_trace': VerifiedProtonationTrace.from_packet(item['verified_protonation_trace']),
                    'expected_trace_spec_sha256': item['expected_trace_spec_sha256'],
                    'expected_admission_sha256': item['expected_admission_sha256'],
                }
            elif any(key in item for key in ('expected_trace_spec_sha256', 'expected_admission_sha256')):
                raise ValueError('verified pH pins without original trace')
            path = ligprep.prepare_ligand_pdbqt(
                item["smiles"], item["candidate_id"], workdir, n_conformers=1, ph=ph, **optional
            )
            if path is None:
                raise ValueError("ligand_preparation_failed")
            receipt = ligprep.read_preparation_receipt(path, ligand_id=item["candidate_id"])
            if receipt is None:
                raise ValueError("ligand_preparation_receipt_invalid")
            if ph is not None and receipt["protonation"].get("applied") is not True:
                row.update(reason="requested_ph_protonation_not_applied", error_type="RequestedProtonationHeld",
                           receipt_path=str(path.with_suffix(".preparation.json")),
                           original_structure_id="raw:" + digest(item["smiles"]))
                rows.append(row)
                continue
            row.update(
                path=str(path), receipt_path=str(path.with_suffix(".preparation.json"))
            )
        except Exception as error:
            row.update(reason=preparation_failure_code(messages.messages), error_type=type(error).__name__,
                       original_structure_id="raw:" + digest(item["smiles"]))
            if 'verified_protonation_trace' in item:
                # Keep authentication and stereo diagnostics distinct without
                # relabelling historical/default failures or inventing holds.
                row['preparation_diagnostics'] = messages.messages or [str(error)[:2048]]
        finally:
            ligprep.log.removeHandler(messages)
        rows.append(row)
    return {"status": "complete", "ligands": rows}


def main():
    input_path, output_path = map(Path, sys.argv[1:])
    try:
        result = execute(json.loads(input_path.read_bytes()))
    except Exception as error:
        result = {"status": "held", "error_type": type(error).__name__}
    atomic_bytes(output_path, canonical_json(result))


if __name__ == "__main__":
    main()
