'Validate the selected structure and retain explicit calculation provenance.'
from __future__ import annotations

import base64
from copy import deepcopy
import hashlib
import json
from pathlib import Path, PurePosixPath
import tempfile
from urllib.parse import unquote, urlparse

from .artifacts import ArtifactIntegrityError
from .campaign import CampaignCoordinator, _sealed_read, _sealed_write, _identifier
from .identity import chemical_identity, candidate_identity, digest
from .planning import plan_target
from .ranking import bind_target_goals, validate_requested_action
from .queue import BudgetExhausted, identifier

INPUT_SCHEMA = "skinpipe-campaign-inputs-v1"
MAX_UPLOAD_BYTES = 8 * 1024 * 1024
MAX_MANIFEST_BYTES = 1024 * 1024


def receptor_interpretation(fields):
    """Carry geometry caveats without confusing identity admission with quality."""
    notes = [note for note in fields.get("notes", []) if isinstance(note, str)]
    evidence = deepcopy(fields.get("preparation_evidence") or {})
    limitations = []
    if evidence.get("discarded_rebuilt_residue_count", 0) or any(
        note.startswith(("RECEPTOR RESIDUES DISCARDED", "REBUILT RESIDUES DISCARDED")) for note in notes
    ):
        limitations.append("discarded_rebuilt_residues")
    if evidence.get("unobserved_gaps") or any("gap preserved" in note.lower()
        or note.startswith("left ") and "unmodelled" in note for note in notes):
        limitations.append("unobserved_residue_gaps")
    if any(note.startswith(("trimmed N-terminal residue", "trimmed C-terminal residue")) for note in notes):
        limitations.append("incomplete_terminal_residues")
    if any(note.startswith(("no SEQRES in the input", "SEQRES declares")) for note in notes):
        limitations.append("receptor_completeness_unverified")
    if evidence.get("observed_heavy_atoms_discarded") is True:
        limitations.append("observed_receptor_atom_loss")
    component_evidence = (evidence.get("receptor_selection") or {}).get("cofactors") or {}
    if (evidence.get("manual_model_status") == "EXPLICIT_MANUAL_MODEL_NOT_CHEMICALLY_QUALIFIED"
            or component_evidence.get("decision_policy") == "source_bound_structural_component_instance_decisions_v1"):
        limitations.append("manual_structural_component_model_unqualified")
    return {"geometry_unverified": bool(limitations), "limitations": limitations,
            "preparation_evidence": evidence, "preparation_notes": notes,
            "scope": "prepared receptor geometry; deposited-structure binding is not established"}


class _OwnerQueue:
    """Use distinct projects as well as queue owner authorization."""
    def __init__(self, queue, owner):
        self.queue, self.owner = queue, owner

    def enqueue(self, kind, request, *, owner, **kwargs):
        if owner != self.owner:
            raise ValueError("Campaign queue owner mismatch")
        return self.queue.enqueue(kind, request, owner=owner,
                                  project="campaign-" + digest(owner)[:32], **kwargs)

    def get(self, job_id, *, owner):
        if owner != self.owner:
            raise ValueError("Campaign queue owner mismatch")
        return self.queue.get(job_id, owner=owner)

    def get_statuses(self, job_ids, *, owner):
        if owner != self.owner:
            raise ValueError("Campaign queue owner mismatch")
        if hasattr(self.queue, "get_statuses"):
            return self.queue.get_statuses(job_ids, owner=owner)
        return {job_id: self.queue.get(job_id, owner=owner)["status"] for job_id in job_ids}

    def cancel_if_queued(self, job_id, *, owner, **kwargs):
        if owner != self.owner:
            raise ValueError("Campaign queue owner mismatch")
        return self.queue.cancel_if_queued(job_id, owner=owner, **kwargs)


class CampaignAPI:
    def __init__(self, root, catalog, artifacts, queue):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.catalog, self.artifacts, self.queue = catalog, artifacts, queue

    def _owner_root(self, owner):
        identifier(owner, "owner")
        path = self.root / digest(owner)
        path.mkdir(parents=True, exist_ok=True)
        if not (path / "owner.json").exists():
            _sealed_write(path / "owner.json", {"owner": owner})
        return path

    def _coordinator(self, owner):
        return CampaignCoordinator(self._owner_root(owner) / "campaigns", self.artifacts,
                                   _OwnerQueue(self.queue, owner))

    @staticmethod
    def _object(raw):
        def reject(value):
            raise ValueError("Nonfinite JSON is not a scientific input")
        value = json.loads(raw, parse_constant=reject)
        if not isinstance(value, dict):
            raise ValueError("Manifest must be a JSON object")
        return value

    def _json_artifact(self, sha256, maximum=MAX_MANIFEST_BYTES):
        path = self.artifacts.verified_path(sha256)
        if path.stat().st_size > maximum:
            raise ValueError("JSON artifact exceeds bounded manifest size")
        return self._object(self.artifacts.get_bytes(sha256))

    def _population(self, snapshot_id, target_id):
        manifest = self.catalog.manifest(snapshot_id)
        if not isinstance(target_id, str) or target_id not in manifest["target_counts"]:
            raise ValueError("Target does not belong to the frozen snapshot")
        rows = self.catalog.selection_pool(snapshot_id, target_id=target_id, limit=10000)
        existing={r['candidate_id'] for r in rows}
        for cid in manifest.get('provenance',{}).get('reference_control_ids',[]):
            if cid not in existing:
                control = self.catalog.get_candidate(snapshot_id,cid,target_id)
                if control is None:
                    raise ValueError("Frozen reference control is absent from the target snapshot: " + cid)
                rows.append(control)
                existing.add(cid)
        return rows, manifest

    def upload_artifact(self, body):
        if set(body) - {"base64", "sha256", "media_type", "name"}:
            raise ValueError("Artifact upload accepts bytes and metadata only")
        value = body.get("base64")
        if not isinstance(value, str) or len(value) > MAX_UPLOAD_BYTES * 4 // 3 + 4:
            raise ValueError("Campaign artifact upload exceeds 8 MiB")
        raw = base64.b64decode(value, validate=True)
        if len(raw) > MAX_UPLOAD_BYTES or hashlib.sha256(raw).hexdigest() != body.get("sha256"):
            raise ValueError("Uploaded bytes do not match the declared SHA256 or size limit")
        return self.artifacts.put_bytes(raw, media_type=body.get("media_type", "application/octet-stream"),
                                        metadata={"name": str(body.get("name", "uploaded input"))[:200]})

    def _prepared(self, candidate, receipt_sha256):
        """Re-run the native byte/provenance validator without preparing a molecule."""
        from ..dock.ligprep import read_preparation_receipt

        receipt = self._json_artifact(receipt_sha256)
        name = receipt.get("pdbqt_file")
        if not isinstance(name, str) or PurePosixPath(name).name != name or "\\" in name:
            raise ValueError("Preparation receipt needs a safe PDBQT basename")
        identity = candidate_identity(candidate, trust_stored=True)
        if identity.get("unassigned_stereo") or identity.get("unassigned_bond_stereo"):
            raise ValueError("undefined_stereochemistry")
        source_identity = chemical_identity(receipt.get("input_smiles"))
        if (identity["exact_structure_id"] is None
                or source_identity["exact_structure_id"] != identity["exact_structure_id"]):
            raise ValueError("Preparation source differs from frozen exact candidate structure")
        if receipt.get("software_snapshot_stable") is not True:
            raise ValueError("Preparation software identity was not stable")
        # Reject this explicit hold before the general reuse validator, which
        # also rejects unapplied-pH receipts. Keep the actionable cause rather
        # than relabelling it as an unspecified receipt/cache failure.
        request, protonation = receipt.get("request"), receipt.get("protonation")
        if (isinstance(request, dict) and request.get("ph") is not None
                and isinstance(protonation, dict) and protonation.get("applied") is not True):
            raise ValueError("requested_ph_protonation_not_applied")
        pdbqt = self.artifacts.get_bytes(receipt["pdbqt_sha256"])
        with tempfile.TemporaryDirectory(prefix="skinpipe-preparation-check-") as directory:
            path = Path(directory) / name
            path.write_bytes(pdbqt)
            path.with_suffix(".preparation.json").write_bytes(self.artifacts.get_bytes(receipt_sha256))
            checked = read_preparation_receipt(path, ligand_id=candidate["candidate_id"], for_reuse=True)
        if checked is None:
            raise ValueError("Native preparation receipt failed byte, chemistry-input or software checks")
        if checked.get("undefined_stereocentres") or checked.get("undefined_bond_stereochemistry"):
            raise ValueError("undefined_stereochemistry")
        prepared_identity = chemical_identity(candidate["smiles"], prepared_smiles=checked["effective_smiles"],
                                               preparation=checked["request"])
        prepared = {"pdbqt_sha256": checked["pdbqt_sha256"], "preparation_receipt_sha256": receipt_sha256,
                "exact_structure_id": identity["exact_structure_id"],
                "prepared_structure_id": prepared_identity["prepared_structure_id"],
                "protonation": deepcopy(checked["protonation"]),
                "undefined_stereocentres": deepcopy(checked.get("undefined_stereocentres", [])),
                "validation_scope": "native preparation input/output binding; physiological protomer and binding not established"}
        # As-written input (ph=None) is distinct from any requested pH, including
        # 7.4. A requested condition requires an actually applied transform.
        requested_ph = checked["request"].get("ph")
        if requested_ph is not None:
            from ..structure.clean import validate_ph
            ph = validate_ph(requested_ph)
            if checked["protonation"].get("applied") is not True:
                raise ValueError("requested_ph_protonation_not_applied")
            prepared["preparation_ph"] = ph
        return prepared

    def _receptor(self, record, snapshot, target_id):
        """Verify a native receptor bundle and observed exact-accession mapping."""
        from ..structure.identity import source_entities, assess
        from ..structure.prepare import _verify_identity_artifacts

        if not isinstance(record, dict) or set(record) - {"identity_card_sha256", "uniprot_record_sha256"}:
            raise ValueError("Receptor context requires native identity card and UniProt record hashes")
        card = self._json_artifact(record["identity_card_sha256"])
        reference = self._json_artifact(record["uniprot_record_sha256"])
        model = card.get("identity", {}).get("target_model") or {}
        accession, chain = model.get("accession"), model.get("chain")
        if not isinstance(accession, str) or not accession or reference.get("primaryAccession") != accession:
            raise ValueError("Receptor accession differs from its exact UniProt source record")
        actual_protein = "uniprot:" + accession
        provenance = snapshot.get("provenance", {})
        expected = provenance.get("actual_protein_id")
        if expected is not None and expected != actual_protein:
            raise ValueError("Receptor differs from snapshot actual protein")
        genes = {entry.get("geneName", {}).get("value") for entry in reference.get("genes", []) if isinstance(entry, dict)}
        target_gene = provenance.get("target_gene", target_id)
        gene_binding = expected is None and target_id not in {actual_protein, accession}
        if gene_binding and target_gene not in genes:
            raise ValueError("Snapshot target is not the exact primary gene/accession in the source record")
        taxon = reference.get("organism", {}).get("taxonId")
        if not isinstance(taxon, int) or isinstance(taxon, bool):
            raise ValueError("UniProt source lacks organism identity")
        # Bare gene names follow the existing SkinPipe human-gene contract.
        # An explicit accession or frozen taxon may name another organism;
        # a same-named homologue must never silently substitute for a human gene.
        expected_taxon = provenance.get("taxon_id", 9606 if gene_binding else None)
        if expected_taxon is not None and expected_taxon != taxon:
            raise ValueError("Receptor organism differs from the frozen target context")
        entries = card.get("artifacts")
        if not isinstance(entries, dict) or len(entries) > 100:
            raise ValueError("Native receptor artifact manifest is absent or too large")
        with tempfile.TemporaryDirectory(prefix="skinpipe-receptor-check-") as directory:
            directory = Path(directory)
            card_path = directory / "receptor_identity.json"
            card_path.write_bytes(self.artifacts.get_bytes(record["identity_card_sha256"]))
            seen = set()
            total_size = 0
            for entry in entries.values():
                name = entry.get("path")
                path = PurePosixPath(name) if isinstance(name, str) else None
                if (path is None or path.is_absolute() or ".." in path.parts or "\\" in name
                        or not name or name == "receptor_identity.json" or name in seen):
                    raise ValueError("Unsafe or duplicated native receptor artifact path")
                seen.add(name)
                if entry.get("present") is True:
                    source = self.artifacts.verified_path(entry["sha256"])
                    total_size += source.stat().st_size
                    if total_size > 64 * 1024 * 1024:
                        raise ValueError("Receptor validation bundle exceeds 64 MiB")
                    output = directory / path
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_bytes(self.artifacts.get_bytes(entry["sha256"]))
                elif entry.get("present") is not False:
                    raise ValueError("Native receptor artifact presence must be explicit")
            _verify_identity_artifacts(card_path, card)
            pocket = card.get("identity", {}).get("explicit_pocket")
            if pocket is not None:
                from ..structure.box import explicit_residue_box
                selected = entries.get("selected_pdb")
                if not selected or selected.get("present") is not True:
                    raise ValueError("Explicit pocket requires the retained deposited source atoms")
                center, size = explicit_residue_box(
                    (directory / selected["path"]).read_text(), pocket, chain,
                    card["identity"]["box_padding"])
                fields = card["receptor"]
                if (list(center) != fields["center"] or list(size) != fields["size"]
                        or fields.get("pocket_source") != "explicit_source_residues"):
                    raise ValueError("Explicit pocket differs from its retained source residue geometry")
            source = entries.get("assembly_source")
            if not source or source.get("present") is not True:
                raise ValueError("Observed deposited source mapping is required for receptor admission")
            entities = source_entities(directory / source["path"], [])
            matching = [entity for entity in entities if any(
                segment.get("chain") == chain and segment.get("accession") == accession
                for segment in entity.get("reference_segments", []))]
            if len(matching) != 1:
                raise ValueError("Receptor chain does not uniquely contain the claimed observed accession")
            reassessed = assess(matching[0], chain, accession, model.get("target_role", "protein"))
            # Native v1 cards predate target_construct_residues. Reassess them
            # with the current policy (including short-construct rejection),
            # then compare every invariant that their schema actually records.
            if model.get("policy") not in {"observed_reference_target_role_v1", "observed_reference_target_role_v2"}:
                raise ValueError("Unrecognized native target assessment policy")
            keys = ["accession", "chain", "reference_segments", "construct_sequence_length", "target_role"]
            if model["policy"] == "observed_reference_target_role_v2" or "target_construct_residues" in model:
                keys.append("target_construct_residues")
            for key in keys:
                if reassessed.get(key) != model.get(key):
                    raise ValueError("Stored target assessment differs from deposited source mapping: " + key)
        fields = card["receptor"]
        receptor = {"receptor_sha256": entries["pdbqt"]["sha256"], "center": fields["center"], "size": fields["size"],
                    "actual_protein_id": actual_protein, "target_identity_valid": True,
                    "identity_card_sha256": record["identity_card_sha256"],
                    "uniprot_record_sha256": record["uniprot_record_sha256"], "taxon_id": taxon,
                    "source_target_assessment_policy": model["policy"], "reassessed_target_model": reassessed,
                    "functional_admissibility": model.get("functional_admissibility", "NOT_ESTABLISHED"),
                    "validation_scope": "frozen source-mapped exact target and artifact integrity; biological binding not established"}
        receptor["interpretation"] = receptor_interpretation(fields)
        from ..structure.clean import validate_ph
        ph = validate_ph(card.get("identity", {}).get("ph", 7.4))
        if ph != 7.4:
            receptor["preparation_ph"] = ph
        # Share finite box and artifact validation with the actual coordinator.
        CampaignCoordinator(self.root / "validation", self.artifacts)._validate_receptor(receptor, actual_protein)
        return receptor

    def _validate_inputs(self, manifest_sha256, *, owner=None):
        manifest = self._json_artifact(manifest_sha256)
        if manifest.get("schema") != INPUT_SCHEMA or set(manifest) - {
                "schema", "snapshot_id", "target_id", "prepared", "receptor", "additional_receptors", "preparation_holds"}:
            raise ValueError("Unsupported campaign input manifest schema or fields")
        snapshot_id, target_id = manifest["snapshot_id"], manifest["target_id"]
        if owner is not None and self.catalog.snapshot_owner(snapshot_id) != owner:
            raise PermissionError("Campaign inputs belong to another project")
        _, snapshot = self._population(snapshot_id, target_id)
        supplied = manifest.get("prepared", {})
        if not isinstance(supplied, dict) or len(supplied) > 2000:
            raise ValueError("Prepared manifest must map at most2000 candidate IDs to receipt hashes")
        prepared, held = {}, []
        ordered = sorted(supplied.items())
        # Full saved payloads can be large. Registration streams the same
        # 32-row blocks as preparation rather than retaining 2,000 payloads.
        candidates = (candidate
                      for offset in range(0, len(ordered), 32)
                      for candidate in self.catalog.get_candidates(
                          snapshot_id, [cid for cid, _ in ordered[offset:offset + 32]], target_id))
        for (cid, receipt_sha256), candidate in zip(ordered, candidates):
            if candidate is None:
                held.append({"candidate_id": cid, "reason": "candidate_not_in_frozen_target"})
                continue
            try:
                prepared[cid] = self._prepared(candidate, receipt_sha256)
            except (ValueError, KeyError, TypeError, OSError) as error:
                held.append({"candidate_id": cid, "reason": str(error)})
        from .native_prepare import PREPARATION_REASONS
        diagnostics = {}
        failures = manifest.get("preparation_holds", [])
        if not isinstance(failures, list) or len(failures) > 2000:
            raise ValueError("Preparation diagnostic ledger must contain at most 2000 records")
        for failure in failures:
            cid = failure.get("candidate_id") if isinstance(failure, dict) else None
            try:
                if (not isinstance(failure, dict) or set(failure) - {"candidate_id", "reason", "error_type", "original_structure_id", "diagnostic_log_sha256"}
                        or failure.get("reason") not in PREPARATION_REASONS or cid in supplied or cid in diagnostics):
                    raise ValueError("invalid preparation diagnostic")
                candidate = self.catalog.get_candidate(snapshot_id, cid, target_id=target_id)
                if candidate is None or failure.get("original_structure_id") != candidate_identity(candidate, trust_stored=True)["original_structure_id"]:
                    raise ValueError("preparation diagnostic structure binding mismatch")
                error_type = failure.get("error_type")
                if not isinstance(error_type, str) or not 1 <= len(error_type) <= 64 or not error_type.replace('_', '').isalnum():
                    raise ValueError("invalid preparation error type")
                log_sha = failure.get("diagnostic_log_sha256")
                if log_sha is not None:
                    self.artifacts.verified_path(log_sha)
                diagnostics[cid] = {**failure, "status": "preparation_held", "scope": "exact source chemistry retained; this preparation attempt did not complete"}
                held.append(deepcopy(failure))
            except (ValueError, KeyError, TypeError, OSError):
                held.append({"candidate_id": cid, "reason": "preparation_diagnostic_quarantined"})
        receptor, receptor_reason = None, "receptor_manifest_not_supplied"
        if manifest.get("receptor") is not None:
            try:
                receptor = self._receptor(manifest["receptor"], snapshot, target_id)
                receptor_reason = None
            except (ValueError, KeyError, TypeError, OSError) as error:
                receptor_reason = str(error)
        if receptor:
            ph = receptor.get("preparation_ph", 7.4)
            for cid, item in list(prepared.items()):
                if item.get("preparation_ph", 7.4) != ph:
                    held.append({"candidate_id": cid, "reason": "ligand_receptor_preparation_ph_mismatch"})
                    del prepared[cid]
        additional, held_contexts = [], []
        extras = manifest.get("additional_receptors", [])
        if not isinstance(extras, list) or len(extras) > 2:
            raise ValueError("At most2 additional receptor contexts may be supplied")
        for index, extra in enumerate(extras):
            try:
                context = self._receptor(extra, snapshot, target_id)
                if receptor is None or context["actual_protein_id"] != receptor["actual_protein_id"]:
                    raise ValueError("Additional receptor requires the same admitted primary target")
                if context.get("preparation_ph", 7.4) != receptor.get("preparation_ph", 7.4):
                    raise ValueError("Different preparation pH requires a separate campaign input manifest")
                additional.append(context)
            except (ValueError, KeyError, TypeError, OSError) as error:
                held_contexts.append({"index": index, "reason": str(error)})
        context = {}
        if receptor:
            context = {"target_id": target_id, "actual_protein_id": receptor["actual_protein_id"],
                       "target_identity_valid": True, "receptor_id": receptor["receptor_sha256"],
                       "pocket_id": digest([receptor["center"], receptor["size"]]),
                       "additional_receptors": [{**item, "receptor_id": item["receptor_sha256"],
                           "pocket_id": digest([item["center"], item["size"]])} for item in additional]}
        return {"manifest_sha256": manifest_sha256, "snapshot_id": snapshot_id, "target_id": target_id,
                "prepared": prepared, "preparation_holds": diagnostics, "receptor": receptor, "context": context,
                "validation": {"status": "ready" if receptor and prepared else "held",
                    "prepared_count": len(prepared), "held_preparations": held, "receptor_reason": receptor_reason,
                    "additional_contexts_held": held_contexts,
                    "receptor_interpretation": receptor.get("interpretation") if receptor else None,
                    "scientific_validation": "NOT_ESTABLISHED"}}

    def register_inputs(self, body, *, owner="researcher"):
        if set(body) != {"manifest_sha256"}:
            raise ValueError("Register an immutable input manifest hash")
        checked = self._validate_inputs(body["manifest_sha256"], owner=owner)
        input_id = digest({"owner": owner, "manifest_sha256": body["manifest_sha256"]})
        directory = self._owner_root(owner) / "inputs"
        directory.mkdir(exist_ok=True)
        _sealed_write(directory / (input_id + ".json"), {"owner": owner, "manifest_sha256": body["manifest_sha256"]})
        return {"input_id": input_id, "snapshot_id": checked["snapshot_id"], "target_id": checked["target_id"],
                "validation": checked["validation"], "target_context": checked["context"]}

    def _inputs(self, input_id, snapshot_id, target_id, owner):
        _identifier(input_id)
        record = _sealed_read(self._owner_root(owner) / "inputs" / (input_id + ".json"))
        if record.get("owner") != owner or digest(record) != input_id:
            raise ArtifactIntegrityError("Input registry owner or content binding changed")
        checked = self._validate_inputs(record["manifest_sha256"], owner=owner)
        if checked["snapshot_id"] != snapshot_id or checked["target_id"] != target_id:
            raise ValueError("Prepared inputs belong to a different frozen snapshot or target")
        return checked

    def target_context(self, input_id, *, snapshot_id, target_id, owner="researcher"):
        checked = self._inputs(input_id, snapshot_id, target_id, owner)
        return {"context": checked["context"], "validation": checked["validation"]}

    def _plan_inputs(self, body, owner):
        if not isinstance(body, dict) or set(body) - {
                "snapshot_id", "target_id", "input_id", "target_context", "compute_budget_seconds", "desired_action",
                "confirmation_policy", "compute_profile"}:
            raise ValueError("Plans accept frozen snapshot/input references, not candidate rows or scores")
        snapshot_id, target_id = body["snapshot_id"], body["target_id"]
        if body.get('compute_profile','balanced-research') not in ('balanced-research','balanced-research-v1','generalized-research','generalized-research-v1'):
            raise ValueError('New API plans require a pinned balanced or generalized compute profile; historical manifests are preserved separately')
        rows, manifest = self._population(snapshot_id, target_id)
        pointer = body.get("target_context", {})
        if not isinstance(pointer, dict) or set(pointer) - {"input_id"}:
            raise ValueError("target_context must reference a registered input_id")
        input_id = body.get("input_id", pointer.get("input_id"))
        if body.get("input_id") and pointer.get("input_id") and body["input_id"] != pointer["input_id"]:
            raise ValueError("Conflicting prepared input references")
        checked = self._inputs(input_id, snapshot_id, target_id, owner) if input_id else {
            "prepared": {}, "receptor": None, "context": {},
            "validation": {"status": "held", "prepared_count": 0, "held_preparations": [],
                           "receptor_reason": "verified_input_manifest_required", "scientific_validation": "NOT_ESTABLISHED"}}
        if self.catalog.snapshot_owner(snapshot_id) != owner:
            raise PermissionError("Campaign snapshot belongs to another project")
        if 'desired_action' in body:
            requested = validate_requested_action(body['desired_action'])
            checked['context'] = {**checked['context'], 'desired_action': requested,
                                  'goal_interpretation': 'user hypothesis; not verified skin biology'}
            rows = bind_target_goals(rows, desired_action=requested)
        if "confirmation_policy" in body and not isinstance(body["confirmation_policy"], dict):
            raise ValueError("confirmation_policy must be an object")
        plan = plan_target(rows, target_id, snapshot_id=snapshot_id,
                           references=manifest.get("provenance", {}).get("references", []),
                           target_context=checked["context"],
                           compute_budget_seconds=body.get("compute_budget_seconds", 28800), trust_stored_identity=True,
                           confirmation_policy=body.get("confirmation_policy"),
                           compute_profile=body.get("compute_profile", "balanced-research"))
        plan.update(input_id=input_id, input_validation=checked["validation"],
                    target_population_count=manifest["target_counts"][target_id],
                    selection_pool_count=len(rows), target_context=checked["context"])
        return plan, rows, checked

    def plan(self, body, *, owner="researcher"):
        return self._plan_inputs(body, owner)[0]

    def _status(self, campaign_id, owner, *, verified_assets=None):
        coordinator = self._coordinator(owner)
        result, plan = coordinator._status_and_plan(campaign_id, verified_assets=verified_assets)
        lifecycle_path = coordinator._directory(campaign_id) / "lifecycle.json"
        lifecycle = _sealed_read(lifecycle_path) if lifecycle_path.exists() else {"enabled": False, "status": "not_started"}
        publication_path = coordinator._directory(campaign_id) / "publication.json"
        publication = _sealed_read(publication_path) if publication_path.exists() else None
        return {**result, "snapshot_id": plan["snapshot_id"], "target_id": plan["target_id"],
                "input_id": plan.get("input_id"), "input_validation": plan.get("input_validation"),
                "target_population_count": plan.get("target_population_count", plan["population_count"]),
                "plan_id": plan["plan_id"], "budget_seconds": plan["budget_seconds"],
                "compute_submitted": result["submitted_chunks"] > 0, "lifecycle": lifecycle,
                "published_snapshot_id": publication.get("snapshot_id") if publication else None}

    def publish(self, campaign_id, *, owner="researcher"):
        from .campaign_results import publish_campaign_results

        coordinator = self._coordinator(owner)
        manifest = publish_campaign_results(coordinator, self.catalog, campaign_id, owner=owner)
        with coordinator._lock(campaign_id) as directory:
            _sealed_write(directory / "publication.json", {"snapshot_id": manifest["snapshot_id"],
                "manifest_sha256": digest(manifest), "scientific_validation": "NOT_ESTABLISHED"})
        return manifest

    def _run_state(self, campaign_id, owner, enabled):
        coordinator = self._coordinator(owner)
        manifest = coordinator._manifest(campaign_id)
        if enabled and (manifest["receptor"] is None or not any(
                value.get("pdbqt_sha256") for value in manifest["prepared"].values())):
            return {"status": "held", "reason": "verified_receptor_and_prepared_candidates_required",
                    "campaign": self._status(campaign_id, owner)}
        with coordinator._lock(campaign_id) as directory:
            _sealed_write(directory / "lifecycle.json", {"enabled": enabled,
                "status": "running" if enabled else "paused",
                "scope": "collect, advance and dispatch frozen stages; automatic deep computation remains disabled"})
        if enabled:
            return self._tick_campaign(campaign_id, owner)
        return {"status": "paused", "running_jobs_unchanged": True, "campaign": self._status(campaign_id, owner)}

    def cancel_run_jobs(self,campaign_id,*,owner='researcher'):
        'Validate the selected structure and retain explicit calculation provenance.'
        errors, ids = [], []
        try:
            self._run_state(campaign_id,owner,False)
        except Exception as error:
            errors.append({'step': 'pause', 'error_type': type(error).__name__})
        try:
            coordinator=self._coordinator(owner)
            with coordinator._lock(campaign_id) as directory:
                state=_sealed_read(directory/'state.json')
                ids=[r['job_id'] for r in state['submissions'].values() if r.get('job_id')]
        except Exception as error:
            errors.append({'step': 'submission_inventory', 'error_type': type(error).__name__})



        try:
            ids.extend(self.queue.campaign_job_ids(campaign_id, owner=owner))
        except Exception as error:
            errors.append({'step': 'queue_inventory', 'error_type': type(error).__name__})
        ids = sorted(set(ids))
        jobs = []
        for job_id in ids:
            try:
                jobs.append(self.queue.cancel(job_id, owner=owner))
            except Exception as error:
                errors.append({'step': 'cancel_job', 'job_id': job_id, 'error_type': type(error).__name__})
        return {'campaign_id':campaign_id,'jobs':jobs,
                'status': 'cancellation_incomplete' if errors else 'cancellation_requested',
                'pause_status': 'unverified' if any(e['step'] == 'pause' for e in errors) else 'paused',
                'errors': errors}

    def resume_run_jobs(self,campaign_id,*,owner='researcher'):
        self._run_state(campaign_id, owner, False)
        coordinator=self._coordinator(owner)
        with coordinator._lock(campaign_id) as directory:
            state=_sealed_read(directory/'state.json')
            ids=[r['job_id'] for r in state['submissions'].values() if r.get('job_id') and not r.get('retired_by_confirmation_policy')]
        failures = []
        resumed = []
        deferred = []
        for jid in ids:
            try:
                job=self.queue.get(jid,owner=owner)
                if job['status'] in {'cancelled','failed','needs_data','budget_deferred','retry_required'}:
                    self.queue.resume(jid,owner=owner)
                    resumed.append(jid)
            except BudgetExhausted as error:
                deferred.append({'job_id': jid, 'reason': 'budget_deferred', 'detail': str(error)})
            except (ValueError, KeyError) as error:
                failures.append({"job_id": jid, "reason": str(error)})
        if failures:
            return {"campaign_id": campaign_id, "status": "resume_incomplete", "resumed_job_ids": resumed,
                    "held_jobs": failures, "deferred_jobs": deferred,
                    "lifecycle": "paused", "population_preserved": True}
        return {**self._run_state(campaign_id,owner,True), 'resumed_job_ids': resumed, 'deferred_jobs': deferred}

    def _research_cancel_requested(self, campaign_id, owner):
        'Validate the selected structure and retain explicit calculation provenance.'
        with self.queue.connection() as conn:
            if not conn.execute("SELECT to_regclass('platform_discovery_runs') AS name").fetchone()['name']:
                return False  # Legacy standalone campaign schema/semantics.
            return conn.execute("""SELECT 1 FROM platform_discovery_runs r
                LEFT JOIN platform_discovery_controls c ON c.run_id=r.id AND c.owner_id=r.owner_id
                WHERE r.campaign_id=%s AND r.owner_id=%s
                AND COALESCE(c.desired_state,r.desired_state)='cancelled' LIMIT 1""",
                (campaign_id, owner)).fetchone() is not None

    def _tick_campaign(self, campaign_id, owner):
        # A separate scheduler lock prevents two daemon processes from making
        # inconsistent lifecycle decisions. Coordinator methods retain their
        # own publication lock and queue submissions retain idempotency keys.
        import fcntl

        if self._research_cancel_requested(campaign_id, owner):
            return {**self.cancel_run_jobs(campaign_id, owner=owner),
                    'status': 'paused_cancel_requested', 'submitted_chunks': 0}
        coordinator = self._coordinator(owner)
        directory = coordinator._directory(campaign_id)
        with (directory / ".scheduler.lock").open("a+b") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            lifecycle_path = directory / "lifecycle.json"
            lifecycle = _sealed_read(lifecycle_path) if lifecycle_path.exists() else {}
            if lifecycle.get("enabled") is not True:
                return {"campaign_id": campaign_id, "status": lifecycle.get("status", "not_started"), "submitted_chunks": 0}
            collected = coordinator.collect(campaign_id)
            advanced = coordinator.advance(campaign_id)
            status = coordinator.status(campaign_id)
            submitted = 0
            for stage in status["stages"]:
                # Pause is a durable boundary checked before each new stage.
                if _sealed_read(lifecycle_path).get("enabled") is not True:
                    break
                if self._research_cancel_requested(campaign_id, owner):
                    return {**self.cancel_run_jobs(campaign_id, owner=owner),
                            'status': 'paused_cancel_requested', 'submitted_chunks': submitted}
                submitted += coordinator.dispatch(campaign_id, stage=stage, owner=owner)["submitted_count"]
                if self._research_cancel_requested(campaign_id, owner):
                    return {**self.cancel_run_jobs(campaign_id, owner=owner),
                            'status': 'paused_cancel_requested', 'submitted_chunks': submitted}
            status = coordinator.status(campaign_id)
            manifest = coordinator._manifest(campaign_id)
            state = _sealed_read(directory / "state.json")
            required_jobs = []
            scope = state.get("confirmation_scope")
            for item in state["submissions"].values():
                if not item.get("job_id") or item.get("retired_by_confirmation_policy"):
                    continue
                request = item["request"]
                if scope and request.get("stage") == "confirmation":
                    active = set(scope["candidate_ids"] + scope["control_ids"])
                    if (request.get("seed") not in scope["seeds"] or not any(
                            row["candidate_id"] in active for row in request["ligands"])):
                        continue
                required_jobs.append(item)
            queued = [coordinator.queue.get(item["job_id"], owner=owner) for item in required_jobs]
            budget_deferred = any(item["status"] == "budget_deferred" for item in queued)
            # Missing unprepared ligands remain explicit holds. They cannot
            # turn a completed subset into a claim of a completed full pool.
            accounted = True
            for stage in status["stages"]:
                _, missing = coordinator._stage_rows(directory, state, stage)
                if any(manifest["prepared"].get(cid, {}).get("pdbqt_sha256") for cid in missing):
                    accounted = False
            finished = "confirmation" in status["stages"] and accounted
            phase = ("budget_deferred" if budget_deferred else
                     "retry_required" if status.get("retry_required") else
                     "finished_with_holds" if finished and manifest["held"] else
                     "finished" if finished else "running")
            settled_budget = budget_deferred and not any(
                item["status"] in {"queued", "running"} for item in queued
            )
            if (finished or settled_budget) and state["outcomes"] and any(state["outcomes"].values()):
                self.publish(campaign_id, owner=owner)
            with coordinator._lock(campaign_id):
                current = _sealed_read(lifecycle_path)
                if current.get("enabled") is True:
                    _sealed_write(lifecycle_path, {**current, "enabled": not (finished or settled_budget),
                        "status": phase, "last_collection": collected, "last_advance": advanced["status"]})
            return {"campaign_id": campaign_id, "status": phase, "submitted_chunks": submitted,
                    "collected": collected, "advanced": advanced, "campaign": self._status(campaign_id, owner)}

    def tick(self, *, owner=None):
        'Validate the selected structure and retain explicit calculation provenance.'
        owners = [owner] if owner is not None else []
        if owner is None:
            for directory in sorted(self.root.iterdir()):
                if directory.is_dir() and (directory / "owner.json").exists():
                    try:
                        record = _sealed_read(directory / "owner.json")
                    except ArtifactIntegrityError:
                        continue
                    if digest(record.get("owner")) != directory.name:
                        continue
                    owners.append(record["owner"])
        results = []
        for name in owners:
            coordinator = self._coordinator(name)
            for directory in sorted(coordinator.root.iterdir()):
                if not directory.is_dir() or len(directory.name) != 64:
                    continue
                try:
                    results.append(self._tick_campaign(directory.name, name))
                except (ValueError, KeyError, OSError) as error:
                    results.append({"campaign_id": directory.name, "status": "needs_attention", "reason": str(error)})
                except Exception as error:


                    import logging
                    logging.getLogger(__name__).exception("Isolated campaign tick failure: %s", directory.name)
                    results.append({"campaign_id": directory.name, "status": "needs_attention",
                                    "reason": "Unexpected campaign error; inspect private operator logs",
                                    "error_type": type(error).__name__})
        return {"campaigns": results, "automatic_deep_count": 0}

    def get(self, url, *, owner="researcher"):
        parts = [unquote(value) for value in urlparse(url).path.strip("/").split("/")]
        if parts == ["v2", "campaigns"]:
            coordinator = self._coordinator(owner)
            result = []
            verified_assets = set()
            for path in sorted(coordinator.root.iterdir()):
                if path.is_dir() and len(path.name) == 64:
                    try:
                        result.append(self._status(path.name, owner, verified_assets=verified_assets))
                    except (ArtifactIntegrityError, ValueError) as error:
                        result.append({"campaign_id": path.name, "status": "quarantined", "reason": str(error)})
            return {"campaigns": result}
        if len(parts) in {3, 4} and parts[:2] == ["v2", "campaigns"] and (len(parts) == 3 or parts[3] == "status"):
            return self._status(parts[2], owner)
        if len(parts) == 3 and parts[:2] == ["v2", "campaign-inputs"]:
            input_id = _identifier(parts[2])
            record = _sealed_read(self._owner_root(owner) / "inputs" / (input_id + ".json"))
            checked = self._validate_inputs(record["manifest_sha256"], owner=owner)
            return {"input_id": input_id, "snapshot_id": checked["snapshot_id"], "target_id": checked["target_id"],
                    **self.target_context(input_id, snapshot_id=checked["snapshot_id"], target_id=checked["target_id"], owner=owner)}
        raise KeyError("Unknown campaign route")

    def post(self, url, body, *, owner="researcher"):
        if not isinstance(body, dict):
            raise ValueError("Campaign request must be an object")
        parts = [unquote(value) for value in urlparse(url).path.strip("/").split("/")]
        if parts == ["v2", "campaign-inputs", "artifacts"]:
            return self.upload_artifact(body)
        if parts == ["v2", "campaign-inputs"]:
            return self.register_inputs(body, owner=owner)
        if parts == ["v2", "campaigns"]:
            plan, rows, checked = self._plan_inputs(body, owner)
            result = self._coordinator(owner).create(plan, rows, checked["prepared"], checked["receptor"],
                                                     trust_stored_identity=True, preparation_holds=checked.get("preparation_holds"))
            return {**self._status(result["campaign_id"], owner), "plan": plan}
        if len(parts) == 4 and parts[:2] == ["v2", "campaigns"]:
            coordinator = self._coordinator(owner)
            _identifier(parts[2])
            if parts[3] == "confirmation-policy":
                if set(body) != {"policy", "reason"}:
                    raise ValueError("Confirmation amendment requires policy and reason")
                import fcntl
                directory = coordinator._directory(parts[2])
                if not directory.is_dir():
                    raise KeyError("campaign not found")
                with (directory / ".scheduler.lock").open("a+b") as handle:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
                    coordinator.amend_confirmation(parts[2], body["policy"], reason=body["reason"], owner=owner)
                    coordinator.collect(parts[2])
                    advanced = coordinator.advance(parts[2])


                    if (directory / "publication.json").exists():
                        self.publish(parts[2], owner=owner)
                return {**advanced, "campaign": self._status(parts[2], owner)}
            if parts[3] == "dispatch":
                if set(body) - {"stage", "chunk_size", "retry_failed"}:
                    raise ValueError("Dispatch cannot replace frozen candidates, settings or scores")
                if "retry_failed" in body and not isinstance(body["retry_failed"], bool):
                    raise ValueError("retry_failed must be boolean")
                return coordinator.dispatch(parts[2], stage=body.get("stage", "broad_docking"),
                    chunk_size=body.get("chunk_size", 32), retry_failed=body.get("retry_failed", False), owner=owner)
            if body:
                raise ValueError("Collect/advance use frozen campaign inputs only")
            if parts[3] == "run":
                return self._run_state(parts[2], owner, True)
            if parts[3] == "pause":
                return self._run_state(parts[2], owner, False)
            if parts[3] == "publish":
                return {"snapshot": self.publish(parts[2], owner=owner), "campaign": self._status(parts[2], owner)}
            if parts[3] == "collect":
                return {**coordinator.collect(parts[2]), "campaign": self._status(parts[2], owner)}
            if parts[3] == "advance":
                coordinator.collect(parts[2])
                return {**coordinator.advance(parts[2]), "campaign": self._status(parts[2], owner)}
        raise KeyError("Unknown campaign route")


def main(argv=None):
    'Validate the selected structure and retain explicit calculation provenance.'
    import argparse
    import signal
    import threading
    from .cli import read_config, service_from_config

    parser = argparse.ArgumentParser(description="Resume explicitly started frozen SkinPipe campaigns")
    parser.add_argument("--config", required=True)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--poll-seconds", type=int, default=30)
    arguments = parser.parse_args(argv)
    if not 1 <= arguments.poll_seconds <= 3600:
        parser.error("poll seconds must be 1..3600")
    service = service_from_config(read_config(arguments.config))
    api = CampaignAPI(service.root / "campaign-api", service.catalog, service.artifacts, service.queue)
    stop = threading.Event()
    for name in (signal.SIGTERM, signal.SIGINT):
        signal.signal(name, lambda *_: stop.set())
    while not stop.is_set():
        try:
            result = api.tick()
        except Exception as error:
            # Operational connection errors must not expose a DSN or kill the
            # restartable dispatcher. Candidate-specific failures stay in receipts.
            result = {"status": "retry_pending", "error_type": type(error).__name__}
        print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
        if arguments.once:
            break
        stop.wait(arguments.poll_seconds)


if __name__ == "__main__":
    main()
