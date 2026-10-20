"""Skin biomarker registry: resolve a marker to something a molecule can bind."""
from __future__ import annotations

import json
import copy
from contextlib import contextmanager
import hashlib
import os
import tempfile
import logging
import time
from pathlib import Path
from typing import Iterable, Optional

import yaml

from .fetch import chembl, target_identity
from .types import Direction, DruggableNode, TargetMode, TargetSpec

log = logging.getLogger(__name__)

_DEFAULT_CONFIG = Path(__file__).resolve().parent.parent / "config" / "skin_targets.yaml"
_NODE_CACHE_SCHEMA = "primary-gene-node-v1"


@contextmanager
def _cache_write_lock(path: Path):
    """Serialize cooperating cache writers before comparing their snapshots."""
    with path.open("a+b") as handle:
        if os.name == "nt":
            import msvcrt
            if path.stat().st_size == 0:
                handle.write(b"\0")
                handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _verified_node(gene: str, entry: dict) -> bool:
    if not isinstance(entry, dict) or entry.get("node_cache_schema") != _NODE_CACHE_SCHEMA:
        return False
    if entry.get("gene") != gene or entry.get("resolved") is not True:
        return False
    identity = entry.get("identity")
    if not target_identity.validate_primary_binding(identity, gene):
        return False
    if any(entry.get(k) != identity.get(k) for k in ("uniprot", "chembl_target")):
        return False
    if any(type(entry.get(k)) is not int or entry[k] < 0 for k in ("n_pdb", "n_ligands_p6")):
        return False
    return identity["chembl_target"] is not None or entry["n_ligands_p6"] == 0


class Registry:
    """Marker table plus a disk cache of resolved node metadata."""

    def __init__(self, config: Optional[Path] = None, cache: Optional[Path] = None,
                 budget_s: Optional[float] = None):













        self._budget_s = budget_s
        self._deadline: Optional[float] = None
        #: Genes the budget stopped us reaching, kept apart from `_failed`'s
        #: API errors. "We ran out of time" and "the API said no" are both
        #: unresolved, but only the first is a statement about *us*, and a
        #: reader deciding whether to re-run needs to know which they have.
        self._skipped_for_budget: list[str] = []
        self.config_path = Path(config or _DEFAULT_CONFIG)
        raw = yaml.safe_load(self.config_path.read_text())
        self._markers: list[dict] = raw["markers"]
        from .paths import user_cache_dir
        config_hash = hashlib.sha256(self.config_path.read_bytes()).hexdigest()[:16]
        self.cache_path = Path(cache or os.environ.get("SKINPIPE_NODE_CACHE") or
                               user_cache_dir() / "registry" / (config_hash + ".json"))
        legacy_path = self.config_path.parent / ".node_cache.json"
        read_path = self.cache_path
        if cache is None and not os.environ.get("SKINPIPE_NODE_CACHE") and not read_path.exists():
            # Read old metadata for migration, but never write inside the package.
            if legacy_path.exists():
                read_path = legacy_path
        self._cache: dict[str, dict] = {}
        #: Genes whose lookup failed in *this* process. In-memory only, so one
        #: outage does not hammer a dead endpoint 31 times, and never written to
        #: disk, so it cannot outlive the outage that produced it.
        self._failed: dict[str, str] = {}
        self._expected_cache_bytes = self.cache_path.read_bytes() if self.cache_path.exists() else None
        self._rejected_cache_bytes = None
        if read_path.exists():
            original = read_path.read_bytes()
            try:
                loaded = json.loads(original)
                if not isinstance(loaded, dict):
                    raise ValueError("node cache must be an object")
            except (ValueError, UnicodeDecodeError):
                self._rejected_cache_bytes = original
                log.warning("node cache at %s is corrupt; retaining bytes before replacement", read_path)
            else:
                self._cache = {g: e for g, e in loaded.items()
                               if _verified_node(g, e)}
                dropped = len(loaded) - len(self._cache)
                if dropped:
                    self._rejected_cache_bytes = original
                    log.warning("rejected %d unverified node cache entries from %s; retaining bytes before replacement",
                                dropped, read_path)

    # ---------------------------------------------------------------- lookup
    def markers(self, category: Optional[str] = None) -> list[str]:
        return [m["marker"] for m in self._markers
                if category is None or m["category"] == category]

    def categories(self) -> list[str]:
        seen: list[str] = []
        for m in self._markers:
            if m["category"] not in seen:
                seen.append(m["category"])
        return seen

    def _raw(self, marker: str) -> dict:
        for m in self._markers:
            if m["marker"].lower() == marker.lower():
                return m
        raise KeyError(f"marker {marker!r} is not in {self.config_path.name}. "
                       f"Known markers: {', '.join(self.markers())}")

    # --------------------------------------------------------------- resolve
    def _resolve_gene(self, gene: str, offline: bool = False) -> dict:
        """Node metadata for a gene, or an explicitly unresolved marker.

        A failed lookup is never cached. ChEMBL returned 500s for hours during
        development, and an entry written in that window is indistinguishable
        from a real answer: zero structures and zero ligands, which triage reads
        as "no viable small-molecule route" and reports with full confidence.
        The cache therefore only ever holds answers the APIs actually gave.
        """
        if gene in self._cache:
            if _verified_node(gene, self._cache[gene]):
                return copy.deepcopy(self._cache[gene])
            # Also refuse entries inserted or mutated by a caller after load.
            del self._cache[gene]
        if gene in self._failed:
            return {"gene": gene, "resolved": False, "resolve_error": self._failed[gene]}
        if offline:
            return {"gene": gene, "resolved": False,
                    "resolve_error": "offline; node counts were never retrieved"}
        if self._out_of_time():
            # Not a finding about this gene. `triage` reads a missing count as
            # `counts_unresolved` and refuses to grade on it, which is the whole
            # point: a budget that let the remaining markers through as zeros
            # would turn "we stopped looking" into "no viable small-molecule
            # route" -- a verdict, reported with full confidence, produced by a
            # clock.
            self._skipped_for_budget.append(gene)
            return self._remember_failure(
                gene, f"the {self._budget_s:.0f}s lookup budget for this run "
                      "expired before this gene was reached; its counts were "
                      "never retrieved. Re-run to resolve it -- what did "
                      "resolve is cached, so a second run starts further along")
        info: dict = {"gene": gene, "resolved": True}
        try:
            identity = target_identity.resolve(gene, primary_gene_only=True)
            if not target_identity.validate_primary_binding(identity, gene):
                raise target_identity.TargetIdentityError(identity.get("chembl_mapping_note") or
                    "Primary gene and complete human target mapping were not verified.")
            info.update(uniprot=identity["uniprot"], protein_name=identity["protein_name"],
                        chembl_target=identity["chembl_target"], identity=copy.deepcopy(identity),
                        node_cache_schema=_NODE_CACHE_SCHEMA)
            info["n_pdb"] = chembl.pdb_count(info["uniprot"])
            tid = info["chembl_target"]
            info["n_ligands_p6"] = chembl.activity_count(tid) if tid else 0
            if not _verified_node(gene, info):
                raise ValueError("Node counts must be nonnegative integer observations for the verified protein.")
        except Exception as e:                       # network / API trouble
            log.warning("could not resolve %s: %s", gene, e)
            return self._remember_failure(gene, f"{type(e).__name__}: {e}")
        self._cache[gene] = copy.deepcopy(info)
        return copy.deepcopy(info)

    def _out_of_time(self) -> bool:
        """Has this run spent its network budget? Starts the clock on first use.

        Lazily, so constructing a Registry costs nothing and a run served
        entirely from cache never starts a clock it cannot exceed.
        """
        if self._budget_s is None:
            return False
        now = time.monotonic()
        if self._deadline is None:
            self._deadline = now + self._budget_s
            return False
        return now >= self._deadline

    @property
    def skipped_for_budget(self) -> list[str]:
        """Genes left unresolved because the run ran out of time, not because
        any API refused them. Disjoint in meaning from the rest of `unresolved`,
        though both appear there."""
        return list(self._skipped_for_budget)

    def _remember_failure(self, gene: str, why: str) -> dict:
        self._failed[gene] = why
        return {"gene": gene, "resolved": False, "resolve_error": why}

    @property
    def unresolved(self) -> dict[str, str]:
        """Genes this process could not resolve, and why."""
        return dict(self._failed)

    def _preserve_rejected_cache(self) -> None:
        if self._rejected_cache_bytes is None:
            return
        original = self._rejected_cache_bytes
        directory = self.cache_path.parent / "rejected"
        directory.mkdir(parents=True, exist_ok=True)
        archive = directory / (hashlib.sha256(original).hexdigest() + ".json")
        if archive.exists():
            if archive.read_bytes() != original:
                raise RuntimeError("Rejected node cache archive changed; refusing replacement")
            return
        fd, temporary = tempfile.mkstemp(dir=directory, prefix=".rejected-node-")
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(original)
                handle.flush()
                os.fsync(handle.fileno())
            # Do not replace an archive concurrently created by another writer.
            try:
                os.link(temporary, archive)
            except FileExistsError:
                if archive.read_bytes() != original:
                    raise RuntimeError("Rejected node cache archive changed; refusing replacement")
        finally:
            Path(temporary).unlink(missing_ok=True)

    def save_cache(self) -> None:
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        with _cache_write_lock(self.cache_path.with_name(self.cache_path.name + ".lock")):
            self._save_cache_locked()

    def _save_cache_locked(self) -> None:
        def check_current():
            current = self.cache_path.read_bytes() if self.cache_path.exists() else None
            if current != self._expected_cache_bytes:
                raise RuntimeError("Node cache changed since it was loaded; reopen it before saving")
        check_current()
        self._preserve_rejected_cache()
        clean = {g: e for g, e in self._cache.items() if _verified_node(g, e)}
        payload = json.dumps(clean, indent=1, ensure_ascii=False).encode()
        fd, temporary = tempfile.mkstemp(dir=self.cache_path.parent, prefix=".node-cache-")
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            check_current()
            os.replace(temporary, self.cache_path)
            self._expected_cache_bytes = payload
            self._rejected_cache_bytes = None
            self._cache = clean
        finally:
            Path(temporary).unlink(missing_ok=True)

    def spec(self, marker: str, resolve: bool = True, offline: bool = False) -> TargetSpec:
        raw = self._raw(marker)
        nodes: list[DruggableNode] = []
        for n in raw.get("nodes") or []:
            info = (self._resolve_gene(n["gene"], offline=offline) if resolve
                    else {"gene": n["gene"], "resolved": False,
                          "resolve_error": "node counts not looked up (resolve=False)"})
            nodes.append(DruggableNode(
                gene=n["gene"], uniprot=info.get("uniprot"),
                chembl_target=info.get("chembl_target"),
                n_pdb=int(info.get("n_pdb") or 0),
                n_ligands_p6=int(info.get("n_ligands_p6") or 0),
                rationale=n.get("rationale", ""),
                # Which way this node has to be pushed. It is declared on
                # `DruggableNode` with ten lines of comment and this is its
                # only producer, so leaving it out made every consumer outside
                # `pipeline._required_action` read None -- and that function
                # only got the right answer by reaching back into this class's
                # private `_raw`, a second source of truth for a field the
                # spec is supposed to carry.
                node_action=n.get("node_action") or None,
                goal_support={**copy.deepcopy(n.get('goal_support') or {}),
                              **({"target_identity": copy.deepcopy(info["identity"])}
                                 if info.get("identity") else {})},
                # The pChEMBL floor this node needs, when it needs one other
                # than the default. Carried on the spec for the same reason
                # `node_action` is: a consumer that has to reach into `_raw`
                # for it is a second source of truth waiting to drift.
                potency_floor=n.get("potency_floor"),
                potency_floor_why=n.get("potency_floor_why"),
                primary=bool(n.get("primary", False)),
                preferred=bool(n.get("preferred", False)),
                metalloenzyme=n.get("metalloenzyme"),
                structure_usable=bool(n.get("structure_usable", True)),
                structure_note=n.get("structure_note", ""),
                resolved=bool(info.get("resolved", True)),
                resolve_error=info.get("resolve_error", "")))
        notes = [raw["note"]] if raw.get("note") else []
        return TargetSpec(marker=raw["marker"], category=raw["category"],
                          direction=Direction(raw["direction"]),
                          mode=TargetMode(raw["mode"]), nodes=nodes, notes=notes)

    def specs(self, markers: Optional[Iterable[str]] = None, category: Optional[str] = None,
              resolve: bool = True, offline: bool = False) -> list[TargetSpec]:
        names = list(markers) if markers else self.markers(category)
        out = [self.spec(m, resolve=resolve, offline=offline) for m in names]
        if resolve and not offline:
            self.save_cache()
        return out

    # ----------------------------------------------------------- explanation
    @staticmethod
    def explain(spec: TargetSpec) -> str:
        """One paragraph a formulator can act on."""
        if spec.mode is TargetMode.ENDPOINT:
            return (f"{spec.marker} is not a molecular target. "
                    f"{spec.notes[0] if spec.notes else 'Use it as an efficacy readout.'}")
        if spec.mode is TargetMode.PROXY:
            genes = ", ".join(n.gene for n in spec.nodes)
            head = (f"{spec.marker} should {spec.direction.value}, but it has no usable "
                    f"binding pocket. Bind {genes} instead.")
            return f"{head} {spec.notes[0]}" if spec.notes else head
        head = f"{spec.marker} can be targeted directly ({spec.direction.value})."
        return f"{head} {spec.notes[0]}" if spec.notes else head
