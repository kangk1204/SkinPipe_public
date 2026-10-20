"""Read-only ChEMBL SQLite adapter preserving individual assay conditions."""
from __future__ import annotations

from contextlib import closing
from pathlib import Path
import re
import sqlite3

from .artifacts import ArtifactIntegrityError, sha256_file


class LocalChEMBL:
    def __init__(self, database, *, source_version, source_sha256):
        self.database = Path(database).resolve()
        if not self.database.is_file() or not source_version or not source_sha256:
            raise ValueError("An existing ChEMBL database and source version/hash are required")
        self.source_version = source_version
        self.source_sha256 = source_sha256
        if not isinstance(source_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", source_sha256):
            raise ValueError("An exact lowercase SHA256 is required")
        self._verified_state = None
        self._verify()

    def _state(self):
        for suffix in ("-wal", "-journal"):
            sidecar = Path(str(self.database) + suffix)
            if sidecar.exists() and sidecar.stat().st_size:
                raise ArtifactIntegrityError("ChEMBL source contains an unpinned live transaction journal")
        info = self.database.stat()
        return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns

    def _verify(self):
        try:
            state = self._state()
            if state != self._verified_state:
                if sha256_file(self.database) != self.source_sha256 or self._state() != state:
                    raise ArtifactIntegrityError("ChEMBL source checksum changed or differs from the pinned database")
                self._verified_state = state
            return state
        except OSError as error:
            raise ArtifactIntegrityError("ChEMBL pinned source is unavailable") from error

    def activities(self, target_chembl_id, *, limit=10000, offset=0):
        """No potency/trust filter and no aggregation across assays or relations."""
        if not isinstance(target_chembl_id, str) or not target_chembl_id:
            raise ValueError("target_chembl_id required")
        if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 20000:
            raise ValueError("ChEMBL limit must be 1..20000")
        if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
            raise ValueError("ChEMBL offset must be nonnegative")
        state = self._verify()
        # URI read-only mode prevents accidental source database creation.
        uri = self.database.as_uri() + "?mode=ro"
        with closing(sqlite3.connect(uri, uri=True)) as connection:
            connection.row_factory = sqlite3.Row
            # Common ChEMBL releases add columns over time. Preserve every
            # activity and assay field under explicit namespaces rather than
            # selecting only a narrow potency summary.
            columns = [row[1] for row in connection.execute("PRAGMA table_info(assays)")]
            select = ",".join('s."' + column.replace('"', '""') + '" AS "assay_' + column.replace('"', '""') + '"' for column in columns)
            # ORDER BY activity_id otherwise encourages a full activities scan
            # in ChEMBL 37. Materialize the target assays first, then pin the
            # outer loop with CROSS JOIN so the assay FK lookup is used.
            sql = """WITH selected_assays AS MATERIALIZED (
                SELECT s.assay_id FROM assays s JOIN target_dictionary t ON t.tid=s.tid
                WHERE t.chembl_id=?)
                SELECT a.*,m.chembl_id AS molecule_chembl_id,t.chembl_id AS target_chembl_id,
                t.organism AS target_organism,c.canonical_smiles,""" + select + """
                FROM selected_assays selected CROSS JOIN activities a ON a.assay_id=selected.assay_id
                JOIN assays s ON s.assay_id=a.assay_id
                JOIN target_dictionary t ON t.tid=s.tid JOIN molecule_dictionary m ON m.molregno=a.molregno
                LEFT JOIN compound_structures c ON c.molregno=a.molregno
                ORDER BY a.activity_id LIMIT ? OFFSET ?"""
            rows = [dict(row) for row in connection.execute(sql, (target_chembl_id, limit, offset))]
        if self._verify() != state:
            raise ArtifactIntegrityError("ChEMBL source changed during the query")
        for row in rows:
            row["source_provenance"] = {"source": "ChEMBL SQLite", "version": self.source_version,
                                        "sha256": self.source_sha256}
        return rows

    def target_identity(self, target_chembl_id):
        """Verify release-specific components before using its measurements."""
        result = self.target_record(target_chembl_id)
        if result['target'] is None or not result['complete']:
            raise ValueError('ChEMBL component identity unavailable')
        target = result['target']
        return {'taxon_id': target['tax_id'], 'target_type': target['target_type'],
                'accessions': [c['accession'] for c in target['target_components']]}

    def _target_records(self, where, parameters):
        """Internal fixed-SQL helper; no caller-supplied SQL or source writes."""
        state = self._verify()
        with closing(sqlite3.connect(self.database.as_uri() + '?mode=ro', uri=True)) as connection:
            rows = connection.execute('''SELECT t.chembl_id,t.tax_id,t.target_type,t.organism,s.accession
                FROM target_dictionary t LEFT JOIN target_components c ON c.tid=t.tid
                LEFT JOIN component_sequences s ON s.component_id=c.component_id WHERE ''' + where +
                ' ORDER BY t.chembl_id,s.accession', parameters).fetchall()
        if self._verify() != state:
            raise ArtifactIntegrityError('ChEMBL mapping changed during query')
        grouped = {}
        for chembl_id, taxon, target_type, organism, accession in rows:
            target = grouped.setdefault(chembl_id, {'target_chembl_id': chembl_id, 'tax_id': taxon,
                'target_type': target_type, 'organism': organism, 'target_components': [], 'identity_complete': True})
            if (target['tax_id'], target['target_type'], target['organism']) != (taxon, target_type, organism):
                raise ArtifactIntegrityError('ChEMBL target ID has conflicting release-specific identities')
            if not accession:
                target['identity_complete'] = False
            elif accession not in {c['accession'] for c in target['target_components']}:
                target['target_components'].append({'accession': accession})
        return list(grouped.values())

    def _target_provenance(self):
        return {'source': 'ChEMBL SQLite', 'version': self.source_version, 'sha256': self.source_sha256}

    def target_record(self, target_chembl_id):
        """An exact local ID can supply components, which still require UniProt verification."""
        if not isinstance(target_chembl_id, str) or not re.fullmatch(r'CHEMBL[0-9]+', target_chembl_id):
            raise ValueError('An exact ChEMBL target ID is required')
        records = self._target_records('t.chembl_id=?', (target_chembl_id,))
        return {'target': records[0] if records else None,
                'complete': not records or records[0]['identity_complete'],
                'source_provenance': self._target_provenance()}

    def target_mappings(self, accessions, *, taxon_id, kind='protein'):
        """Exact verified accession SET, taxon and target scope; never a first match.

        ChEMBL complexes describe component sets, not an independently verified
        assembly/stoichiometry. The receptor still needs its own assembly checks.
        """
        from .targets import ACCESSION
        if (kind not in {'protein', 'complex'} or isinstance(taxon_id, bool)
                or not isinstance(taxon_id, int) or taxon_id < 1
                or not isinstance(accessions, (list, tuple, set)) or not 1 <= len(accessions) <= 16
                or any(not isinstance(a, str) or not ACCESSION.fullmatch(a) for a in accessions)):
            raise ValueError('Verified accessions, NCBI taxonomy and protein/complex kind are required')
        expected = sorted(set(accessions))
        if len(expected) != len(accessions) or kind == 'protein' and len(expected) != 1:
            raise ValueError('Exact unique target component set is required')
        target_type = 'SINGLE PROTEIN' if kind == 'protein' else 'PROTEIN COMPLEX'
        placeholders = ','.join('?' for _ in expected)
        records = self._target_records('''t.tid IN (SELECT DISTINCT x.tid
            FROM target_dictionary x JOIN target_components c ON c.tid=x.tid
            JOIN component_sequences s ON s.component_id=c.component_id
            WHERE x.tax_id=? AND x.target_type=? AND s.accession IN (''' + placeholders + '))',
            (taxon_id, target_type, *expected))
        matches, unresolved = [], []
        for target in records:
            if not target['identity_complete'] or any(not ACCESSION.fullmatch(c['accession']) for c in target['target_components']):
                unresolved.append(target['target_chembl_id'])
            elif [c['accession'] for c in target['target_components']] == expected:
                matches.append(target)
        return {'targets': matches, 'complete': not unresolved, 'unresolved_target_ids': unresolved,
                'query': {'accessions': expected, 'taxon_id': taxon_id, 'target_type': target_type},
                'source_provenance': self._target_provenance()}

    def activities_page(self, target_chembl_id, *, limit=10000, offset=0):
        """Expose pagination completeness without treating missing values as inactive."""
        rows = self.activities(target_chembl_id, limit=limit, offset=offset)
        state = self._verify()
        with closing(sqlite3.connect(self.database.as_uri() + "?mode=ro", uri=True)) as connection:
            total, recorded = connection.execute("""SELECT count(*),count(a.standard_value)
                FROM activities a JOIN assays s ON s.assay_id=a.assay_id
                JOIN target_dictionary t ON t.tid=s.tid JOIN molecule_dictionary m ON m.molregno=a.molregno
                LEFT JOIN compound_structures c ON c.molregno=a.molregno
                WHERE t.chembl_id=?""", (target_chembl_id,)).fetchone()
        if self._verify() != state:
            raise ArtifactIntegrityError("ChEMBL source changed during count query")
        return {"rows": rows, "target_chembl_id": target_chembl_id, "total_records": total,
                "records_with_standard_value": recorded, "offset": offset, "limit": limit,
                "next_offset": offset + len(rows) if offset + len(rows) < total else None,
                "complete": offset + len(rows) >= total,
                "count_scope": "raw source records; value presence is not measurement validity or activity",
                "source_provenance": {"source": "ChEMBL SQLite", "version": self.source_version,
                                      "sha256": self.source_sha256}}
