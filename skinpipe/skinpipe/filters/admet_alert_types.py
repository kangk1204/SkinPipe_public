"""The alert record and the compile/match machinery shared by the alert tables.

Separated from `admet_alerts` only so the tables can be split across files
without a circular import: `admet_alerts` and `admet_photo_alerts` both import
from here, and `admet_alerts` re-exports everything so existing importers
(`admet`, `admet_irritation`) keep working unchanged.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable, Optional

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)


@dataclass(frozen=True)
class StructuralAlert:
    """One named substructure alert."""
    name: str
    smarts: str
    severity: str        # critical | high | moderate | low | context
    mechanism: str       # why this structure causes the endpoint
    action: str = ""     # the assay or restriction it points at


def compile_alerts(alerts: Iterable[StructuralAlert]
                   ) -> tuple[tuple[StructuralAlert, Chem.Mol], ...]:
    """Compile SMARTS, dropping (loudly) any pattern RDKit refuses.

    A pattern that fails to compile must never fail silently: it would look
    exactly like a compound with no alert.
    """
    out: list[tuple[StructuralAlert, Chem.Mol]] = []
    for a in alerts:
        patt = Chem.MolFromSmarts(a.smarts)
        if patt is None:
            log.error("admet: SMARTS for %r failed to compile: %s", a.name, a.smarts)
            continue
        out.append((a, patt))
    return tuple(out)


def match_alerts(mol: Chem.Mol,
                 compiled: Iterable[tuple[StructuralAlert, Chem.Mol]],
                 skip: Optional[set[str]] = None) -> list[dict]:
    """Return one record per firing alert, preserving table order."""
    skip = skip or set()
    hits: list[dict] = []
    for alert, patt in compiled:
        if alert.name in skip:
            continue
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue
        hits.append({"name": alert.name, "severity": alert.severity,
                     "mechanism": alert.mechanism, "action": alert.action,
                     "n_matches": len(matches)})
    return hits
