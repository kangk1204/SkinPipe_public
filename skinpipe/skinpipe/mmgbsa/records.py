"""The two public shapes a worker record turns into.

``_run_single`` returns a plain dict so it can be pickled across a process
boundary and cached as JSON.  Callers never see that dict: ``mmgbsa_single``
hands back a :class:`~skinpipe.types.GBSAResult` and ``mmgbsa_batch`` a
DataFrame.  Both conversions live here so the column list and the field mapping
sit next to each other -- a column added to one and forgotten in the other is
how a reported number goes missing.

The frame is always built with the full column set even when a record lacks a
key, because a failed pose still occupies a row and a caller indexing
``frame["dG_bind"]`` must not get a ``KeyError`` on a batch where everything
failed.  Failures sort last (``na_position="last"``): a null dG at the top of a
table sorted best-first reads as the best hit.
"""
from __future__ import annotations

import math
import json
from typing import Any, Optional

from ..types import GBSAResult

#: Column order of the ``mmgbsa_batch`` DataFrame; also its complete schema.
DATAFRAME_COLUMNS = [
    "ligand_id", "pose_index", "dG_bind", "dG_std", "e_gas", "e_solv",
    "e_polar", "e_nonpolar", "e_complex", "e_receptor", "e_ligand", "converged",
    "n_frames", "rms_force_kj_nm", "max_force_kj_nm", "n_atoms_complex",
    "n_atoms_ligand", "wall_s", "sdf_path", "notes",
    "execution", "parameterization",
]

_RECEIPT_FIELDS = ("execution", "parameterization")


def _receipt(value: Any, name: str) -> tuple[Optional[dict[str, Any]], Optional[str]]:
    """JSON object or missing evidence, with no Python-expression decoding.

    Public objects carry dictionaries; batch/CSV columns carry canonical JSON
    text. Copies prevent a caller mutating the returned evidence from changing
    the underlying worker record. Malformed evidence is explicitly unverified.
    """
    if value is None or _is_nan(value):
        return None, None
    try:
        def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
            result = {}
            for key, item in pairs:
                if key in result:
                    raise ValueError("duplicate JSON key")
                result[key] = item
            return result
        if isinstance(value, str):
            def invalid_constant(token: str) -> None:
                raise ValueError(f"non-finite JSON constant {token}")
            value = json.loads(value, parse_constant=invalid_constant, object_pairs_hook=unique_object)
        if value is None:
            return None, None
        if not isinstance(value, dict):
            raise ValueError("expected a JSON object")
        # Parser recursion limits depend on the embedding host (including
        # pytest). The receipt contract has its own bound, independent of it.
        pending = [(value, 0)]
        nodes = 0
        while pending:
            item, depth = pending.pop()
            nodes += 1
            if depth > 64 or nodes > 100000:
                raise ValueError("receipt nesting or item budget exceeded")
            if isinstance(item, dict):
                pending.extend((child, depth + 1) for child in item.values())
            elif isinstance(item, (list, tuple)):
                pending.extend((child, depth + 1) for child in item)
        copied = json.loads(json.dumps(value, allow_nan=False), object_pairs_hook=unique_object)
        return copied, None
    except (TypeError, ValueError, OverflowError, RecursionError):
        return None, f"unverified {name} receipt: invalid JSON object"


def _with_receipt_notes(notes: Any, issues: list[Optional[str]]) -> str:
    text = _text(notes)
    for issue in issues:
        if issue and issue not in text:
            text = "; ".join(filter(None, [text, issue]))
    return text


#: Every ``GBSAResult`` float the worker record carries. Listed once and
#: splatted, rather than written out argument by argument, because the last
#: five of these were simply missing from the call and the single-pose API
#: reported ``None`` for numbers the batch frame had all along.
_FLOAT_FIELDS = (
    "dG_bind", "e_complex", "e_receptor", "e_ligand", "e_gas", "e_solv",
    "e_polar", "e_nonpolar", "dG_std", "wall_s",
)

#: The subset of :data:`_FLOAT_FIELDS` that is an *energy*. A record claiming
#: convergence must carry finite numbers in these; see
#: :func:`_non_finite_energies`. ``dG_std`` and ``wall_s`` are deliberately
#: absent: ``dG_std`` is ``None`` on every single-point run and ``NaN`` in that
#: run's frame row, which is normal rather than corrupt.
_ENERGY_FIELDS = (
    "dG_bind", "e_gas", "e_solv", "e_polar", "e_nonpolar",
    "e_complex", "e_receptor", "e_ligand",
)


def _number(value: Any) -> Optional[float]:
    """A float, or ``None`` when the record has nothing usable.

    ``NaN`` is never a measurement here: it is what a DataFrame puts in a
    failed pose's cell, and letting it through hands the best slot in the
    screen to a pose that did not compute. Missing and not-a-number are the
    same answer -- no value -- and both become ``None``.

    ``+/-Inf`` is the same answer again, and used not to be (D-11). ``NaN``
    was rejected here while ``Inf`` came through as a measurement, so a pose
    whose energy had blown up -- ``mean_std([-45.0, -44.0, inf])`` is
    ``(inf, nan)`` -- was carried as a converged result all the way to
    :func:`skinpipe.score.score_screen`, whose ``_finite`` rejects both and
    raised out of ``run_screen`` *after* docking and MM/GBSA had been paid
    for. This gate is now the same gate as ``score._finite``: ``None``, a
    bool, anything ``float()`` refuses, and any non-finite float are all
    "no value". ``tests/test_fix_c01_d11.py`` asserts the two agree value by
    value, so neither can drift away from the other again.
    """
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return out if math.isfinite(out) else None


def _non_finite_energies(record: dict[str, Any]) -> list[str]:
    """Energy fields the record carries but that are not finite numbers.

    "Carries but is not finite" and "does not carry" are different states: a
    record that never had ``e_polar`` is not corrupt, while one whose
    ``e_polar`` is ``Inf`` is.
    """
    return [name for name in _ENERGY_FIELDS
            if record.get(name) is not None and _number(record.get(name)) is None]


def _frames(value: Any) -> int:
    """Snapshot count, defaulting to the single-point 1.

    ``int(record.get("n_frames") or 1)`` looks equivalent and is not: ``NaN``
    is truthy, so it reached ``int()`` and raised ``ValueError``, which turned
    one failed pose in a batch into a dead mapping for every pose.
    """
    n = _number(value)
    if n is None or n < 1:
        return 1
    return int(n)


def _pose_index(value: Any) -> int:
    n = _number(value)
    return 0 if n is None else int(n)


def _flag(value: Any) -> bool:
    """``converged``, with ``NaN`` read as "did not converge".

    ``bool(nan)`` is ``True``. A failed pose's cell in the batch frame is
    ``NaN``, and reading it as converged would put a pose that never computed
    into the scoring term reserved for ones that did.
    """
    if isinstance(value, bool):
        return value
    if _is_nan(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1"}
    # CSV/numpy scalars can encode a boolean as 0 or 1. Other objects and
    # arbitrary positive numbers do not establish convergence.
    return _number(value) == 1.0


def _text(value: Any) -> str:
    """``notes``, with ``NaN`` read as empty rather than printed as "nan"."""
    if value is None or _is_nan(value):
        return ""
    return str(value)


def _is_nan(value: Any) -> bool:
    return isinstance(value, float) and math.isnan(value)


def to_result(record: dict[str, Any]) -> GBSAResult:
    """One worker record -- or one row of the batch frame -- as the public type.

    Every field the dataclass has is carried across. ``dG_std``, ``n_frames``,
    ``e_polar``, ``e_nonpolar`` and ``wall_s`` used to be dropped here, so the
    same computation read through ``mmgbsa_single`` and through
    ``mmgbsa_batch`` disagreed -- and, worse, ``GBSAResult.resolved`` is
    computed from sampling metadata: with the spread thrown away a 2.4
    kcal/mol within-run SD was misreported as below the legacy threshold.
    That flag does not measure ligand separation or repeated-run stability.

    A record that claims convergence and carries a non-finite energy is
    reported as a **failed pose** (D-11). That is what it is: there is no
    number to rank it by. Left as converged-with-no-dG it broke the contract
    ``score_screen`` enforces and raised out of ``run_screen``, and
    ``pipeline._attach_gbsa_records`` counted it among the poses it had
    "scored" on the way past.
    """
    values = {name: _number(record.get(name)) for name in _FLOAT_FIELDS}
    converged = _flag(record.get("converged", False))
    notes = _text(record.get("notes"))
    receipts = {name: _receipt(record.get(name), name) for name in _RECEIPT_FIELDS}
    notes = _with_receipt_notes(notes, [issue for _value, issue in receipts.values()])
    # Gated on ``dG_bind`` specifically, although every corrupt term is named:
    # this function also converts one *row of the batch frame*, where a column
    # another row introduced reads ``NaN`` here. That is an absent value, not a
    # corrupt one, and demoting a pose with a perfectly good dG_bind over it
    # would throw away a real result. ``dG_bind`` cannot be absent from a
    # record claiming convergence, and it is what the whole screen ranks on.
    corrupt = _non_finite_energies(record) if converged else []
    if "dG_bind" in corrupt:
        converged = False
        reported = ", ".join(f"{name}={record[name]!r}" for name in corrupt)
        notes = "; ".join(filter(None, [
            f"FAILED: non-finite MM/GBSA energy ({reported}); there is no "
            f"number to rank this pose by, so it is reported as a failed pose "
            f"rather than a converged one", notes]))
    return GBSAResult(
        ligand_id=str(record["ligand_id"]),
        pose_index=_pose_index(record.get("pose_index", 0)),
        n_frames=_frames(record.get("n_frames")),
        converged=converged,
        notes=notes,
        **{name: value for name, (value, _issue) in receipts.items()},
        **values,
    )


def to_frame(rows: list[dict[str, Any]]) -> "Any":
    """Batch rows with canonical JSON receipt columns; failures sort last."""
    import pandas as pd

    normalized = []
    for row in rows:
        copied = dict(row)
        issues = []
        for name in _RECEIPT_FIELDS:
            value, issue = _receipt(row.get(name), name)
            issues.append(issue)
            if value is not None:
                copied[name] = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
            elif issue is None:
                copied[name] = None
            else:
                # Keep malformed supplied material inspectable in the table,
                # while to_result refuses it as verified object evidence.
                copied[name] = row.get(name)
        copied["notes"] = _with_receipt_notes(row.get("notes"), issues)
        normalized.append(copied)
    frame = pd.DataFrame(normalized)
    for column in DATAFRAME_COLUMNS:
        if column not in frame.columns:
            frame[column] = None
    return frame[DATAFRAME_COLUMNS].sort_values(
        "dG_bind", na_position="last").reset_index(drop=True)
