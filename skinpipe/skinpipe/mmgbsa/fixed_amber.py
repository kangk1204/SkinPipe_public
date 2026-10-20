"""Prospective, explicitly specified SQM charge production.

This producer is not yet the default charge path. Its numerical protocol must
pass the frozen realistic replay before replacing the adaptive OpenFF route.
The caller supplies the generated conformer as an SDF, independently of the
docked pose. Raw inputs, outputs and the normalized vector remain inspectable.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import time
from typing import Any

from .amber_identity import ambertools_manifest
from ..processes import run_calculation

PROTOCOL = "sqm-am1-fixed-diag6-pseudo0-v1"
SQM_KEYWORDS = ("qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, "
                "ndiis_attempts=700, diag_routine=6, pseudo_diag=0, ")
_METHODS = {"am1bcc": "bcc", "am1-mulliken": "mul"}


class FixedChargeError(RuntimeError):
    """No verified charge vector was produced."""


def vector_sha256(values: list[float]) -> str:
    """SHA-256 of atom-ordered little-endian float64 elementary charges."""
    if not values or not all(math.isfinite(v) for v in values):
        raise FixedChargeError("charges must be a nonempty finite vector")
    return hashlib.sha256(struct.pack("<" + "d" * len(values), *values)).hexdigest()


def normalize_charges(charges: list[float], formal_charge: int) -> tuple[list[float], float]:
    """Preserve the inspected OpenFF normalization's sequential float sum.

    A NumPy reduction uses a different operation order. Although the offset
    difference can be tiny, it must not be introduced into this controlled
    numerical protocol unnoticed.
    """
    vector_sha256(charges)
    current = 0.0
    for charge in charges:
        current += charge
    offset = (formal_charge - current) / len(charges)
    return [charge + offset for charge in charges], offset


def _run_command(argv: list[str], directory: Path, deadline: float,
                 name: str) -> dict[str, Any]:
    """Own SQM descendants on completion, timeout and Linux parent death."""
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise FixedChargeError("charge production exceeded its total timeout")
    started = time.time()
    with (directory / (name + ".log")).open("wb") as stream:
        try:
            result = run_calculation(argv, cwd=directory, stdout=stream,
                                     stderr=subprocess.STDOUT, timeout=remaining)
            code = result.returncode
        except BaseException as exc:
            # The supervisor owns cleanup. A timeout supplies no observed
            # command exit status, so retain null rather than inventing one.
            failed = {"argv": argv, "exit_code": None,
                      "started_at": started, "finished_at": time.time(),
                      "exception_type": type(exc).__name__,
                      "owned_process_group_terminated": os.name == "posix"}
            (directory / (name + "_receipt.json")).write_text(json.dumps(failed, indent=2) + "\n")
            raise
    row = {"argv": argv, "exit_code": code, "started_at": started,
           "finished_at": time.time()}
    (directory / (name + "_receipt.json")).write_text(json.dumps(row, indent=2) + "\n")
    if code:
        raise FixedChargeError(f"{name} exited with status {code}; see its saved log")
    return row


def _validate_atom_order(sdf: str, mol2: str) -> int:
    """Require sequential IDs, preserved coordinate order and connectivity.

    Amber's MOL2 writer rounds coordinates to three decimal places. Checking
    within that stated resolution plus every indexed edge prevents unnoticed
    permutations, rather than assuming that charge-file order is meaningful.
    """
    lines = sdf.splitlines()
    if len(lines) < 4 or "V2000" not in lines[3] or sdf.count("$$$$") != 1:
        raise FixedChargeError("exactly one V2000 SDF record is required")
    try:
        count, bonds = int(lines[3][:3]), int(lines[3][3:6])
        expected = [tuple(float(row[i:i+10]) for i in (0, 10, 20))
                    for row in lines[4:4+count]]
        symbols = [row[31:34].strip() for row in lines[4:4+count]]
        expected_edges = {tuple(sorted((int(row[:3]), int(row[3:6]))))
                          for row in lines[4+count:4+count+bonds]}
        sections = re.split(r"^@<TRIPOS>", mol2, flags=re.MULTILINE)
        atom_text = next(s for s in sections if s.startswith("ATOM\n"))
        bond_text = next(s for s in sections if s.startswith("BOND\n"))
        atoms = [row.split() for row in atom_text.splitlines()[1:] if row.strip()]
        bond_rows = [row.split() for row in bond_text.splitlines()[1:] if row.strip()]
        edges = {tuple(sorted((int(row[1]), int(row[2])))) for row in bond_rows}
        if count <= 0 or len(expected) != count or len(atoms) != count:
            raise ValueError("atom count")
        if len(expected_edges) != bonds or len(bond_rows) != bonds or edges != expected_edges:
            raise ValueError("indexed connectivity")
        for index, (atom, coords) in enumerate(zip(atoms, expected), 1):
            if int(atom[0]) != index or len(atom) < 9:
                raise ValueError("atom IDs")
            if not re.fullmatch(re.escape(symbols[index-1]) + r"\d+", atom[1], flags=re.IGNORECASE):
                raise ValueError("element label order")
            observed = [float(x) for x in atom[2:5]]
            if not all(math.isfinite(x) and abs(x-y) <= 0.00051
                       for x, y in zip(observed, coords)):
                raise ValueError("coordinate order")
    except (ValueError, IndexError, StopIteration) as exc:
        raise FixedChargeError(f"cannot verify the generated atom order: {exc}") from exc
    return count


def charge_sdf(input_sdf: Path, formal_charge: int, method: str,
               workdir: Path, timeout_s: float = 600.0) -> dict[str, Any]:
    """Produce a charge vector and retain all producer evidence in a new folder.

    This implements the same AM1/BCC model and thresholds as the inspected
    OpenFF AmberTools wrapper, with diag_routine=6 and pseudo_diag=0 explicitly
    specified. It does not perform ELF10 selection or charge the docking pose.
    """
    if method not in _METHODS:
        raise FixedChargeError(f"unsupported fixed SQM method: {method}")
    if isinstance(formal_charge, bool) or not isinstance(formal_charge, int):
        raise FixedChargeError("formal_charge must be an integer")
    if not math.isfinite(timeout_s) or timeout_s <= 0:
        raise FixedChargeError("timeout_s must be positive and finite")
    before = ambertools_manifest()
    if not before["verified"]:
        raise FixedChargeError("AmberTools producer identity is unverified")
    binary = shutil.which("antechamber")
    if binary is None:
        raise FixedChargeError("antechamber is unavailable")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=False)
    raw_input = Path(input_sdf).read_bytes()
    (workdir / "molecule.sdf").write_bytes(raw_input)
    protocol = {"protocol": PROTOCOL, "method": method, "formal_charge": formal_charge,
                "producer_before": before, "input_sdf_sha256": hashlib.sha256(raw_input).hexdigest(),
                "timeout_s": timeout_s, "scope": "prospective numerical charge protocol; no affinity validation"}
    (workdir / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    deadline = time.monotonic() + timeout_s
    argv = [binary, "-i", "molecule.sdf", "-fi", "sdf", "-o", "charged.mol2",
            "-fo", "mol2", "-pf", "no", "-dr", "n", "-c", _METHODS[method],
            "-nc", str(formal_charge), "-ek", SQM_KEYWORDS + f"qmcharge={formal_charge},"]
    commands = [_run_command(argv, workdir, deadline, "charge")]
    sqm_input = (workdir / "sqm.in").read_text()
    sqm_output = (workdir / "sqm.out").read_text()
    if (not re.search(r"\bdiag_routine\s*=\s*6\b", sqm_input) or
            not re.search(r"\bpseudo_diag\s*=\s*0\b", sqm_input) or
            "diag_routine=6)." not in sqm_output or
            "Auto diagonalization routine selection is in use." in sqm_output or
            "Calculation Completed" not in sqm_output):
        raise FixedChargeError("SQM output does not confirm the fixed completed protocol")
    count = _validate_atom_order(raw_input.decode(), (workdir / "charged.mol2").read_text())
    commands.append(_run_command([binary, "-dr", "n", "-i", "charged.mol2", "-fi", "mol2",
                                  "-o", "charges2.mol2", "-fo", "mol2", "-c", "wc",
                                  "-cf", "charges.txt", "-pf", "no"], workdir, deadline, "extract"))
    try:
        charges = [float(x) for x in (workdir / "charges.txt").read_text().split()]
    except (OSError, ValueError) as exc:
        raise FixedChargeError("no readable charge vector was produced") from exc
    if len(charges) != count:
        raise FixedChargeError(f"received {len(charges)} charges for {count} verified atoms")
    raw_digest = vector_sha256(charges)
    normalized, offset = normalize_charges(charges, formal_charge)
    after = ambertools_manifest()
    if not after["verified"] or after != before:
        raise FixedChargeError("AmberTools installation changed during charge production")
    receipt = {**protocol, "producer_after": after, "commands": commands,
               "atom_order_verified": True, "n_atoms": count,
               "raw_charges_e": charges, "raw_vector_sha256": raw_digest,
               "normalization_offset_e": offset, "normalized_charges_e": normalized,
               "normalized_vector_sha256": vector_sha256(normalized),
               "files": [{"path": p.name, "bytes": p.stat().st_size,
                          "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                         for p in sorted(workdir.iterdir()) if p.is_file()]}
    (workdir / "producer_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return receipt
