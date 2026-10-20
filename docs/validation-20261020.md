# Installation, browser and optional-runtime verification — 2026-10-20

This follow-up exercised the 35 previously skipped and 32 previously excluded
tests in the maintainer development suite. All **67 distinct original cases
passed**, with no remaining skip, exclusion, failure or error in that requested
set. Original scientific assertions and acceptance thresholds were unchanged.

## Execution scope

| Check | Observed result |
| --- | --- |
| Originally omitted cases | 35 skipped + 32 excluded → 67 executed and passed |
| Full development-suite selection | 6,027 distinct cases passed across the required environments and exact reruns |
| Full selection plus the 32 excluded cases | 6,059 distinct cases with final passing results |
| Public Python regression | 61 passed, including 22 passing subtests, on Ubuntu with isolated PostgreSQL |
| New browser progress tests | 18 passed in each source tree, including actual renderer state transitions |
| Development web checks | Passed, including 144 Python bridge tests; offline bundle and module graph rebuilt and checked |
| Fresh starter installation | 333 structures; 322 successful CPU ADMET predictions; 11 unsupported multicomponent structures retained as such |
| Browser workflow | Nicotinamide input confirmation, assessment, actual Excel ZIP download, and a separate similarity-search progress-to-completion transition |
| Downloaded assessment workbook | Five sheets, two embedded structure images and 755 nonempty cells |

The 6,059-case figure refers to the broader development suite, which includes
private historical-fixture checks and is not distributed in this public tree.
It is an environment-matrix result, **not one all-pass pytest invocation**.
The first 6,027-case run recorded 5,874 passes, 23 skips, 17 failures and 113
setup errors. Exact reruns resolved missing QA fixtures, an omitted test-DSN
port, the fixture database-creation permission, a required RDKit version and
the separate OpenFF environment. The original logs and per-case attempt
history were retained. No test outcome was converted to a pass without execution.

The fresh Ubuntu 24.04.5 x86_64 guest had eight virtual CPUs, 32 GiB memory
and no GPU. Separate environments were used for Linux lifecycle, catalog,
chemistry/ML and OpenFF calculations. The migrated-catalog checks used the
required RDKit 2025.9.5; the chemistry environment used RDKit 2026.3.6.
The existing application, research records and host services were preserved.
The temporary test role's database-creation permission was restored.

## Findings and repairs

- The installer could treat a partially written preparation launcher as a
  resumable installation. It now checks all six launchers and publishes the
  preparation entry point only after the other launchers are executable.
- “Ready” could appear before the server had finished checking its database.
  Installation preparation and the actual browser-ready message are now
  distinguished in the terminal and README.
- Input-only assessment displayed a contradictory 200-candidate routine
  description. That description is now hidden in assessment mode.
- Submission now displays a moving indicator immediately. Waiting and running
  use a blue progress panel; completion changes the entire panel to green with
  a check mark and “연구 완료”. Partial, failed, cancelled and unverified-cleanup
  states have distinct text and presentation.
- A failed status request or 20 seconds without a new response stops the
  spinner and marks the current state as unknown. No percentage or completion
  time is fabricated. Reduced-motion settings disable rotation.

A fresh reviewed-source installation completed all package, model and starter
catalog steps. The browser workflow and workbook download were checked against
that installation. The final three browser assets were then taken from the
reviewed wheel and atomically applied to the same QA installation while its
process and database remained active. Their installed hashes matched the wheel.
The updated screenshots in the README were captured from this application.

## Scientific limits

The optional runs included live public-service lookups, actual CPU Vina
redocking, OpenMM/OpenFF MM/GBSA, and chemistry/ML libraries. The 1T46/STI
redocking fixture yielded a symmetry-corrected pose RMSD of **0.352 Å**.
The MM/GBSA defaults selected **NAGL `openff-gnn-am1bcc-1.0.0.pt`, GAFF 2.11
and OpenMM 8.6.1 on CPU**; this does not establish an AmberTools QM charge run.
Large crystal-complex and truncated-ligand tests were executed. Some separate,
intentionally short minimization fixtures retained `converged=false`.

ML tests use synthetic software fixtures; they do not measure generalization
on independent biological data. These checks establish the tested software
behaviour and selected protocol reproducibility. They do not establish
biological efficacy, safety, experimental affinity, full-library ADMET coverage
or GPU performance.

## Evidence identity

The full local ledger preserves the original and repeated JUnit outcomes.
The parent verifier independently parsed every recorded attempt and checked
the source XML hashes.

- 6,027-case attempt ledger SHA-256:
  `4559bdf417760927287c40e60d8e8c23498456ebba177269e469d575da0a4d3a`
- Exact original 67-case verification ledger SHA-256:
  `8961e2b2045711f598b6a26d1b273bdf7e58f5836d5b7f5036a7c9d9a2504690`
- Final UI wheel used in the QA installation SHA-256:
  `54d8c238717a03a0ca3f7fc69173e4d55e6f6b8bbf78c7b168cf072ac9b22525`

Raw private fixtures, study databases, installation state and credentials are
not part of this public verification summary.
