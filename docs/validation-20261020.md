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
| Public Python regression at the dated verification | 61 passed, including 22 passing subtests, on Ubuntu with isolated PostgreSQL; this is a historical count, not the current tree's count |
| New browser progress tests | 18 passed in each source tree, including actual renderer state transitions |
| Development web checks | Passed, including 144 Python bridge tests; offline bundle and module graph rebuilt and checked |
| Fresh starter installation | 333 structures; 322 successful CPU ADMET predictions; 11 unsupported multicomponent structures retained as such |
| Browser workflow | Nicotinamide input confirmation, assessment, actual Excel ZIP download, and a separate similarity-search progress-to-completion transition |
| Downloaded assessment workbook in that environment | Five sheets, two embedded structure images and 755 nonempty cells; exports without endpoint definitions omit the ADMET-definition sheet and may have four sheets |

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

## Follow-up: published-source installation and operating UI

The README's one-line command was run unchanged in another pristine Ubuntu
24.04.5 x86_64 guest, as the regular `researcher` user with sudo available.
It fetched public source `3f1174656be3e64d7103c2002b3824650c18ddaa` and completed
the starter installation without a repair or retry. The default browser port
4317 returned HTTP 200. Preparation took approximately 5 minutes 45 seconds
in this guest; this is an observed duration, not an installation-time promise.
It again produced 333 structures, 322 successful CPU predictions and 11
unsupported multicomponent structures. Installed UI files matched the source.

Following the README, the browser showed installation readiness, resolved
Nicotinamide to PubChem CID 936, ran input-only assessment, displayed completion,
and saved the Excel ZIP using the documented preparation button and download
link. The actual downloaded workbook opened with five sheets, two structure
images and one candidate row. All 52 ADMET numeric values and their complete
JSON metadata matched the source row included in that ZIP; name and SMILES
also matched. No formula cells were present.

The new app was then stopped and restarted through its launchers. Its study,
snapshot, export and database rows were unchanged. Research/result files,
catalog files, model profiles and weights retained their hashes, sizes and
modification times. The browser reopened the completed study from its history.
For parallel browser QA only, this guest used the supported `--port 4396`
option after the default-port check; the existing user's port 4317 installation
was preserved. The sudo password prompt itself was outside this test because
the disposable guest's account had passwordless sudo.

The same progress UI was also published to the existing owner-private operating
Site. Its access policy and zero-schedule configuration were unchanged.
Source comparison found only the intended progress assets and generated
offline bundle changed; all other existing public assets retained their hashes.
Existing study rows and stored files were checked before and after publication
and a separately labelled test study. The operating API, coordinator and
PostgreSQL processes were unchanged.

At the operating URL, a preexisting result reopened with the green completion
panel. A new Nicotinamide assessment visibly transitioned from submission to
blue running status and green completion, then its Excel ZIP was downloaded.
The retained operating exporter produces four sheets, while the fresh public
installation produces five including the reading guide. Both layouts were
checked explicitly; both contain the candidate sheet, individual compound tab,
source contract, ADMET definitions and two structure images. Each ZIP's 52
ADMET values and full metadata matched its own exported source row.

The Site build and existing web checks passed, including 136 bridge tests and
the new 18 progress tests. Initial checkout-only checks required canonical
Python sources that are absent from the standalone Site repository; those
sources were supplied outside the Site checkout for testing. They were not
added to the deployment. These checks cover the stated Nicotinamide workflow,
not a study of first-time users or all possible research inputs.
