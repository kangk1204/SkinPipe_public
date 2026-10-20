# Candidate docking workflow verification — 2026-10-21

The operating SkinPipe was tested through its real browser and CPU runtime.
The README now gives a CRABP2 / UniProt P29373 — 1CBS / chain A / REA example,
with representative target 3, catalog limit 50 and a 1,800-second study budget.
Purchase, structural-alert, toxicity and opposing-action filters were retained;
new-structure generation remained off. These are workflow-check limits, not
validated sampling or biological recommendations. Starter-catalog yields can
be smaller than those from the existing operating catalog.

## Why earlier studies showed zero candidate docking

Replaying the actual deployed selection code on 69 stored observation rows
reproduced all 13 prior selection results. Causes differed by input: no explicit
protein receptor, structural alerts, the same parent structure as a measured
reference, unspecified stereochemistry, procurement requirements, or a
compound-only mode that never requested receptor docking. No missing candidate
was replaced with a control or an unassessed structure.

A 2XDX / HSP90AA1 / WOE trial acquired 50 catalog rows and selected 21 for
comparison, but receptor preparation correctly stopped before Vina. The
original deposited terminal LYS A:224 contains only an observed N atom;
Meeko could not type that incomplete residue while retaining the atom. All
1,678 selected heavy atoms matched their deposited coordinates before typing.
The preservation guard rejected the loss of N. The original input, failure,
checkpoint and downloaded workbook were retained. No residue was trimmed and
no QC threshold was relaxed. Further HSP90 preparation/redocking failures were
retained as failures; a better lower-ranked pose was not substituted for a
failing top pose.

## Actual candidate docking

| Browser input | Catalog rows | New candidate Vina calculations | Candidate pose QC passed | Representative count | Observed-ligand redocking |
|---|---:|---:|---:|---:|---:|
| Protein CRABP2/P29373 with explicit 1CBS/A/REA | 50 | 3 | 3 | 3 of target 3 | 1, QC passed |
| Unmodified original 1CBS mmCIF, A/REA, confirmed ligand SMILES | 50 | 2 | 2 | 2 of target 3 | 1, QC passed |

Both runs recorded zero reused docking results. The protein run stopped at its
representative target; the bound-file run exhausted its eligible pool. Five
candidate calculations represent four distinct exact structures/parents across
the two studies, not five independent molecules. Their seeds, observed ligand
and measured references are excluded from candidate completion counts.

The original receptor's 1,091 heavy-atom coordinates were unchanged. Vina used
seed 42, exhaustiveness 16 and nine poses in the same frozen receptor/pocket.
The top observed-ligand pose had symmetry-corrected RMSD
1.9997109215 Å against the unchanged 2.0 Å threshold. This is close to the
threshold. The existing RDKit resonance-equivalent terminal-atom symmetry is
part of that metric; a stricter bond-order-only atom mapping gives about
2.093765 Å and is a different metric. The policy and comparison were not
changed for this check. Passing this example does not establish a robust
redocking benchmark, receptor protonation accuracy or biological binding.
Candidate QC used the installed PoseBusters 0.6.5. Missing or unassessed items
remain explicit; calculation/QC completion does not establish efficacy or
safety.

## Browser, Excel and changes

The browser showed live processing and green completion, then each study's
actual Excel ZIP was saved using the prepared-download link. The two ZIPs
contain 52 and 51 stored observation rows, including their different input and
reference roles. Every candidate identifier, name, SMILES and complete ADMET
JSON was compared against its workbook; each row has an individual compound
tab. Raw docking, protocol and QC records are retained alongside readable
summaries. These are workflow observations, not an independent usability study.

The result page now distinguishes candidate calculations from observed-ligand
redocking and its QC. It explains evidenced not-requested, held and empty-pool
outcomes without inferring failure from zero alone. Counts and the scientific
selection policy are unchanged.

One actual protein submission took 77 seconds because official structure and
ligand sources retried. The browser previously timed out after 20 seconds even
though the server later accepted the study. The exact submission route now
uses the same bounded 600-second allowance as input resolution; automatic POST
retry remains disabled. Other status/result-read deadlines stay unchanged.
Opening a verified saved result also clears an obsolete lost-submission notice. The 3D viewer also reads a routine calculation's recorded protocol seed when a legacy top-level seed is absent, preserving legacy precedence and a recorded zero.

The actual candidate-table walkthrough exposed a second delay: its snapshot
list transferred 52,626,939 bytes before rendering. An optional
`/v2/snapshots?view=summary` now returns the same verified list metadata in
17,633 bytes for 46 snapshots. Existing full-response bytes are unchanged.
The compact response preserves ownership filtering, integrity/quarantine state,
counts, targets and protocol selection; omitted detail stays in the original
manifest. The list UI uses this response, while full-manifest consumers retain
their existing endpoint. It still fetches the separate study list, which was
not compacted in this change.

Focused regressions cover separate counts/QC, uncertain or future summaries,
near-threshold RMSD display, escaping, the measured 77-second delay, a
240-second response, deadline expiry, lost headers/body, cancellation and no
duplicate POST. Development web checks and the offline build were rerun. No
scientific thresholds, engine settings, eligibility rules or ranking weights
were changed.

## Preservation and publication scope

All pre-existing 42 studies, 71 jobs and 1,830 research files retained their
original row/file hashes. The new checks added three studies and six jobs;
45 studies, 77 completed jobs and 2,305 files remained afterward. The research runs used the original backend and scientific code. After all
77 jobs completed with no live leases, only the API was restarted once to
apply the minimal snapshot-list projection. The installed patch passed all
seven focused tests in the actual Python environment before replacement, and
its full response was byte-identical afterward. Coordinator and PostgreSQL
process identities, configuration, model profile, frozen docking policy and
scientific execution sources were retained. Original API source and rollback
instructions were preserved. The existing backend version/base commit remains
recorded separately from this scoped source patch. The separate installed QA app also retained
its four studies, five completed jobs and 99 files while only the intended
browser assets were replaced with recoverable backups.

The public repository contains code, synthetic tests, documentation and prior
public-example screenshots. It does not contain the operating database,
original 31 studies, new study exports, candidate pose packages or credentials.
The new studies remain in the existing owner-private installation.
