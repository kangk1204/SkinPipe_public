# Operating browser mode verification — 2026-10-21

The existing operating SkinPipe was exercised through its browser, with new,
labelled example studies. Its database, API and coordinator were retained.
The earlier pristine Ubuntu installation is recorded in
[the installation verification](validation-20261020.md); it was not repeated
for this UI-only follow-up. This is a bounded workflow check, not an independent
usability study or an exhaustive parameter-combination evaluation.

## Executed browser cases

| Input and analysis | Separate studies | Result |
| --- | ---: | --- |
| Nicotinamide name: similarity, scaffold comparison, input-only assessment | 3 | Completed; results reopened |
| Exact Nicotinamide SMILES: the same three analysis choices | 3 | Completed; results reopened |
| Protein: unspecified direction and inhibition (HSP90AA1), activation request (CRABP2 with matching 1CBS/A/REA receptor) | 3 | Completed; results reopened; unavailable candidate-level evidence retained |
| Bound complex: 1CBS PDB ID, original PDB upload, original mmCIF upload | 3 | Completed; results reopened; uploaded ligand chemistry explicitly confirmed |
| Compound scaffold comparison with new-structure generation enabled | 1 | Completed; results reopened |

Generation was also enabled in the mmCIF case. The six product choices and
four skin-location goals, including unspecified conditions, were distributed
across these examples. Their full Cartesian product was not evaluated.
For these small workflow checks the representative target was 2, catalog limit
3, and compute budget 360 seconds. Generation used BRICS limit 2, 3D comparison
limit 2 and one conformer per molecule. These are not adequate sampling
recommendations for a research campaign. Normal defaults were unchanged.

A newly created similarity study was cancelled through the confirmation dialog.
The browser reported process/resource cleanup, then resumed the same study
through its own confirmation flow. The second attempt completed. Existing
studies were not cancelled or resubmitted.

## Actual downloads and content verification

Protein and PDB-ID studies were individually downloaded as Excel ZIP files.
The 13 example studies were then selected in the updated operating UI and
exported together. The actual browser download contained a study index and
13 study workbooks with 69 stored observation rows and 138 embedded structure
images. These are stored observations, including seeds and references, not
69 independent new molecules.

Every workbook opened with openpyxl. Each study identity, candidate identifier,
name, SMILES and full ADMET JSON matched the source included in the same ZIP.
All 13 index links resolved to archive members, and the 454 manifest-listed
files matched their recorded SHA-256 and size. No formula cells were present
in the candidate workbooks. The batch ZIP SHA-256 was
`0c15940c5f9b695d92ea5be5e3954f58fd70633fc02f7552d473eb45d54b80a1`.
The operating exporter retains its four base sheet types; the newer fresh
installation additionally provides the reading guide. Individual C-number
tabs are added for each stored row.

## Findings repaired in both source trees and the operating UI

- A finished export notice could remain when another input or result opened.
  Notices and ready links now belong to their result; delayed responses cannot
  annotate a newly opened study. History-row exports retain their own context.
- Export execution waiting is distinguished from running file preparation.
- The generation control now says **BRICS recombination limit**. Functional-group
  substitutions are additional products, so the total can exceed that limit.
- Assembly labels describe the structure's recorded assembly without implying
  that an uploaded mmCIF was obtained through the public-PDB lookup path.
- README/browser instructions now cover protein direction, bound structure and
  file selection, confirmed ligand SMILES, generation, held/zero results and
  the two-click Excel download.

The deployed UI was reloaded and the stale-notice fix and BRICS label were
checked in the browser. Full development web checks passed, including 144
Python bridge tests; the operating Site checks passed with its 136 bridge
tests. The public export-state regression has 17 passing cases and the progress
regression has 18. The generated offline bundles and built module graph passed.

## Observed limitations

Two name-lookup attempts failed with connection/timeout errors. The browser
retained the failure and blocked unverified identity. The official service
subsequently returned Nicotinamide/CID 936, and all three name-input cases
completed. The documented exact-SMILES alternative was also exercised. No
identity was invented and no response was mocked.

The bound examples executed observed-ligand Vina redocking. The recorded RMSD
was 1.9997109 angstrom under the unchanged 2.0-angstrom gate. New-candidate
docking remained zero: structure alerts, undefined stereochemistry, procurement
requirements and other recorded eligibility conditions held the candidates.
Protein runs without an explicit receptor retained the corresponding docking
hold. The activation request is not evidence that any candidate activates the
protein. Generated structures and representative shortfalls remain visible
in results and exports; they were not converted into scientific successes.

These checks verify the tested input-to-export workflow. They do not establish
efficacy, safety, affinity, independent predictive performance, exhaustive
catalog coverage, optional GPU scorers, or all possible input/parameter cases.
Private existing research, database files, exports, raw run IDs, account/server
configuration and credentials are excluded from this public summary.

Before/after fingerprints confirmed every preexisting study/job row and stored
file was unchanged. API, coordinator and PostgreSQL process identities, settings
and model-profile hashes were unchanged. New labelled checks and their original
cancellation/retry/export records remain available. The operating Site access
policy was unchanged, and only the intended UI/test/build source files changed;
its existing stored research assets were preserved.
