# Public tutorial images

These images illustrate public tutorial examples in isolated Ubuntu 24.04 installations. Study databases, generated prediction files and result bundles are not included in this repository. Each installation starts with an empty study history.

## High-resolution 3D previews, 2026-10-21

| Files | Source and resolution |
|---|---|
| `21-molstar-bound-pocket.png` | Native Mol* 5.12.0 PNG export, **3840 × 2160**. Binding-site focus on the observed REA ligand, PDB **1CBS / model 1 / chain A / REA 200 / assembly 1**. The close view shows part of the receptor. |
| `22-molstar-bound-overview.png` | Native Mol* PNG export, **3840 × 2160**. The same observed structure, using **전체 구조** (whole structure). |

The installed backend was **0.3.0rc1**; the public UI came from commit `31e3b53cc3cccee9309b904fd314742cb6881d62`. The images show observed coordinates, not candidate docking or evidence of efficacy. They were exported directly at the stated resolution, without upscaling, retouching or generative editing.

## Browser walkthrough

JPEGs `01`–`20` provide the step-by-step walkthrough at approximately **1265 × 712** or **1253 × 705**. Most are the original 2026-10-20 captures. Four views were replaced on 2026-11-18 using public source `f1a77923bfb0671dc581b0d00e4560bca5831772` on **Ubuntu 24.04.5 LTS**, with **Chrome 154 headless**.

| Updated files | State shown |
|---|---|
| `08-structure-and-skin.jpg` | Input molecule and 2D structure, before catalog registration. |
| `09-admet-results.jpg` | Unassessed ADMET before catalog registration. Available stored predictions appear after starter preparation. |
| `19-candidate-filters.jpg` | Nicotinamide scaffold-search filters; one stored seed observation and zero new unique candidates. |
| `20-empty-filter-results.jpg` | Purchase-only filtering hides that observation while the saved count stays one. |

The updated capture environment used **Python 3.12** and **PostgreSQL 16** from conda because apt required an unavailable sudo password. This capture verifies the tutorial screens and catalog preparation; it is not an unchanged rerun of the documented sudo/apt installation. Candidate yields depend on catalog coverage and search conditions.

Every approved image is hash-pinned in [manifest.json](manifest.json) and listed in `PUBLIC_FILES.txt`. The source check rejects unreviewed binary files, study/runtime formats and private-content patterns. Prediction examples are generated and exported within each user's installation.
