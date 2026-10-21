# Public tutorial images

These images illustrate public tutorial examples in an isolated Ubuntu 24.04 starter installation. They do not publish a user's study database, result bundle or private research history. The software starts with an empty study history; screenshots illustrate studies submitted in that tutorial installation.

## High-resolution 3D and radar previews, 2026-10-21

| Files | Source and resolution |
|---|---|
| `21-molstar-bound-pocket.png` | Native Mol* 5.12.0 PNG export, **3840 × 2160**. Binding-site focus on the observed REA ligand, PDB **1CBS / model 1 / chain A / REA 200 / assembly 1**. The close view shows part of the receptor. |
| `22-molstar-bound-overview.png` | Native Mol* PNG export, **3840 × 2160**. The same observed structure, using **전체 구조** (whole structure). |
| `23-nicotinamide-properties.svg` | Native SkinPipe molecular-property SVG export for **Nicotinamide / PubChem CID 936**. Six axes retain raw values, scales and reference-policy metadata. |
| `24-nicotinamide-toxicity.svg` | Native SkinPipe toxicity SVG export for the same example, retaining six raw model outputs and their source. |

The installed backend is **0.3.0rc1**; the displayed public UI was copied from public commit `31e3b53cc3cccee9309b904fd314742cb6881d62`. The tutorial's saved prediction source is **ADMET-AI 2.0.1**, identified in each SVG's metadata. No new docking or model inference was performed to create these images. The Mol* images show observed coordinates, not candidate docking or evidence of efficacy.

PNG files were exported directly at the stated resolution, without upscaling, retouching or generative editing. SVGs preserve the app's exported geometry, text, values and metadata; one white-background CSS rule was added so labels remain readable in both light and dark README themes. Native-export and published-file SHA-256 hashes are recorded in [manifest.json](manifest.json).

The radar layout resembles a SwissADME-style profile, but its axes, ranges, property policy and model outputs are SkinPipe's. It is not an export from the SwissADME service. Amber is the highest component-score range in SkinPipe's skin-property policy for logP, molecular weight, TPSA and HBD; HBA and QED have no reference band. Display normalization does not modify scoring, QC or filters.

## Earlier browser walkthrough

Files `01`–`20` are the original 2026-10-20 browser JPEGs, approximately **1265 × 712**. They remain unchanged for step-by-step UI context; they are not described as 4K captures. The new PNG/SVG originals above provide the high-resolution 3D and plot previews.

Every approved image is individually hash-pinned in the manifest and listed in `PUBLIC_FILES.txt`. The public-source check continues to reject unreviewed binary files, study/runtime formats and private-content patterns.
