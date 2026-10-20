# Browser walkthrough

[한국어](browser.ko.md) · [README](../README.md)

Open http://127.0.0.1:4317 after the server prints **`SkinPipe: http://127.0.0.1:4317`**. If the browser cannot connect yet, check the running terminal; database verification can take several minutes. The current interface has Korean labels. This guide pairs those labels with English descriptions.


Screenshots use public tutorial examples in a local **0.3.0rc1 starter installation**: Nicotinamide (compound), HSP90AA1 (protein) and PDB 1CBS (bound structure). Click an image to enlarge it. These are walkthrough captures, not bundled research datasets.

## 1. Choose an input

Open **새 피부 연구** (*New skin study*, top navigation). The four steps are input → confirmation → execution → results.

| Starting point | What to enter | What to check |
|---|---|---|
| **화합물을 알고 있어요** — Compound | A compound name, or choose **SMILES 분자 구조** and paste SMILES | Exact 2D structure, stereochemistry, salt and charge |
| **단백질을 알고 있어요** — Protein | A protein name, gene symbol or identifier | Official identity and species; human is the visible default |
| **결합 구조가 있어요** — Bound structure | PDB ID, or choose **PDB/mmCIF 파일** and upload a file | Model, protein chain, ligand instance and observed pocket |

Choose the input type **before** entering the value. Its input-method menu changes with the selected type. Click **입력 확인하기 →** (*Check input*). The button shows an in-progress message while the lookup runs.

If a name refers to a family, mixture or several structures, choose an exact item or provide an exact structure. A compound is not replaced automatically by a related protein. A lookup failure keeps you at the input step with a suggested next action. The original uploaded file is retained; a filesystem path is not treated as a command.

![Compound name input and Check input button](images/02-compound-input.jpg)

### SMILES example

Select **SMILES 분자 구조 (선택)** in the compound input-method menu and paste `NC(=O)c1cccnc1` (Nicotinamide). Use the same confirmation steps as the name example.

![SMILES input selected with Nicotinamide structure](images/11-smiles-input.jpg)

### Protein example

Choose **단백질을 알고 있어요**, enter `HSP90AA1`, and click **입력 확인하기 →**. The compound-name menu is replaced by a protein input field. Confirm **Homo sapiens / UniProt P07900** before continuing.

![Protein input field with HSP90AA1](images/12-protein-input.jpg)

![Confirmed HSP90AA1 species and UniProt identifier](images/13-confirm-protein.jpg)

### Bound-structure example

Choose **결합 구조가 있어요**, leave **PDB ID** selected, and enter `1CBS`. To use your own file, select **PDB/mmCIF 파일** in that menu and use the file chooser, then **입력 확인하기 →**.

![Bound structure input with public PDB ID 1CBS](images/14-bound-structure-input.jpg)

## 2. Confirm the target

Read **확인한 입력** (*Confirmed input*). Check the displayed identifiers and source records. For compounds, verify the 2D diagram and SMILES. For proteins, verify the species and identity. For bound structures, check the model, chain and bound ligand; use **구성이나 정확한 분자 구조를 다시 선택** (*Change components or chemical identity*) if necessary. If ligand chemistry is incomplete, provide the requested SMILES or SDF rather than accepting an unknown structure.

The observed bound ligand coordinates define the pocket in structure mode. A free conformer is shown separately from an observed or docked pose. The protein can be used as a validated structure even when its biological name is still unidentified.

Click **← 입력 수정** (*Edit input*) to go back, or **이 입력으로 연구하기 →** (*Use this input*) to continue. This does not start a calculation yet.

![Confirmed compound and 2D structure](images/03-confirm-compound.jpg)

For the public `1CBS` example, inspect **model 1 / chain A / REA / residue 200** and the assembly. Click **결합 화합물과 주변 포켓을 3D로 확인** to open the observed pocket. Use **구성이나 정확한 분자 구조를 다시 선택** to change the selection when several ligands are present.

![Selected ligand and chain in the public 1CBS structure](images/15-confirm-bound-structure.jpg)

![MolStar display of observed ligand and receptor, with pocket and whole-structure controls](images/16-observed-pocket-3d.jpg)

**3D controls:** drag to rotate, use the wheel/pinch to zoom, **결합 부위** to focus on the pocket, **전체 구조** to restore the whole protein, and **×** to close. This screenshot is an observed PDB ligand, not a new docking result.

## 3. Choose the analysis and limits

For a compound, choose:

- Similar candidates: search the installed catalog using the confirmed seed.
- Scaffold hopping: compare other scaffolds; same-scaffold analogues remain distinguishable. Acyclic compounds use their connectivity rather than an invented ring scaffold.
- Assess input only: evaluate the seed without a protein or docking step.

Protein mode searches target evidence and candidates; docking needs a confirmed structure and pocket. Bound-structure mode uses the selected ligand for hopping and keeps the receptor, coordinate frame and pocket fixed. A failed or unsupported redocking check holds the docking stage and appears in the results.

Under **피부 평가 조건** (*Skin conditions*), choose a research goal and optionally formulation, concentration and application site. Leave unknown values empty or select **아직 정하지 않음** (*Not decided*). Missing exposure conditions do not become a hidden default lotion/concentration. All input modes evaluate discovered compounds for skin properties and available toxicity predictions.

The routine default aims for **up to 200 representative candidates**, evaluates at most **1,000 catalog candidates**, and has a **2-hour study time limit**. These are limits, not guaranteed yields. Open **후보 수·필터·최대 처리 시간 조정** (*Candidate/filter/time settings*) to change them. Source availability, filtering, duplicates and docking/QC holds can reduce the final count. New structure generation is optional and off by default.

The additional-toxicity panel shows whether a local model is registered. Exact-structure DB predictions are reused first. A new structure can use the registered model; unsupported structures, model failures and exhausted budgets remain unassessed. They do not receive a zero/safe value.

Open **ChatGPT 설명과 자료 전달 설정 (선택)** (*Optional ChatGPT and data permissions*) to inspect the permitted data categories. Finally click **새 연구 실행** (*Run new study*) **once**. Submission has a busy state. If a request times out, check history first; do not immediately submit it again.

![First tutorial assessment mode and optional skin conditions](images/04-analysis-options.jpg)

![Calculation scope and the final Run study button](images/05-run-study.jpg)

## 4. Follow progress, cancel or resume

Once accepted, the study is saved in **내 연구 실행 기록** (*My study history*). A moving indicator appears during submission. Waiting and running use a blue panel with a spinner; completion changes the whole panel to green with a check mark and **연구 완료** (*Study complete*). Partial results and held work use amber, and failure uses red; text also names each state. Stage, elapsed time, completed counts and held/failed stages are reported as work proceeds.

Status refreshes every five seconds. A failed request or 20 seconds without a new response stops the spinner and marks the current status as unknown. Use **결과 갱신** (*Refresh results*) to check again. The system reduced-motion preference disables the rotation. **진행·결과 보기** (*View progress/results*) in the study history opens that study's result panel.

Changing menu pages or closing the browser does **not** cancel an accepted study. Before submission, changing the input type or leaving the input screen can discard the unsent form; return to step 1 and check the input again. Use the page's back buttons when editing. A browser Back action changes the page, not the saved study.

To cancel computation, use that study's cancel control and follow the displayed confirmation/status. Wait until cancellation/resource cleanup is acknowledged; a queued request is not yet a completed cleanup. Completed observations and logs remain. The worker stops its managed child processes; other users' jobs are not cancelled. A process that cannot be reclaimed is reported as a cleanup problem rather than silently marked finished.

Resume appears for resumable interrupted/cancelled work. It validates the saved input, catalog, model and calculation conditions before using completed work. Changed conditions may require a separate study. Do not start a duplicate study just because you changed screens. Use the installation's **Stop-SkinPipe.sh** when intentionally shutting down the whole local app.

![Blue progress panel with spinner, current stage and refresh time](images/06-study-progress.jpg)

The progress screenshot comes from a separate Nicotinamide similarity search. The input-only assessment may finish before you see an intermediate stage.

![Cancel confirmation with Continue and Cancel controls](images/17-cancel-confirmation.jpg)

**계속 분석하기** keeps the study running. **이 연구 분석 취소** confirms cancellation. The screenshot shows the confirmation prompt only; wait for the subsequent cleanup status after confirming.

## 5. Read results and structures

The summary separates new unique candidates, representatives, generated structures, seed/reference/control rows, excluded/held items and saved observations. A reference row is not another newly found candidate. A completed task can still contain a held docking stage or unassessed endpoints; inspect the stage details.

Open a candidate to see:

1. Name/ID, SMILES and its **2D** chemical diagram.
2. Search and scaffold comparison evidence.
3. Docking values and available **3D** pose, with its pose type.
4. Skin properties, warnings and toxicity/ADMET endpoints, with model/source status.
5. Purchasing/synthesis information and source records.

The 2D diagram shows bonding; the 3D viewer shows spatial coordinates. Rotate/zoom in the viewer. Where no saved pose is available, the interface states that instead of displaying a generated conformer as a bound pose. Close the viewer before opening another candidate on a small screen.

Click **전체 후보 표 보기** (*Open complete candidate table*) in a completed study to compare rows. Filters can be toggled independently. **All enabled conditions must match** (AND). **제외 필터 해제** (*Clear exclusion filters*) resets exclusion criteria; the snapshot, search and sorting controls are separate. The count distinguishes stored rows, matches and the display limit. Sorting changes order, not the underlying saved observations. Hide/unhide unassessed values deliberately; absent activity or toxicity is not a measured zero.

![Candidate comparison filters and counts from a separate Nicotinamide similarity search](images/19-candidate-filters.jpg)

This similarity-search example has **4 stored observations**, including the seed and a repeated structure; its study summary reports **2 new unique candidates**. The first “assess input only” example does not search for candidates.

![An enabled purchase-only filter gives zero matches while four observations remain saved](images/20-empty-filter-results.jpg)

If the table becomes empty, compare **전체 모집단** (*Stored population*) with **필터 적용 후** (*After filters*). Here the purchase-only condition hides all four rows. Turn it off or click **제외 필터 해제** to see them again. Filters do not delete saved results.

![Completed Nicotinamide assessment and result actions](images/07-completed-study.jpg)

![Two-dimensional compound structure and skin properties](images/08-structure-and-skin.jpg)

![ADMET outputs with scales and model version](images/09-admet-results.jpg)

## 6. Download Excel and raw data

Use **Excel 받기** (*Download Excel*) in the study history/result panel. Wait for export preparation and click the ready download. For several studies, mark their **Excel에 포함** (*Include in Excel*) checkboxes and click **선택한 연구를 Excel로 받기** (*Export selected studies*).

An **Excel 묶음 ZIP** (*Excel bundle ZIP*) contains the summary, complete candidate tables, compound worksheets, conditions and source data. Individual compound tabs are split into groups of at most 200; download/extract the whole bundle and follow its index. The same compound in different studies retains each observation and protocol. Excluded/failed/unassessed records and zero-candidate studies remain in the export. A visible table limit does not truncate the stored export scope.

Use the summary worksheet to choose a study and the candidate index to find a compound tab. The compound tabs include a structure image, key numbers, selection/exclusion reasons, skin/toxicity information and sources. Keep the JSON/SDF/conditions included with the export when sharing results privately. An export failure does not delete the study; correct the displayed issue and retry the export.

![The prepared ZIP save link must be clicked to save the export](images/10-excel-ready.jpg)

**Download check:** look in your browser’s Downloads folder for `SkinPipe-…-xlsx.zip`. Extract it, open `review-0001.xlsx`, then **후보 검토** → `C000001`. In the first assessment example there is one input row and one compound tab; searching for new candidates is a separate mode.

## Optional ChatGPT

Open **ChatGPT 연결 (선택)** (*Connect ChatGPT*) → **Continue with ChatGPT**. Complete sign-in in the account provider's page and return. Never paste passwords, session files or API keys into a study input. **연결 상태 갱신** (*Refresh connection*) checks the connection; **ChatGPT 연결 해제** disconnects it. Credentials are stored in the protected local installation, not exported in Excel.

- **Input guidance:** explains the verified identity and choices; it does not choose chemical identity for the resolver.
- **문헌 검토 — Literature review:** enter the subject, species and action direction; inspect retrieved papers and run the selected model. Claims link to source records/evidence sentences.
- **Result explanation:** summarizes stored candidate numbers, reasons and missing information when you permit those results to be sent.

Choose a model from the account's returned model list. The preference is `gpt-6-luna` with `max` reasoning; an explicitly returned provider default is the fallback when Luna is absent. If neither is available, choose manually. Check the reasoning option supported by that model. A listed model is not a successful call: access/limits are confirmed by an actual completed response. A timeout, quota or validation error affects the AI step; screening and exports remain available. There is no automatic conversion to a paid API request.

Public literature review is separate from permission to send private SMILES, uploaded structures or calculated results. Inspect the data-permission choices before running that stage. ChatGPT does not modify the input identity, pocket, toxicity values, scores or task status. Without ChatGPT, the platform provides the saved numerical tables and standard result text.

![Optional ChatGPT connection screen before account sign-in](images/18-optional-chatgpt.jpg)
