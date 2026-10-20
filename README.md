# SkinPipe

[한국어 안내](README.ko.md) · [Browser guide](docs/browser.en.md) · [Installation and data](docs/catalog.en.md)

SkinPipe is a locally installed platform for compound screening in skin research. Start with a compound name or SMILES, a protein name, or a protein–ligand structure. Confirm the input, run a study, compare candidates, inspect structures and ADMET predictions, and download Excel workbooks.

This repository contains the software and installation instructions. Research inputs, study results, chemical databases, calculated predictions and model weights are kept in each user's installation folder. The platform starts with an empty study history and no predefined study targets.

## Install and start

Use a **regular user account on Ubuntu 24.04 x86_64**, with sudo access, Python 3 and an Internet connection. Allow **20 GiB of free space** for the starter installation. The installer requests the Ubuntu password in the terminal when installing system packages.

Paste this single command into the Ubuntu terminal:

```sh
python3 -c 'import urllib.request; exec(compile(urllib.request.urlopen("https://raw.githubusercontent.com/kangk1204/SkinPipe_public/main/install.py", timeout=30).read(), "SkinPipe installer", "exec"))'
```

It downloads the public source, installs the app and a separate ADMET model environment, obtains a small catalog from official providers, calculates ADMET, builds the search index and starts SkinPipe. Open **http://127.0.0.1:4317** in a browser when the terminal says **Ready**. Keep that terminal open.

**The first preparation takes longer than later launches.** Models and dependencies are downloaded once; ADMET and the index are calculated before the research screen becomes ready. Follow the terminal messages and log paths. With `--device auto` (the default), a compatible NVIDIA GPU is tested and used when available; otherwise ADMET runs on CPU. Docking uses the installed CPU Vina engine. The installer does not install or change GPU drivers.

The default catalog is a small **starter** for learning the workflow. It does not provide the search breadth of the full library. Full-library preparation requires at least **300 GiB** of free space and can take **days on CPU**. Read [data sources, terms, progress and full preparation](docs/catalog.en.md) before choosing it. There is no bundled precomputed database to download from this repository.

For macOS or Windows, install on an Ubuntu server and use your browser through the [SSH tunnel instructions](docs/catalog.en.md#use-a-browser-on-macos-or-windows). Native macOS/Windows installation is not provided by this installer.

## Your first study, in pictures

The screenshots below come from an installed **0.3.0rc1 starter environment**. Try **Nicotinamide** first: it needs no protein or docking setup. The interface currently uses Korean labels; the instructions give the exact button text. Click any image to enlarge it. Counts and times will depend on your catalog and computer.

### 1. Check that the installation is ready

After the terminal prints the app address, open **http://127.0.0.1:4317** and click **설치 상태** (*Installation status*, top right). Look for **PostgreSQL 연결됨** and a catalog structure count. An unconnected optional ChatGPT account is normal. Click **연구 화면으로** (*Go to research*) to continue.

![Installation status: database connected, starter catalog and optional ChatGPT](docs/images/01-installation-ready.jpg)

### 2. Enter a compound name

Choose **화합물을 알고 있어요** (*I have a compound*). Leave the input method as **화합물 이름** (*Compound name*), type `Nicotinamide`, then click **입력 확인하기 →** (*Check input*). Wait for the identity lookup; this has not started the study yet.

![Compound input: Nicotinamide and the Check input button](docs/images/02-compound-input.jpg)

### 3. Check the molecule before proceeding

Check **Nicotinamide**, **PubChem CID 936** and the **2D structure**. If it is the wrong molecule, use **← 입력 수정** (*Edit input*). Otherwise click **이 입력으로 연구하기 →** (*Use this input*).

![Confirmed Nicotinamide identity and molecular structure](docs/images/03-confirm-compound.jpg)

### 4. Choose a short first analysis and run it

Select **입력 화합물만 평가** (*Assess input only*) for this first example. Leave unknown skin/formulation/concentration fields empty. This evaluates one input molecule; it does not search for new candidates. Later, choose **비슷한 화합물 찾기** (*Find similar compounds*) or **다른 골격 후보 찾기** (*Find other scaffolds*) to screen the installed catalog.

![Assess input only selected, with optional skin conditions](docs/images/04-analysis-options.jpg)

Scroll to the bottom, check the calculation scope, then click **새 연구 실행** (*Run new study*) **once**. The button changes to **연구 접수 중…** while the request is being accepted.

![Calculation scope and the Run new study button](docs/images/05-run-study.jpg)

### 5. Follow the study and open its results

The saved study checks its status every five seconds. A short assessment may reach **처리 완료** (*Complete*) almost immediately. Expect **1 input molecule** and **no new-candidate search** in this example. For a longer search, the panel reports the current stage. You can leave the page and reopen the study using **진행·결과 보기** (*View progress/results*) in **내 연구 실행 기록** (*My study history*).

![Completed Nicotinamide assessment, one input molecule and result actions](docs/images/07-completed-study.jpg)

### 6. Inspect structures, skin properties and ADMET

Scroll down to the **2D structure**, skin properties and **기존 DB 독성·ADMET 예측** (*Cached toxicity/ADMET predictions*). The table shows each model output, scale and model version. Missing values remain unassessed. For a structure-based study, use the [3D pocket walkthrough](docs/browser.en.md#bound-structure-example) to inspect the observed ligand and receptor.

![ADMET outputs with endpoint names, units and model version](docs/images/09-admet-results.jpg)

### 7. Save Excel — there are two clicks

First click **Excel · 전체 후보·개별 물질 탭** (*Excel: complete table and compound tabs*), or **Excel 받기** in the history. After preparation finishes, click the new **준비된 Excel 묶음 ZIP 저장** (*Save prepared Excel bundle ZIP*) link. The first click prepares the file; the second saves it to your browser's Downloads folder.

![Prepared Excel bundle ZIP save link beside the Excel button](docs/images/10-excel-ready.jpg)

Extract the ZIP and open `review-0001.xlsx`. Start with **후보 검토** (*Candidate review*) and follow its links to individual compound tabs such as `C000001`. Keep the whole ZIP when sharing a study. A larger study may have several workbooks; its index links them together.

You have completed the first workflow when the study says **처리 완료** and the downloaded workbook opens. Continue with the [illustrated browser guide](docs/browser.en.md) for **SMILES, protein and PDB/file inputs**, running/cancelling/resuming, candidate comparison, 2D/3D and multiple-study exports.

## Start, stop and check progress

After installation:

```sh
bash ~/SkinPipe/runtime/Start-SkinPipe.sh
```

In a second terminal, to stop your local app and its managed jobs:

```sh
bash ~/SkinPipe/runtime/Stop-SkinPipe.sh
```

To inspect first-time catalog preparation:

```sh
bash ~/SkinPipe/runtime/Catalog-Status.sh
```

If catalog preparation was interrupted after `Prepare-Catalog.sh` was created, rerun the original installation command with the same options. It resumes the saved catalog rather than starting a second build. Dependency installation failures before that point retain their logs; the error explains how to retry with a new installation folder. Starting an already prepared installation does not recalculate its full ADMET cache.

In the browser, use **설치 상태** — *Installation status* to check the database, search index, model and docking engine. See [troubleshooting and file locations](docs/catalog.en.md).

## ChatGPT is optional

Core screening, local prediction, result viewing and Excel export work without ChatGPT. **ChatGPT 연결 (선택)** opens the account connection flow; **문헌 검토** reviews public literature. Input guidance and result explanations use the data categories you permit for that study.

The model selector prefers `gpt-6-luna` with `max` reasoning when the account advertises that choice. If Luna is absent, it uses an explicitly advertised provider default; if no default is returned, select an available model yourself. A failed request does not silently switch models or turn into a separately billed API request. Actual model access and usage limits depend on the connected account. See [ChatGPT setup and data controls](docs/browser.en.md#optional-chatgpt).

## Development

```sh
python3 scripts/check-public-tree.py
python3 -m unittest discover -s tests -v
python3 -m pip wheel --no-deps --wheel-dir /tmp/skinpipe-public-wheel ./skinpipe
```

Put installations and test outputs outside the checkout. Before contributing, read [CONTRIBUTING.md](CONTRIBUTING.md). Distribution licensing for SkinPipe is being decided; no project license has been selected. Bundled third-party components retain their own notices.
