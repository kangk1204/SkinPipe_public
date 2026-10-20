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

## Your first study

The current browser interface uses Korean labels. The English guide shows the exact labels beside their meanings.

1. Open **새 연구 시작** — *Start a new study*.
2. Choose **화합물을 알고 있어요** — *I have a compound*, **단백질을 알고 있어요** — *I have a protein*, or **결합 구조가 있어요** — *I have a bound structure*.
3. Enter your name, SMILES, PDB ID or structure file. Click **입력 확인하기** — *Check input*. Review the molecular structure, species and selected ligand/pocket before continuing. Ambiguous names require a choice.
4. Choose the search method and skin research conditions. Leave unknown formulation/concentration fields empty. Check the calculation limits, then click **새 연구 실행** — *Run study* once.
5. Watch the stage and elapsed-time messages. You may browse other pages after submission; the study is saved. Open **진행·결과 보기** — *View progress/results* from the study history to return to it.
6. Compare candidates, inspect **2D** structures and available **3D** poses, and open skin/toxicity details. Click **Excel 받기** — *Download Excel*. Large exports are delivered as an **Excel bundle ZIP**; extract it and open the workbooks.

See [the complete browser walkthrough](docs/browser.en.md) for all three input modes, cancel/resume, result counts, filtering, structures, multi-study exports and optional ChatGPT assistance.

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

In the browser, use **설치 상태 확인** — *Check installation* to check the database, search index, model and docking engine. See [troubleshooting and file locations](docs/catalog.en.md).

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
