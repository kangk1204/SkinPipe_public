# Installation, catalog and model preparation

[한국어](catalog.ko.md) · [README](../README.md)

## Where files go

The one-command installer uses `~/SkinPipe`, outside the source repository. It refuses unrecognized existing folders and symbolic-link paths.

| Location | Contents |
|---|---|
| `~/SkinPipe/installation.json` | Source revision, wheel hash, installation choices and completed preparation state |
| `~/SkinPipe/source` | Public software source used for this installation |
| `~/SkinPipe/runtime/app/software` | Application Python environment |
| `~/SkinPipe/runtime/model-runtime` | Separate ADMET model environment |
| `~/SkinPipe/runtime/model-profile` | Model profile and identity |
| `~/SkinPipe/runtime/catalog` | Official downloads, normalized structures, predictions, audit and search indexes |
| `~/SkinPipe/runtime/app/state` | Private local configuration, PostgreSQL files and research outputs |
| `~/SkinPipe/*.log` and `~/SkinPipe/runtime/*.log` | Installation and model-runtime diagnostics |

The catalog and weights are obtained during installation, not through a GitHub data archive. Do not add these folders to a Git repository or upload them as a public release.

## Small starter or full library

The default **skin-starter** prepares the pinned **huskinDB** and **SkinPiX** source snapshots. It is meant for learning the workflow and quick local checks. The full source plan also includes ChEMBL, ChEBI, COCONUT, LOTUS, NPAtlas, NTP ICE collections, ZINC and Mcule. Exact source URLs, hashes and declared terms are in [sources.json](../skinpipe/skinpipe/platform/_catalog_assets/sources.json); ZINC artifact URLs are in the adjacent [artifact index](../skinpipe/skinpipe/platform/_catalog_assets/zinc-artifacts.json).

Each installer obtains data directly from those providers. Review the linked provider's conditions before using a source. Some providers restrict redistribution; installing their files does not grant permission to pass a combined database to someone else. Keep provider attribution with locally produced datasets. SkinPipe does not grant licenses for third-party data or model weights.

To choose the full library, first download the installer to your home directory, read its source and the source conditions, then run it with explicit options:

```sh
python3 -c 'import urllib.request; urllib.request.urlretrieve("https://raw.githubusercontent.com/kangk1204/SkinPipe_public/main/install.py", "SkinPipe-install.py")'
python3 SkinPipe-install.py --prefix ~/SkinPipe-full --library full --device auto --accept-source-terms
```

Allow at least **300 GiB** free space and a persistent terminal session. Full-library ADMET can take days on CPU; runtime varies with downloaded source size and available hardware. The default starter installer checks for **20 GiB**. Neither choice copies another user's existing research/database.

## What happens automatically

1. Resolve and save an immutable public GitHub revision; build the application wheel from that source.
2. Install Ubuntu packages, including PostgreSQL and CPU Vina, then create the private app environment.
3. Install pinned ADMET-AI dependencies and model files in a separate environment. Their identity is checked against the packaged model metadata.
4. Download source snapshots, verify hashes, ingest structures, calculate ADMET and retain unsupported/failed records distinctly.
5. Audit stored predictions, build global search indexes and register the completed local catalog/model with the app.
6. Start the loopback browser server. Open the displayed URL when it is ready.

With **auto**, a compatible NVIDIA driver is checked, then a real CUDA calculation is tested. If the GPU/runtime cannot be used, a separate CPU environment is installed; logs and any failed CUDA environment are retained. Explicit **cuda** fails with an explanation rather than claiming GPU use. Explicit **cpu** skips the GPU probe. No driver replacement or other-user process termination is performed.

The full-library ADMET runtime and a study's additional-toxicity model are recorded separately. Exact-structure cached values can be reused only with their input/model/preprocessing identity. New unsupported structures stay unassessed. A candidate is not assigned another molecule's toxicity by name or similarity.

## Progress and resuming

The terminal prints numbered stages and log paths. In a second terminal:

```sh
bash ~/SkinPipe/runtime/Catalog-Status.sh
```

The catalog status records download → ingest → database → ADMET → audit → indexes. After `Prepare-Catalog.sh` exists, re-run the same installation command with the **same prefix, library and device** to resume interrupted catalog stages. Already acknowledged stages are checked and reused. The catalog has its own exclusive build lock, and the installer refuses a second owner of the same folder.

If package installation fails before that launcher exists, read the reported log. The partial environment remains for diagnosis. Correct the dependency/network/storage error and choose a new `--prefix`; no partial environment or completed research is overwritten. If source/model/code settings differ, use a new catalog folder instead of modifying an in-progress plan.

Once installed, launch with `Start-SkinPipe.sh`; do not repeat full preparation for every study. The browser's **설치 상태** (*Installation status*) reports app/database/index/model/engine readiness. For details, inspect the installation logs and `model-runtime-selection.json` under `runtime`.

## Use a browser on macOS or Windows

SkinPipe listens on the Ubuntu server's loopback interface. Keep its start terminal open. On your Mac/Windows terminal, replace the placeholders with your own SSH account and host:

```sh
ssh -N -L 4317:127.0.0.1:4317 YOUR_USER@YOUR_UBUNTU_SERVER
```

Keep this tunnel open and browse to **http://127.0.0.1:4317** on the Mac/Windows computer. Do not enter the Ubuntu server's private filesystem path into your browser. Keep local port **4317**: the app validates its host and port, so forwarding local 4318 to server 4317 is rejected. If the tunnel reports “Address already in use”, use a browser on Ubuntu or free port 4317 only after identifying the conflicting local app. Do not disable authentication or expose the app's HTTP port publicly to work around a tunnel problem.

## Stop and update

Stop using your installation's `Stop-SkinPipe.sh`. Wait for its process cleanup message before copying or moving its state folder. Keep the catalog, model profile and research state together when backing up. Do not remove a PostgreSQL directory while the app is running.

This public installer creates a versioned source snapshot. Re-running it starts/resumes **that installation**; it does not silently replace its software with newer main. For a later software release, use a new prefix and follow that release's migration instructions. Preserve the old installation and its research until the new one is checked. Copying a new wheel over a live environment is not an update procedure.

## Common problems

| Message or symptom | Next action |
|---|---|
| Unsupported operating system/root user | Use regular-user Ubuntu 24.04 x86_64; use a browser tunnel from other systems |
| Folder exists/unrecognized | Choose a new `--prefix`; do not delete another installation |
| Download/hash error | Check network/provider availability; preserve the log and retry the unchanged preparation |
| Not enough space | Free space outside active research or choose a larger new volume/path |
| Another build is running | Follow the existing terminal/status; do not start another catalog build |
| GPU check fails under auto | CPU fallback is reported; inspect CUDA logs if you need GPU acceleration |
| Model or search index unavailable | Open installation status; complete catalog/model preparation before the study |
| No candidate survives | Inspect saved exclusions/holds and search scope; change filters or input in a separate study |
| Excel preparation fails | Read the export error and retry the export; saved study results remain |

The source guard and tests are described in [CONTRIBUTING.md](../CONTRIBUTING.md). No study inputs, example result files or model weights belong in a public issue attachment.
