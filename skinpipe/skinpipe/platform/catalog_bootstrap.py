"""One-time official downloads -> owner-local ADMET -> verified read-only catalog.

Every mutation is under a new, marked folder. This CLI never resumes legacy
studies, writes to registered catalogs, or transfers the original mixed DB.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from .artifacts import atomic_bytes, canonical_json, sha256_file
from .identity import digest

ASSETS = Path(__file__).with_name("_catalog_assets")
STAGES = ("download", "ingest", "database", "admet", "audit", "indexes")
LABELS = {
    "download": "공식 제공처에서 다운로드",
    "ingest": "구조 검증·중복 정리",
    "database": "새 DB 구성",
    "admet": "로컬 ADMET 계산",
    "audit": "저장 예측·구조·출처 대조",
    "indexes": "검색 인덱스 구성",
}


def write(path, value):
    atomic_bytes(path, canonical_json(value))


def code_identity():
    parent = Path(__file__).parents[1]
    # Freeze every packaged Python dependency used by any build stage.
    files = list(parent.rglob("*.py"))
    return {str(p.relative_to(parent)): sha256_file(p) for p in sorted(files)}


def make_plan(library, profile, batch_size=None, device_policy="auto"):
    if device_policy not in {"auto", "cpu", "cuda"}:
        raise ValueError("Choose auto, cpu or cuda")
    if batch_size is None:
        batch_size = 10000 if library == "full" else 512
    from .research_toxicity import _profile
    from rdkit import rdBase

    if rdBase.rdkitVersion != "2025.09.5":
        raise ValueError(
            "Install skinpipe[catalog-build] with the pinned RDKit 2025.09.5"
        )
    specs = json.loads((ASSETS / "sources.json").read_bytes())
    if library == "skin-starter":
        specs = [s for s in specs if s["id"] in {"huskindb", "skinpix"}]
    elif library != "full":
        raise ValueError("Choose skin-starter or full")
    if not 1 <= batch_size <= 10000:
        raise ValueError("Batch size must be 1..10000")
    profile = Path(profile).absolute()
    value, identity, _, _, _ = _profile(json.loads(profile.read_bytes()))
    return {
        "schema": "skinpipe-owner-catalog-v1",
        "library": library,
        "sources": specs,
        "zinc_artifacts_sha256": sha256_file(ASSETS / "zinc-artifacts.json")
        if library == "full"
        else None,
        "profile_path": str(profile),
        "profile_sha256": sha256_file(profile),
        "model_identity": identity,
        "batch_size": batch_size,
        "device_policy": device_policy,
        "code": code_identity(),
        "download_space_gib": 4 if library == "skin-starter" else 300,
        "scope": "These source snapshots only; not a copy of the legacy DB or its assays/annotations",
        "scientific_validation": "NOT_ESTABLISHED",
    }


@contextmanager
def owned_root(root, plan):
    root = Path(root).absolute()
    # Reject existing donor folders and symlink aliases, including symlink parents.
    if root.resolve() != root or root.is_symlink():
        raise ValueError("Use a real new folder, not a symlink")
    marker = root / "owner-plan.json"
    if root.exists():
        if not marker.is_file() or marker.is_symlink():
            raise ValueError(
                "Existing unowned folder preserved. Choose a new catalog directory."
            )
        if json.loads(marker.read_bytes()) != plan:
            raise ValueError(
                "Source/model/code/settings changed. Preserve this catalog and use a new directory."
            )
    else:
        root.mkdir(parents=True, mode=0o700)
        write(marker, plan)
    with (root / ".build.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # Defense against redirecting a resumed owned folder into another DB.
        if any(p.is_symlink() for p in root.rglob("*")):
            raise ValueError("Symlinks inside a build folder are unsupported")
        yield root


def run_child(argv, log, progress=None, completion=None, *, allow_cuda=False):
    """Reap only this subprocess session, including its download/model children."""
    from .engines import terminate_owned

    environment = {
        **os.environ,
        **({} if allow_cuda else {"CUDA_VISIBLE_DEVICES": ""}),
        "OMP_NUM_THREADS": "2",
        "OPENBLAS_NUM_THREADS": "2",
        "MKL_NUM_THREADS": "2",
        "PYTHONUNBUFFERED": "1",
    }
    with Path(log).open("ab") as output:
        child = subprocess.Popen(
            [str(x) for x in argv],
            stdout=output,
            stderr=output,
            env=environment,
            start_new_session=True,
        )
        try:
            while child.poll() is None:
                if completion and completion[0].is_file():
                    if json.loads(completion[0].read_bytes()) != completion[1]:
                        raise ValueError(
                            "Worker completion differs from requested stage"
                        )
                    return
                if progress:
                    progress()
                time.sleep(0.5)
            if child.returncode:
                raise RuntimeError(f"Stage failed (exit {child.returncode}); see {log}")
            # Fast workers may exit before the polling loop reads their receipt.
            # Exit status alone must not commit an unacknowledged build stage.
            if completion:
                if not completion[0].is_file():
                    raise RuntimeError(f"Worker completion receipt is missing; see {log}")
                if json.loads(completion[0].read_bytes()) != completion[1]:
                    raise ValueError("Worker completion differs from requested stage")
        finally:
            terminate_owned(child)


def file_manifest(root, paths):
    return {
        str(p.relative_to(root)): {"bytes": p.stat().st_size, "sha256": sha256_file(p)}
        for p in sorted(paths)
        if p.is_file()
    }


def verify_manifest(root, manifest):
    for name, item in manifest.items():
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("Unsafe checkpoint file path")
        path = root / name
        if (
            not path.is_relative_to(root)
            or path.is_symlink()
            or not path.is_file()
            or path.stat().st_size != item["bytes"]
            or sha256_file(path) != item["sha256"]
        ):
            raise ValueError("Committed catalog file changed: " + name)


def stage_files(root, stage):
    if stage == "download":
        return [
            *(root / "raw").rglob("*"),
            *(root / "manifests/downloads").glob("*.json"),
        ]
    if stage == "ingest":
        return [
            *(root / "data/records").rglob("*"),
            *(root / "manifests/ingest").glob("*.json"),
        ]
    if stage == "database":
        return [root / "data/skinscout.duckdb", root / "reports/db_build.json"]
    if stage == "admet":
        return [
            *(root / "data/predictions").rglob("*"),
            *(root / "models").rglob("*"),
            root / "reports/admet-local.json",
            root / "manifests/catalog-execution.json",
        ]
    if stage == "audit":
        return [
            root / "reports/audit.json",
            *(root / "data/published_predictions").glob("*"),
            *(root / "data/models").glob("*"),
            *(root / "data/endpoints").glob("*"),
        ]
    return [*(root / "serving").rglob("*"), root / "local-catalog.json"]


def build(root, plan, app_root=None):
    with owned_root(root, plan) as root:
        marker_sha = digest(plan)
        progress_path = root / "build-status.json"
        for name in ["logs", "reports", "config", "src/skinscout", "checkpoints"]:
            (root / name).mkdir(parents=True, exist_ok=True)
        from .._catalog_source import audit

        shutil.copyfile(audit.__file__, root / "src/skinscout/audit.py")
        write(root / "config/sources.json", plan["sources"])
        write(
            root / "config/admet_policy.json",
            {"use": "auxiliary_metrics", "require_skin_validation": False},
        )
        completed = []
        active = "preflight"
        try:
            for number, stage in enumerate(STAGES, 1):
                active = stage
                checkpoint = root / "checkpoints" / (stage + ".json")
                if checkpoint.exists():
                    receipt = json.loads(checkpoint.read_bytes())
                    if receipt["plan_sha256"] != marker_sha:
                        raise ValueError("Checkpoint plan differs")
                    verify_manifest(root, receipt["files"])
                    completed.append(stage)
                    print(
                        f"{number}/6 {LABELS[stage]}: 기존 완료분 검증·재사용",
                        flush=True,
                    )
                    continue
                minimum = (
                    plan["download_space_gib"] if stage == "download" else 2
                ) * 1024**3
                if shutil.disk_usage(root).free < minimum:
                    raise ValueError(
                        f"Insufficient free space for {stage}: need {minimum // 1024**3} GiB"
                    )
                status = {
                    "status": "running",
                    "stage": stage,
                    "label": LABELS[stage],
                    "completed_stages": completed,
                    "plan_sha256": marker_sha,
                }
                write(progress_path, status)
                print(
                    f"{number}/6 {LABELS[stage]} 진행 중 — Ctrl+C로 중단, 같은 명령으로 재개",
                    flush=True,
                )
                last = [0.0]

                def report():
                    if time.monotonic() - last[0] < 15:
                        return
                    last[0] = time.monotonic()
                    detail_path = root / "reports" / (
                        "download-progress.json" if stage == "download" else "admet-local-progress.json"
                    )
                    detail = (
                        json.loads(detail_path.read_bytes())
                        if stage in {"admet", "download"} and detail_path.exists()
                        else {}
                    )
                    write(
                        progress_path,
                        {**status, "updated_at": time.time(), "detail": detail},
                    )
                    print(
                        f"  {LABELS[stage]} 실행 중 "
                        + (json.dumps(detail, ensure_ascii=False) if detail else ""),
                        flush=True,
                    )

                python = (
                    json.loads(Path(plan["profile_path"]).read_bytes())["python"][
                        "path"
                    ]
                    if stage == "admet"
                    else sys.executable
                )
                completion_path = root / "reports" / (stage + "-worker-complete.json")
                completion_path.unlink(missing_ok=True)
                run_child(
                    [
                        python,
                        "-m",
                        "skinpipe.platform.catalog_build_worker",
                        str(root),
                        stage,
                    ],
                    root / "logs" / (stage + ".log"),
                    report,
                    (completion_path, {"plan_sha256": marker_sha, "stage": stage}),
                    allow_cuda=stage == "admet" and plan["device_policy"] != "cpu",
                )
                write(
                    checkpoint,
                    {
                        "plan_sha256": marker_sha,
                        "files": file_manifest(root, stage_files(root, stage)),
                    },
                )
                completed.append(stage)
            audit_result = json.loads((root / "reports/audit.json").read_bytes())
            result = {
                "status": "complete_with_failed_predictions"
                if audit_result["failed_prediction_rows"]
                else "complete",
                "completed_stages": completed,
                "compound_count": audit_result["compound_count"],
                "prediction_counts": audit_result["prediction_counts"],
                "scope": plan["scope"],
                "catalog_config": str(root / "local-catalog.json"),
                "scientific_validation": "NOT_ESTABLISHED",
            }
            write(progress_path, result)
            if app_root:
                active = "registration"
                # Existing personal-local lock rejects running apps; no stop or job mutation.
                from .personal_local import main as research_main

                research_main(
                    [
                        "database",
                        "--root",
                        str(Path(app_root).absolute()),
                        "--scout-config",
                        str(root / "local-catalog.json"),
                    ]
                )
                research_main(
                    [
                        "toxicity",
                        "register",
                        "--root",
                        str(Path(app_root).absolute()),
                        "--profile",
                        plan["profile_path"],
                    ]
                )
                result["app_registration"] = "complete"
                write(progress_path, result)
            print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
            return result
        except BaseException as error:
            write(
                progress_path,
                {
                    "status": "interrupted"
                    if isinstance(error, KeyboardInterrupt)
                    else "needs_attention",
                    "stage": active,
                    "completed_stages": completed,
                    "error": str(error),
                    "next_action": "Fix the reported issue and rerun the same command. Completed inputs/predictions are retained.",
                },
            )
            raise


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="공식 출처 직접 다운로드 + 최초 1회 로컬 ADMET DB 구축"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ["plan", "build"]:
        p = sub.add_parser(name)
        p.add_argument(
            "--library", choices=["skin-starter", "full"], default="skin-starter"
        )
        p.add_argument("--profile", required=True, type=Path)
        p.add_argument("--batch-size", type=int)
        p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
        if name == "build":
            p.add_argument("--directory", required=True, type=Path)
            p.add_argument("--app-root", type=Path)
            p.add_argument("--accept-source-terms", action="store_true")
    p = sub.add_parser("status")
    p.add_argument("--directory", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.command == "status":
        print((args.directory / "build-status.json").read_text())
        return
    plan = make_plan(args.library, args.profile, args.batch_size, args.device)
    if args.command == "plan":
        print(json.dumps(plan, ensure_ascii=False, indent=2))
        return
    if not args.accept_source_terms:
        parser.error(
            "Review docs/local-catalog-build.ko.md and source licenses; then add --accept-source-terms for local use."
        )

    def interrupted(*_):
        raise KeyboardInterrupt("사용자 중단: 완료된 다운로드·예측을 보존했습니다.")

    signal.signal(signal.SIGTERM, interrupted)
    try:
        build(args.directory, plan, args.app_root)
    except KeyboardInterrupt:
        print("중단·자원 정리 완료. 같은 명령으로 이어갈 수 있습니다.", file=sys.stderr)
        return 130
    except (ValueError, RuntimeError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
