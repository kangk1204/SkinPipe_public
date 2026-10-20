#!/usr/bin/env bash
set -euo pipefail
# Fresh app; official data/model downloads; no original DB redistribution.
# Functions are sourceable for isolated installer branch tests.

task_new_install_directory() {
  local task_directory
  task_directory="$(realpath -ms "$1")" || return
  if [[ -e "$task_directory" || -L "$task_directory" || "$(realpath -m "$1")" != "$task_directory" ]]; then
    echo '기존 폴더와 심볼릭 링크를 보존합니다. 실제 새 폴더를 지정하세요. 중단한 DB 구축에는 Prepare-Catalog.sh를 사용하세요.' >&2
    return 2
  fi
  printf '%s\n' "$task_directory"
}

task_nvidia_driver_query() {
  command -v nvidia-smi >/dev/null 2>&1 || return 1
  timeout 15s nvidia-smi --query-gpu=driver_version --format=csv,noheader,nounits
}

task_select_model_runtime() {
  local task_requested="$1" task_drivers task_driver task_major
  task_runtime=cpu
  task_device=cpu
  task_runtime_reason='CPU 실행을 선택했습니다.'
  if [[ "$task_requested" == cpu ]]; then return 0; fi
  if [[ "$task_requested" != auto && "$task_requested" != cuda ]]; then
    echo '장치 선택은 auto, cpu, cuda 중 하나여야 합니다.' >&2
    return 2
  fi
  task_runtime_reason='사용할 수 있는 NVIDIA GPU 드라이버가 없어 CPU로 진행합니다.'
  if task_drivers="$(task_nvidia_driver_query 2>/dev/null)" && [[ -n "$task_drivers" ]]; then
    task_runtime_reason='CUDA 12.8에는 NVIDIA 드라이버 570 이상이 필요하여 CPU로 진행합니다.'
    while IFS= read -r task_driver; do
      task_driver="${task_driver//[[:space:]]/}"
      task_major="${task_driver%%.*}"
      if [[ ! "$task_major" =~ ^[0-9]+$ ]] || (( 10#$task_major < 570 )); then
        break
      fi
      task_runtime=cuda
    done <<< "$task_drivers"
    # Every reported device must have a compatible driver.
    if [[ ! "$task_major" =~ ^[0-9]+$ ]] || (( 10#$task_major < 570 )); then task_runtime=cpu; fi
    if [[ "$task_runtime" == cuda ]]; then
      task_device="$task_requested"
      task_runtime_reason='NVIDIA GPU와 호환 드라이버를 확인했습니다. CUDA 설치 후 실제 계산 가능 여부를 확인합니다.'
      return 0
    fi
  fi
  if [[ "$task_requested" == cuda ]]; then
    echo 'CUDA 실행을 요청했지만 호환되는 NVIDIA GPU를 확인하지 못했습니다. nvidia-smi와 드라이버 570 이상을 확인하거나 마지막 인수를 cpu로 바꿔 새 폴더에 설치하세요.' >&2
    return 2
  fi
}

task_install_model_environment() {
  local task_kind="$1" task_requirements task_constraints task_index
  task_requirements="$task_assets/requirements-model-hashes.txt"
  task_constraints="$task_assets/model-constraints.txt"
  task_index=https://download.pytorch.org/whl/cpu
  if [[ "$task_kind" == cuda ]]; then
    task_requirements="$task_assets/requirements-model-cu128-hashes.txt"
    task_constraints="$task_assets/model-cu128-constraints.txt"
    task_index=https://download.pytorch.org/whl/cu128
  fi
  python3.12 -m venv "$task_root/model-runtime" || return
  "$task_model_python" -I -m pip install --only-binary=:all: --require-hashes \
    --extra-index-url "$task_index" -r "$task_requirements" || return
  "$task_model_python" -I -m pip install --only-binary=:all: --require-hashes \
    -c "$task_constraints" -r "$task_assets/app-cli-supplement-hashes.txt" || return
  "$task_model_python" -I -m pip install --only-binary=:all: -c "$task_constraints" \
    duckdb==1.5.0 pyarrow==23.0.1 || return
  "$task_model_python" -I -m pip install --no-deps "$task_wheel" || return
  "$task_model_python" -I -m pip check || return
}

task_probe_cuda_runtime() {
  # A driver listing alone does not prove that this process can use the GPU.
  timeout 90s "$task_model_python" -I - <<'PY'
import json
import torch

if torch.__version__ != '2.8.0+cu128':
    raise RuntimeError('Pinned PyTorch CUDA 12.8 runtime is not installed')
if not torch.cuda.is_available():
    raise RuntimeError('CUDA is not available to this process')
x = torch.tensor([1.0, 2.0], device='cuda')
if float((x * x).sum().item()) != 5.0:
    raise RuntimeError('CUDA numerical smoke test failed')
torch.cuda.synchronize()
print(json.dumps({'torch': torch.__version__, 'cuda': torch.version.cuda,
                  'device': torch.cuda.get_device_name(0), 'probe': 'passed'}))
PY
}

task_prepare_model_runtime() {
  local task_requested="$1" task_failure=''
  task_preserved_cuda_attempt=''
  if [[ "$task_runtime" == cuda ]]; then
    echo "CUDA 모델 패키지를 설치합니다. 진행 기록: $task_root/model-cuda-install.log"
    if task_install_model_environment cuda 2>&1 | tee "$task_root/model-cuda-install.log"; then
      if task_probe_cuda_runtime > "$task_root/model-cuda-probe.log" 2>&1; then
        task_runtime_reason='실제 GPU 계산 확인을 통과했습니다. 전체 라이브러리 ADMET 계산에서 GPU를 사용할 수 있습니다.'
        echo "$task_runtime_reason"
        return 0
      fi
      task_failure='CUDA 패키지는 설치됐지만 실제 GPU 계산 확인에 실패했습니다.'
    else
      task_failure='CUDA 패키지 설치에 실패했습니다.'
    fi
    if [[ "$task_requested" == cuda ]]; then
      echo "$task_failure 설치 기록: $task_root/model-cuda-install.log / GPU 확인 기록: $task_root/model-cuda-probe.log. NVIDIA 드라이버·네트워크·여유 공간을 확인하거나 cpu로 새 폴더에 설치하세요. 기존 설치 자료는 보존했습니다." >&2
      return 1
    fi
    task_preserved_cuda_attempt="$(mktemp -d "$task_root/model-runtime-cuda-attempt.XXXXXXXX")" || return
    if [[ -e "$task_root/model-runtime" || -L "$task_root/model-runtime" ]]; then
      mv "$task_root/model-runtime" "$task_preserved_cuda_attempt/environment" || return
    fi
    task_runtime=cpu
    task_device=cpu
    task_runtime_reason="$task_failure 새 CPU 환경으로 전환합니다. CUDA 설치 자료와 기록은 보존했습니다."
    echo "$task_runtime_reason"
    echo "CUDA 환경 보존 위치: $task_preserved_cuda_attempt"
  fi
  echo "CPU 모델 패키지를 설치합니다. 진행 기록: $task_root/model-cpu-install.log"
  if ! task_install_model_environment cpu 2>&1 | tee "$task_root/model-cpu-install.log"; then
    echo "CPU 모델 설치에 실패했습니다. $task_root/model-cpu-install.log에서 원인을 확인하세요. 설치 자료는 보존했으며, 네트워크·여유 공간을 확인한 뒤 새 폴더에 설치하세요." >&2
    return 1
  fi
}

task_write_prepare_catalog() {
  printf '#!/usr/bin/env bash\nexec %q build --directory %q --profile %q --library %q --device %q --app-root %q --accept-source-terms\n' \
    "$task_root/app/software/bin/skinpipe-catalog" "$task_root/catalog" "$task_root/model-profile/profile.json" "$task_library" "$task_device" "$task_root/app/state" > "$task_root/Prepare-Catalog.sh"
}

task_write_app_launchers() {
  # Keep top-level launchers bound to the app launchers when the app is updated.
  printf '#!/usr/bin/env bash\nexec %q "$@"\n' "$task_root/app/Start-SkinPipe.sh" > "$task_root/Start-SkinPipe.sh"
  printf '#!/usr/bin/env bash\nexec %q "$@"\n' "$task_root/app/Stop-SkinPipe.sh" > "$task_root/Stop-SkinPipe.sh"
}

task_install_main() {
  if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo '사용법: bash install-from-sources-ubuntu.sh SkinPipe.whl 새설치폴더 [skin-starter|full] [auto|cpu|cuda]' >&2
    return 2
  fi
  task_wheel="$(realpath "$1")"
  task_root="$(task_new_install_directory "$2")" || return
  task_library="${3:-skin-starter}"
  task_requested_device="${4:-auto}"
  task_scripts="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  source /etc/os-release
  if [[ "${ID:-}" != ubuntu || "${VERSION_ID:-}" != 24.04 || "$(uname -m)" != x86_64 || "$EUID" == 0 ]]; then
    echo 'Ubuntu 24.04 x86_64의 일반 사용자로 실행하세요. sudo는 패키지 설치 시만 요청합니다.' >&2
    return 2
  fi
  if [[ "$task_library" != skin-starter && "$task_library" != full ]]; then
    echo 'skin-starter 또는 full을 선택하세요.' >&2
    return 2
  fi
  if [[ -e "$task_root" || -L "$task_root" ]]; then
    echo '기존 폴더는 보존합니다. 신규 설치에는 새 폴더를, 중단한 DB 구축에는 Prepare-Catalog.sh를 사용하세요.' >&2
    return 2
  fi
  if [[ ! -f "$task_wheel" || "$task_wheel" != *.whl ]]; then
    echo 'SkinPipe wheel 파일이 필요합니다.' >&2
    return 2
  fi
  task_select_model_runtime "$task_requested_device" || return
  cat <<'NOTICE'
공식 제공처에서 각 설치자에게 데이터를 직접 내려받아 로컬에서 ADMET을 계산합니다.
DB·예측 파일은 재배포하지 않습니다. 각 출처 이용조건은 docs/local-catalog-build.ko.md를 확인하세요.
skin-starter: 작은 피부 연구용 DB / full: 15개 출처의 대규모 구조 목록, 최소 300 GiB 여유 공간 권장.
auto: 사용 가능한 NVIDIA GPU를 먼저 확인하며, 없거나 설치·계산 확인이 실패하면 CPU로 진행합니다.
전체 CPU 계산은 컴퓨터에 따라 수일 이상 걸릴 수 있습니다. 중단 후 같은 폴더에서 재개할 수 있습니다.
NOTICE
  echo "$task_runtime_reason"
  umask 077
  mkdir -p "$task_root"
  bash "$task_scripts/install-research-ubuntu.sh" "$task_wheel" "$task_root/app"
  task_python="$task_root/app/software/bin/python"
  "$task_python" -m pip install "${task_wheel}[catalog-build,catalog-index]"
  # Separate immutable environments: a failed CUDA attempt never becomes a CPU
  # environment through in-place dependency replacement.
  task_model_python="$task_root/model-runtime/bin/python"
  task_assets="$("$task_python" -c 'from skinpipe.platform.catalog_bootstrap import ASSETS; print(ASSETS)')"
  task_prepare_model_runtime "$task_requested_device" || return
  "$task_model_python" -I -m skinpipe.platform.catalog_model "$task_root/model-profile"
  "$task_model_python" -I -m pip freeze > "$task_root/model-environment.lock.txt"
  "$task_python" - "$task_root/model-runtime-selection.json" "$task_requested_device" "$task_runtime" "$task_device" "$task_runtime_reason" "$task_preserved_cuda_attempt" <<'PY'
import datetime
import json
from pathlib import Path
import sys

path, requested, runtime, policy, reason, preserved = sys.argv[1:]
Path(path).write_text(json.dumps({
    'requested_device': requested, 'installed_runtime': runtime,
    'catalog_device_policy': policy, 'reason': reason,
    'preserved_cuda_attempt': preserved or None,
    'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
}, ensure_ascii=False, indent=2) + '\n')
PY
  task_write_prepare_catalog
  printf '#!/usr/bin/env bash\nexec %q status --directory %q\n' "$task_root/app/software/bin/skinpipe-catalog" "$task_root/catalog" > "$task_root/Catalog-Status.sh"
  task_write_app_launchers
  chmod 700 "$task_root/Prepare-Catalog.sh" "$task_root/Catalog-Status.sh" "$task_root/Start-SkinPipe.sh" "$task_root/Stop-SkinPipe.sh"
  "$task_root/Prepare-Catalog.sh"
  echo "DB 준비·등록 완료. $task_root/Start-SkinPipe.sh 실행 후 표시되는 주소를 여세요."
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  task_install_main "$@"
fi
