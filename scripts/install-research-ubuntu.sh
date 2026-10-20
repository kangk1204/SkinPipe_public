#!/usr/bin/env bash
set -euo pipefail
# One-time installation by the workstation owner. Never opens legacy jobs.
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo '사용법: bash install-research-ubuntu.sh SOFTWARE.whl 새설치폴더 [DATABASE.tar.gz]' >&2
  exit 2
fi
task_wheel="$(realpath "$1")"
task_root="$(realpath -ms "$2")"
task_database="${3:-}"
if [[ ! -f "$task_wheel" || "$task_wheel" != *.whl ]]; then
  echo '제공된 SkinPipe wheel 파일을 선택하세요.' >&2; exit 2
fi
if [[ -e "$task_root" || -L "$task_root" || "$(realpath -m "$2")" != "$task_root" ]]; then
  echo '기존 폴더를 덮어쓰지 않습니다. 새 설치 폴더를 선택하세요.' >&2; exit 2
fi
source /etc/os-release
if [[ "${ID:-}" != ubuntu || "${VERSION_ID:-}" != 24.04 ]]; then
  echo '이 설치기는 Ubuntu 24.04 전용입니다.' >&2; exit 2
fi
sudo apt-get update
sudo apt-get install -y python3.12 python3.12-venv postgresql-16 autodock-vina curl ca-certificates libxrender1 libxext6 libsm6 libgomp1
umask 077
mkdir -p "$task_root"
python3.12 -m venv "$task_root/software"
"$task_root/software/bin/python" -m pip install "${task_wheel}[research,catalog-index]"
task_options=(setup --root "$task_root/state" --vina /usr/bin/vina)
if [[ -n "$task_database" ]]; then
  task_options+=(--database-package "$(realpath "$task_database")")
fi
"$task_root/software/bin/skinpipe-research" "${task_options[@]}"
"$task_root/software/bin/python" -m pip freeze > "$task_root/environment.lock.txt"
sha256sum "$task_wheel" > "$task_root/software.sha256.txt"
printf '#!/usr/bin/env bash\nexec %q serve --root %q "$@"\n' "$task_root/software/bin/skinpipe-research" "$task_root/state" > "$task_root/Start-SkinPipe.sh"
printf '#!/usr/bin/env bash\nexec %q stop --root %q\n' "$task_root/software/bin/skinpipe-research" "$task_root/state" > "$task_root/Stop-SkinPipe.sh"
chmod 700 "$task_root/Start-SkinPipe.sh" "$task_root/Stop-SkinPipe.sh"
echo "설치 완료. $task_root/Start-SkinPipe.sh 실행 후 http://127.0.0.1:4317 을 여세요."
