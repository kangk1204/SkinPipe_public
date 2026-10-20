# 설치·DB·모델 준비

[English](catalog.en.md) · [README](../README.ko.md)

## 파일은 어디에 저장되나요?

기본 설치 폴더는 Git 밖의 `~/SkinPipe`입니다. 출처를 알 수 없는 기존 폴더와 심볼릭 링크 경로는 덮어쓰지 않습니다.

| 위치 | 내용 |
|---|---|
| `~/SkinPipe/installation.json` | 소스 리비전·wheel 해시·설치 선택·준비 완료 상태 |
| `~/SkinPipe/source` | 설치에 사용한 공개 소스 |
| `~/SkinPipe/runtime/app/software` | 앱 Python 환경 |
| `~/SkinPipe/runtime/model-runtime` | 별도 ADMET 모델 환경 |
| `~/SkinPipe/runtime/model-profile` | 모델 프로필·동일성 정보 |
| `~/SkinPipe/runtime/catalog` | 공식 다운로드·검증 구조·예측·감사·검색 인덱스 |
| `~/SkinPipe/runtime/app/state` | 비공개 설정·PostgreSQL·연구 결과 |
| `~/SkinPipe/*.log`, `~/SkinPipe/runtime/*.log` | 설치·모델 환경 로그 |

DB·모델은 설치 중 공식 출처에서 받습니다. 위 실행 자료를 Git에 넣거나 공개 release로 올리지 마세요.

## 작은 DB와 전체 라이브러리

기본 **skin-starter**는 버전이 고정된 **huskinDB·SkinPiX**로 구성하며 사용법 학습과 작은 검증에 적합합니다. 전체 계획은 여기에 ChEMBL·ChEBI·COCONUT·LOTUS·NPAtlas·NTP ICE·ZINC·Mcule을 포함합니다. 정확한 URL·해시·명시된 이용조건은 [sources.json](../skinpipe/skinpipe/platform/_catalog_assets/sources.json)과 [ZINC 파일 인덱스](../skinpipe/skinpipe/platform/_catalog_assets/zinc-artifacts.json)에 있습니다.

설치자마다 각 제공처에서 직접 받습니다. 출처별 이용조건을 읽고 사용하세요. 일부 출처는 재배포를 제한하며, 설치했다고 통합 DB를 다른 사람에게 전달할 권한이 생기지는 않습니다. 로컬 데이터의 출처 표기를 유지하세요. SkinPipe가 외부 데이터·모델의 이용권을 부여하지 않습니다.

전체 라이브러리를 고르려면 설치 스크립트를 받고 내용·이용조건을 검토한 뒤 명시적으로 실행합니다.

```sh
python3 -c 'import urllib.request; urllib.request.urlretrieve("https://raw.githubusercontent.com/kangk1204/SkinPipe_public/main/install.py", "SkinPipe-install.py")'
python3 SkinPipe-install.py --prefix ~/SkinPipe-full --library full --device auto --accept-source-terms
```

**300 GiB 이상**의 여유 공간과 유지되는 터미널 세션이 필요합니다. CPU 전수 ADMET은 수일 이상 걸릴 수 있으며 자료 크기·장비에 따라 달라집니다. 기본 starter는 **20 GiB**를 검사합니다. 다른 사용자의 기존 연구·DB를 복사하지 않습니다.

## 자동으로 수행하는 작업

1. 공개 GitHub의 특정 리비전을 고정·저장하고 wheel을 빌드합니다.
2. PostgreSQL·CPU Vina 등 시스템 패키지와 비공개 앱 환경을 설치합니다.
3. 별도 환경에 버전 고정 ADMET-AI와 모델을 받고 동일성을 검사합니다.
4. 출처 스냅샷 다운로드·해시 확인·구조 정리·ADMET 계산을 합니다. 미지원·실패를 구분해 남깁니다.
5. 예측 감사를 거쳐 전역 검색 인덱스를 만들고 DB·모델을 앱에 등록합니다.
6. 로컬 브라우저 서버를 시작합니다. 준비가 끝나면 표시된 주소를 엽니다.

**auto**는 NVIDIA 드라이버와 실제 CUDA 계산을 검사합니다. 사용할 수 없으면 별도 CPU 환경으로 전환하고 실패 환경·로그를 보존합니다. **cuda**를 명시했는데 실패하면 이유를 표시하며 GPU 사용으로 처리하지 않습니다. **cpu**는 GPU 검사를 생략합니다. 드라이버 교체나 다른 사용자의 프로세스 종료는 하지 않습니다.

전체 DB의 ADMET과 연구별 추가 독성 모델 기록은 구분합니다. 구조·모델·전처리 조건이 맞는 예측만 재사용합니다. 새 구조의 미지원 값은 미평가로 남기며 이름·유사도만으로 다른 물질의 독성을 가져오지 않습니다.

## 진행 확인과 재개

설치 터미널에 단계와 로그 경로가 나옵니다. 다른 터미널에서:

```sh
bash ~/SkinPipe/runtime/Catalog-Status.sh
```

다운로드 → 구조 정리 → DB → ADMET → 감사 → 인덱스 상태를 확인합니다. `Prepare-Catalog.sh`가 생성된 뒤 중단됐다면 **동일한 prefix·library·device**의 원래 명령으로 재개합니다. 완료를 확인한 단계는 대조 후 재사용합니다. 설치기와 DB 구축에 각각 잠금이 있어 같은 폴더의 중복 구축을 막습니다.

그 전에 패키지 설치가 실패하면 표시된 로그를 읽습니다. 부분 환경은 진단용으로 남습니다. 원인을 해결한 뒤 새 `--prefix`를 사용하세요. 소스·모델·설정이 바뀌었을 때도 진행 중인 계획을 바꾸지 말고 새 폴더를 사용합니다.

준비 후에는 매번 `Start-SkinPipe.sh`로 시작합니다. 연구마다 전체 DB를 재구축할 필요가 없습니다. 브라우저 **설치 상태**에서 앱·DB·인덱스·모델·엔진을 확인하고, 세부 장치 선택은 runtime의 `model-runtime-selection.json`과 로그에서 봅니다.

## macOS·Windows의 브라우저에서 접속

Ubuntu에서 앱 실행 터미널을 열어 둡니다. Mac/Windows 터미널에서는 아래 계정·주소를 본인 서버로 바꿔 실행합니다.

```sh
ssh -N -L 4317:127.0.0.1:4317 YOUR_USER@YOUR_UBUNTU_SERVER
```

터널을 유지하고 Mac/Windows 브라우저에서 **http://127.0.0.1:4317** 을 엽니다. 로컬 포트도 **4317**을 유지하세요. 앱이 호스트와 포트를 검사하므로 로컬 4318 → 서버 4317 연결은 거부됩니다. “Address already in use”가 나오면 Ubuntu의 브라우저를 사용하거나, 4317을 사용하는 앱을 먼저 확인한 뒤 해당 포트를 비웁니다. 서버 파일 경로를 브라우저 주소에 넣지 마세요. 터널 문제를 피하려고 인증을 끄거나 HTTP 포트를 공개하지 마세요.

## 종료와 업데이트

해당 설치의 `Stop-SkinPipe.sh`로 종료하고 자원 정리 완료를 기다린 뒤 상태 폴더를 복사·이동합니다. 백업할 때 DB·모델 프로필·연구 상태를 함께 보관합니다. 앱 실행 중 PostgreSQL 폴더를 삭제하지 않습니다.

설치기는 소스를 특정 리비전으로 고정합니다. 같은 명령은 **그 설치를 재개·시작**하며 최신 main으로 몰래 바꾸지 않습니다. 추후 버전은 새 prefix와 해당 버전의 이전 안내를 사용하고 이전 연구를 보존하세요. 실행 중인 환경에 wheel을 덮어쓰는 방식으로 업데이트하지 않습니다.

## 자주 생기는 문제

| 증상 | 다음 행동 |
|---|---|
| OS·root 오류 | Ubuntu 24.04 x86_64 일반 계정 사용 |
| 기존 폴더 오류 | 새 prefix 선택. 다른 설치 삭제 금지 |
| 다운로드·해시 오류 | 네트워크·제공처를 확인하고 같은 구축 재시도 |
| 공간 부족 | 진행 중인 연구 밖에서 공간을 확보하거나 새 볼륨 선택 |
| 다른 구축 실행 중 | 기존 터미널·상태를 확인. 중복 구축하지 않음 |
| auto GPU 검사 실패 | CPU 전환 표시·CUDA 로그 확인 |
| 모델·인덱스 미설정 | 설치 상태에서 DB·모델 준비 완료 |
| 후보 0개 | 제외·보류·탐색 범위를 확인하고 별도 연구에서 조건 변경 |
| Excel 실패 | 오류를 확인하고 내보내기만 재요청. 연구는 보존 |

공개 기여 규칙은 [CONTRIBUTING.md](../CONTRIBUTING.md)에 있습니다. 연구 입력·결과·모델은 공개 issue에도 첨부하지 마세요.
