"""Owned single-researcher installation; no legacy database or job recovery.

The HTTP service deliberately binds loopback. Remote browsers use an SSH tunnel
to the same port, including the optional OAuth callback relay.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager, nullcontext
import fcntl
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import hmac
import importlib.util
import json
import mimetypes
import os
from pathlib import Path
import secrets
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import threading
from urllib.parse import urlsplit, unquote
from urllib.request import Request, build_opener, ProxyHandler, HTTPRedirectHandler

from .artifacts import atomic_bytes, canonical_json
from .cli import _outside_git, read_config, service_from_config
from .export_assets import export_asset_root

SCHEMA = "skinpipe-personal-install-v1"
OWNER = "researcher"


def _database_wait_notice():
    print("화합물 DB의 원본·인덱스 무결성을 확인하고 있습니다. 대형 DB는 수분이 걸립니다. "
          "완료 또는 오류가 표시될 때까지 기다리고 등록 명령을 반복하지 마세요.", file=sys.stderr, flush=True)


def _runtime_state(root):
    with (root / "server.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return "running" if (root / "server.json").exists() else "preparing_or_cleanup"
        if (root / "server.json").exists() or (root / "postgres/postmaster.pid").exists():
            return "cleanup_confirmation_required"
        return "stopped"


@contextmanager
def _configuration_change(root):
    # The same lifetime lock covers cold startup, registration and shutdown;
    # a missing HTTP marker alone never authorizes a configuration mutation.
    with (root / "server.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("서버 준비·실행·정리 또는 설정 변경 중입니다. stop 요청만으로 설정을 바꾸지 말고 "
                             "원래 Start 터미널 종료와 status=stopped를 확인하세요.") from error
        if (root / "server.json").exists() or (root / "postgres/postmaster.pid").exists():
            raise ValueError("DB·모델 설정을 바꾸기 전에 stop을 실행하고 원래 Start 터미널 종료와 "
                             "status=stopped를 확인하세요.")
        yield installation(root)[1]


def _pg_binary(name):
    executable = shutil.which(name)
    if executable:
        return executable
    candidates = sorted(Path("/usr/lib/postgresql").glob("*/bin/" + name), reverse=True)
    candidates += sorted(Path("/opt/homebrew/opt").glob("postgresql*/bin/" + name), reverse=True)
    if not candidates:
        raise ValueError("PostgreSQL 실행 파일이 없습니다. 설치 안내의 PostgreSQL 설치를 먼저 실행하세요.")
    return str(candidates[0])


def _pg(root, command):
    return subprocess.run([_pg_binary("pg_ctl"), "-D", str(root / "postgres"), *command],
                          check=True, capture_output=True, text=True)


def _start_database(root, config):
    if not (root / "postgres/postmaster.pid").exists():
        _pg(root, ["-l", str(root / "postgres.log"), "-w", "start"])
    import psycopg
    with psycopg.connect(**config["database"]) as conn:
        conn.execute("SELECT 1")


def setup(root, *, database_package=None, scout_config=None, vina=None):
    """Create a new cluster, never adopt an existing PostgreSQL ledger."""
    if database_package and scout_config:
        raise ValueError("DB 패키지와 기존 DB 등록 중 하나를 선택하세요.")
    root = _outside_git(root)
    if root.exists() and any(root.iterdir()):
        raise ValueError("설치 폴더가 비어 있지 않습니다. 새 폴더를 선택하세요.")
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    root.chmod(0o700)
    sock = Path(tempfile.mkdtemp(prefix="skinpipe-pg-"))
    sock.chmod(0o700)
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    log = root / "install.log"
    with log.open("w") as output:
        subprocess.run([_pg_binary("initdb"), "-D", str(root / "postgres"),
                        "--auth-local=trust", "--auth-host=reject", "--no-locale", "--encoding=UTF8"],
                       check=True, stdout=output, stderr=subprocess.STDOUT)
    # Only this owner's private Unix socket is reachable; TCP stays disabled.
    with (root / "postgres/postgresql.conf").open("a") as stream:
        stream.write("\nlisten_addresses = ''\nport = " + str(port)
                     + "\nunix_socket_directories = '" + str(sock).replace("'", "''")
                     + "'\nunix_socket_permissions = 0700\n")
    config = {"schema_version": 1, "root": str(root / "research-data"),
              "database": {"host": str(sock), "port": port, "dbname": "skinpipe_research"},
              "worker_tokens": {}, "legacy": {},
              "discovery": {"literature": {"enabled": True, "max_papers": 2, "timeout": 20, "review_timeout": 300,
                  "review_provider": "chatgpt", "model": "gpt-6-luna", "reasoning_effort": "max",
                  "owner_id": OWNER, "deployment_mode": "local_single_user",
                  "state_dir": str(root / "chatgpt")}}}
    if database_package or scout_config:
        from .database_bundle import install_database, register_database
        _database_wait_notice()
        if database_package:
            receipt = install_database(database_package, root / "databases")
        else:
            receipt = register_database(json.loads(Path(scout_config).read_text()), root / "databases")
        config["scout"] = receipt["config"]
    binary = vina or shutil.which("vina")
    if binary:
        config["discovery"]["research_docking"] = {
            "engine": "vina", "engine_paths": {"vina": str(Path(binary).resolve())}}
    atomic_bytes(root / "config.json", canonical_json(config))
    atomic_bytes(root / "personal-install.json", canonical_json({
        "schema": SCHEMA, "installation_id": secrets.token_hex(16), "root": str(root),
        "owner_id": OWNER, "os_uid": os.getuid(),
        "database": config["database"], "legacy_recovery": False}))
    try:
        _pg(root, ["-l", str(root / "postgres.log"), "-w", "start"])
        import psycopg
        with psycopg.connect(**{**config["database"], "dbname": "postgres"}, autocommit=True) as conn:
            conn.execute("CREATE DATABASE skinpipe_research")
        service_from_config(config)
    finally:
        if (root / "postgres/postmaster.pid").exists():
            _pg(root, ["-m", "fast", "-w", "stop"])
    return {"status": "installed", "root": str(root), "catalog": bool(config.get("scout"))}


def installation(root):
    root = _outside_git(root)
    marker = json.loads((root / "personal-install.json").read_text())
    config = read_config(root / "config.json")
    if (marker.get("schema") != SCHEMA or marker.get("root") != str(root)
            or marker.get("owner_id") != OWNER or marker.get("os_uid") != os.getuid()
            or marker.get("database") != config["database"] or config.get("legacy")
            or Path(config["root"]) != root / "research-data"):
        raise ValueError("개인 설치 기록과 실행 설정이 일치하지 않습니다. 기존 분석 DB는 자동으로 열지 않습니다.")
    return root, config


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


class PersonalServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, service, config, *, port=4317, asset_root=None):
        self.service, self.config = service, config
        self.asset_root = Path(asset_root or export_asset_root()) / "public"
        self.session = secrets.token_urlsafe(32)
        self.signin = {"status": "idle"}
        self.signin_lock = threading.Lock()
        self.cancel_signin = threading.Event()
        self.callback_base = None
        self.stopping = threading.Event()
        self.coordinator_error = None
        self.export_coordinator_error = None
        from .research_toxicity import toxicity_profile_status
        self.additional_toxicity = toxicity_profile_status(config.get("discovery", {}).get("toxicity_profile"))
        super().__init__(("127.0.0.1", port), PersonalHandler)

    def connection_status(self):
        status = self.service.personal_literature.status(owner=OWNER)
        status.update(connect_path="/v2/setup/chatgpt/connect", signin_status=self.signin["status"],
                      login_in_progress=self.signin["status"] == "pending")
        if self.signin.get("error"):
            status["error"] = self.signin["error"]
        return status

    def connect(self):
        with self.signin_lock:
            if self.signin["status"] == "pending":
                return dict(self.signin)
            self.signin = {"status": "pending"}
            self.cancel_signin.clear()
            ready = threading.Event()

            def authorize(url):
                self.signin["authorization_url"] = url
                ready.set()

            def work():
                from ..literature.chatgpt import ChatGPTClient
                try:
                    value = self.config["discovery"]["literature"]
                    with ChatGPTClient(value["state_dir"], model=value["model"],
                          authorization_ready=authorize,
                          callback_redirect=f"http://127.0.0.1:{self.server_port}/auth/callback",
                          callback_listener_ready=lambda url: setattr(self, "callback_base", url)) as client:
                        client.sign_in(open_browser=False, cancel_event=self.cancel_signin)
                    self.signin = {"status": "connected"}
                except Exception as error:
                    self.signin = {"status": "failed", "error": getattr(error, "code", type(error).__name__)}
                finally:
                    self.callback_base = None
                    ready.set()
            threading.Thread(target=work, daemon=True, name="optional-chatgpt-signin").start()
        ready.wait(10)
        return dict(self.signin)

    def coordinate(self, *, export_only=False):
        while not self.stopping.is_set():
            try:
                worked = self.service.discovery.run_once(shutdown_event=self.stopping, export_only=export_only)
                if not export_only and not self.stopping.is_set():
                    self.service.discovery.reconcile()
                setattr(self, "export_coordinator_error" if export_only else "coordinator_error", None)
            except Exception as error:
                setattr(self, "export_coordinator_error" if export_only else "coordinator_error", type(error).__name__)
                worked = False
            self.stopping.wait(0.2 if worked else 1)


class PersonalHandler(BaseHTTPRequestHandler):
    def log_message(self, *_args):
        pass  # No user input, authorization URL, callback code or credential logs.

    def _headers(self, status, content_type, length, *, filename=None, cookie=False):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(length))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; worker-src 'self' blob:; frame-ancestors 'none'; object-src 'none'; base-uri 'none'")
        if cookie:
            self.send_header("Set-Cookie", "skinpipe_local=" + self.server.session + "; HttpOnly; SameSite=Strict; Path=/")
        if filename:
            self.send_header("Content-Disposition", 'attachment; filename="' + filename.replace('"', '') + '"')
        self.end_headers()

    def _json(self, value, status=200):
        data = canonical_json(value)
        self._headers(status, "application/json; charset=utf-8", len(data))
        self.wfile.write(data)

    def _host(self):
        return self.headers.get("Host") in {f"127.0.0.1:{self.server.server_port}", f"localhost:{self.server.server_port}"}

    def _authenticated(self):
        cookie = SimpleCookie()
        try:
            cookie.load(self.headers.get("Cookie", ""))
            value = cookie.get("skinpipe_local")
            return value is not None and hmac.compare_digest(value.value, self.server.session)
        except Exception:
            return False

    def _guard(self, *, mutation=False):
        if not self._host() or not self._authenticated():
            raise PermissionError("이 설치의 연구 화면에서 요청해 주세요.")
        if self.headers.get("Sec-Fetch-Site") in {"cross-site", "same-site"}:
            raise PermissionError("다른 사이트에서 보낸 요청은 허용되지 않습니다.")
        if mutation:
            origin = self.headers.get("Origin")
            allowed = {f"http://127.0.0.1:{self.server.server_port}", f"http://localhost:{self.server.server_port}"}
            if origin is not None and origin not in allowed:
                raise PermissionError("연구 화면의 출처가 일치하지 않습니다.")
            if self.headers.get("X-SkinPipe-Local") != "1":
                raise PermissionError("연구 화면의 요청 확인 값이 없습니다.")

    def do_GET(self):
        try:
            if not self._host():
                raise PermissionError("허용되지 않은 호스트입니다.")
            path = urlsplit(self.path).path
            if path == "/auth/callback":
                if not self.server.callback_base or self.server.signin["status"] != "pending":
                    raise ValueError("진행 중인 연결 요청이 없습니다.")
                opener = build_opener(ProxyHandler({}), _NoRedirect())
                with opener.open(Request(self.server.callback_base + self.path), timeout=5) as response:
                    data = response.read(16384)
                self._headers(200, "text/plain; charset=utf-8", len(data))
                self.wfile.write(data)
                return
            if path.startswith("/v2/"):
                self._guard()
                if path == "/v2/setup/status":
                    catalog = self.server.service.get("/v2/catalog-status", owner=OWNER)
                    return self._json({"status": "ready", "database": "PostgreSQL 연결됨", "legacy_available": False,
                        "catalog": catalog, "docking": self.server.config["discovery"].get("research_docking", {"status": "not_configured"}),
                        "chatgpt": self.server.connection_status(), "coordinator_error": self.server.coordinator_error,
                        "export_coordinator_error": self.server.export_coordinator_error,
                        "next_actions": [] if catalog.get("configured", True) else ["로컬 DB를 설치하면 전체 후보 탐색과 기존 독성 예측을 사용할 수 있습니다."]})
                if path == "/v2/me":
                    return self._json({"owner_id": OWNER, "role": "researcher"})
                if path == "/v2/local/profile":
                    return self._json({"owner_id": OWNER, "additional_toxicity": self.server.additional_toxicity})
                parts = path.strip("/").split("/")
                if len(parts) == 4 and parts[:2] in (["v2", "exports"], ["v2", "result-exports"]) and parts[3] == "download":
                    file, mime, name = self.server.service.exports.download(parts[2], owner=OWNER)
                    self._headers(200, mime, file.stat().st_size, filename=name)
                    with file.open("rb") as stream:
                        shutil.copyfileobj(stream, self.wfile, 1024 * 1024)
                    return
                return self._json(self.server.service.get(self.path, owner=OWNER))
            requested = "local-research.html" if path == "/" else "index.html" if path == "/legacy.html" else unquote(path).lstrip("/")
            file = (self.server.asset_root / requested).resolve()
            if not file.is_relative_to(self.server.asset_root.resolve()) or not file.is_file():
                raise KeyError("화면 파일을 찾을 수 없습니다.")
            self._headers(200, mimetypes.guess_type(file.name)[0] or "application/octet-stream", file.stat().st_size,
                          cookie=requested in {"local-research.html", "index.html"})
            with file.open("rb") as stream:
                shutil.copyfileobj(stream, self.wfile, 1024 * 1024)
        except Exception as error:
            self._error(error)

    def do_POST(self):
        try:
            self._guard(mutation=True)
            if self.headers.get_content_type() != "application/json":
                raise ValueError("JSON 입력이 필요합니다.")
            length = int(self.headers.get("Content-Length", "0"))
            limit = 72 * 1024**2 if self.path == "/v2/research-inputs/uploads" else 18 * 1024**2
            if not 0 < length <= limit:
                raise ValueError("입력 크기가 허용 범위를 초과했습니다.")
            body = json.loads(self.rfile.read(length))
            if self.path == "/v2/setup/shutdown":
                if body != {}:
                    raise ValueError("종료 요청 형식이 잘못되었습니다.")
                self._json({"status": "stop_requested"})
                threading.Thread(target=self.server.shutdown, daemon=True).start()
                return
            if self.path.startswith("/v2/setup/chatgpt/"):
                if body != {}:
                    raise ValueError("연결 요청에 비밀 값이나 추가 입력은 필요하지 않습니다.")
                action = self.path.rsplit("/", 1)[-1]
                if action == "connect":
                    return self._json(self.server.connect())
                self.server.cancel_signin.set()
                if action == "logout":
                    with self.server.service.personal_literature._client() as client:
                        client.disconnect()
                elif action != "cancel":
                    raise KeyError("알 수 없는 연결 요청입니다.")
                return self._json({"status": "cancel_requested" if action == "cancel" else "disconnected"})
            self._json(self.server.service.post(self.path, body, owner=OWNER))
        except Exception as error:
            self._error(error)

    def _error(self, error):
        status = 403 if isinstance(error, PermissionError) else 404 if isinstance(error, KeyError) else 400 if isinstance(error, (ValueError, TypeError)) else 500
        message = str(error) if status != 500 else "처리를 완료하지 못했습니다. 설치 상태를 확인한 뒤 다시 시도해 주세요."
        self._json({"error": message, "error_type": type(error).__name__}, status)


def serve(root, *, port=4317):
    root, config = installation(root)
    lock = (root / "server.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        lock.close()
        raise ValueError("이 설치의 서버가 이미 시작 준비 중이거나 실행 중입니다. "
                         "기존 Start 터미널과 접속 주소를 확인하고 시작 명령을 반복하지 마세요.") from error
    except BaseException:
        lock.close()
        raise
    service = server = None
    threads = []
    database_attempted = False
    stopping_requested = False

    def stop(_signal, _frame):
        nonlocal stopping_requested
        if not stopping_requested:
            stopping_requested = True
            raise KeyboardInterrupt

    previous_handlers = {s: signal.signal(s, stop) for s in (signal.SIGTERM, signal.SIGINT)}
    try:
        root, config = installation(root)
        # The lifetime lock proves a previous HTTP record is no longer active.
        (root / "server.json").unlink(missing_ok=True)
        print("시작 준비 중입니다. 설치한 화합물 DB가 크면 원본·인덱스 확인에 수분이 걸립니다. "
              "접속 주소가 표시될 때까지 기다리고 새 시작 명령을 반복하지 마세요.", file=sys.stderr, flush=True)
        # A reboot may remove the private temporary socket directory.
        Path(config["database"]["host"]).mkdir(mode=0o700, parents=True, exist_ok=True)
        database_attempted = True
        _start_database(root, config)
        service = service_from_config(config)
        server = PersonalServer(service, config, port=port)
        for name, export_only in [("personal-research-coordinator", False), ("personal-export-coordinator", True)]:
            thread = threading.Thread(target=server.coordinate, kwargs={"export_only": export_only}, daemon=True, name=name)
            threads.append(thread)
            thread.start()
        atomic_bytes(root / "server.json", canonical_json({"pid": os.getpid(), "port": server.server_port,
                                                           "local_session": server.session}))
        print(f"SkinPipe: http://127.0.0.1:{server.server_port}", flush=True)
        server.serve_forever(poll_interval=0.2)
    except KeyboardInterrupt:
        pass
    finally:
        stopping_requested = True
        print("종료·자원 정리 중입니다. 원래 시작 터미널이 종료될 때까지 기다려 주세요.",
              file=sys.stderr, flush=True)
        try:
            try:
                try:
                    (root / "server.json").unlink(missing_ok=True)
                finally:
                    if server is not None:
                        server.stopping.set()
                        server.cancel_signin.set()
                        try:
                            server.server_close()
                        finally:
                            # Native loops finish their owned process groups before
                            # the private ledger can be closed.
                            for thread in threads:
                                if thread.ident is not None:
                                    thread.join()
            finally:
                if service is not None:
                    service.close_catalog_sources()
        finally:
            try:
                try:
                    (root / "server.json").unlink(missing_ok=True)
                finally:
                    if database_attempted and (root / "postgres/postmaster.pid").exists():
                        _pg(root, ["-m", "fast", "-w", "stop"])
            finally:
                lock.close()
                for sig, handler in previous_handlers.items():
                    signal.signal(sig, handler)


def main(argv=None):
    parser = argparse.ArgumentParser(description="SkinPipe 개인 연구 설치 및 실행")
    parser.add_argument("command", choices=["setup", "serve", "status", "stop", "database", "toxicity"])
    parser.add_argument("toxicity_action", nargs="?", choices=["register", "status", "disable"])
    parser.add_argument("--root", default="~/.local/share/skinpipe-research")
    parser.add_argument("--database-package")
    parser.add_argument("--scout-config")
    parser.add_argument("--vina")
    parser.add_argument("--profile", help="설치 담당자가 검증한 추가 독성 모델 프로필 JSON")
    parser.add_argument("--port", type=int, default=4317)
    args = parser.parse_args(argv)
    if args.command != "toxicity" and (args.toxicity_action or args.profile):
        parser.error("모델 프로필은 toxicity 명령에서만 지정할 수 있습니다.")
    if args.command == "setup":
        result = setup(args.root, database_package=args.database_package, scout_config=args.scout_config, vina=args.vina)
    elif args.command == "serve":
        return serve(args.root, port=args.port)
    else:
        root, config = installation(args.root)
        if args.command == "stop":
            if not (root / "server.json").exists():
                raise ValueError("접속 주소가 아직 없거나 서버가 종료 중·중지 상태입니다. status를 확인하세요. "
                                 "준비 중 중단은 원래 Start 터미널에서 Ctrl+C를 누르고 종료를 기다립니다.")
            record = json.loads((root / "server.json").read_text())
            req = Request(f"http://127.0.0.1:{int(record['port'])}/v2/setup/shutdown", data=b"{}",
                          headers={"Cookie": "skinpipe_local=" + record["local_session"],
                                   "X-SkinPipe-Local": "1", "Content-Type": "application/json"})
            with build_opener(ProxyHandler({}), _NoRedirect()).open(req, timeout=10) as response:
                result = json.load(response)
        elif args.command == "toxicity":
            from .research_toxicity import toxicity_profile_status, validate_toxicity_profile
            action = args.toxicity_action or "status"
            with (_configuration_change(root) if action != "status" else nullcontext(config)) as config:
                if action == "register":
                    if not args.profile:
                        raise ValueError("검증한 로컬 모델 프로필을 --profile로 지정하세요.")
                    source = Path(args.profile).expanduser().resolve()
                    validate_toxicity_profile(source)
                    private = root / "optional-toxicity"
                    private.mkdir(mode=0o700, exist_ok=True)
                    private.chmod(0o700)
                    target = private / "profile.json"
                    atomic_bytes(target, source.read_bytes())
                    target.chmod(0o600)
                    config.setdefault("discovery", {})["toxicity_profile"] = str(target)
                    atomic_bytes(root / "config.json", canonical_json(config))
                elif action == "disable":
                    config.setdefault("discovery", {}).pop("toxicity_profile", None)
                    atomic_bytes(root / "config.json", canonical_json(config))
            result = toxicity_profile_status(config.get("discovery", {}).get("toxicity_profile"))
        elif args.command == "database":
            from .database_bundle import install_database, register_database
            with _configuration_change(root) as config:
                if bool(args.database_package) == bool(args.scout_config):
                    raise ValueError("DB 패키지 또는 기존 DB 설정을 하나 선택하세요.")
                _database_wait_notice()
                receipt = (install_database(args.database_package, root / "databases") if args.database_package
                           else register_database(json.loads(Path(args.scout_config).read_text()), root / "databases"))
                config["scout"] = receipt["config"]
                atomic_bytes(root / "config.json", canonical_json(config))
                result = {"status": "database_selected", "catalog": True}
        else:
            result = {"status": _runtime_state(root),
                      "catalog": bool(config.get("scout")), "python": __import__("sys").version.split()[0],
                      "postgres": bool((root / "postgres/PG_VERSION").exists()),
                      "chatgpt_dependency": importlib.util.find_spec("jwt") is not None}
            if result["status"] == "preparing_or_cleanup":
                result["message"] = "서버 준비·설정 변경·자원 정리 중입니다. 원래 터미널의 안내와 종료를 확인하세요."
            elif result["status"] == "cleanup_confirmation_required":
                result["message"] = "이전 서버 또는 PostgreSQL 기록의 종료 확인이 필요합니다. 설정을 바꾸지 말고 설치 담당자에게 확인하세요."
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
