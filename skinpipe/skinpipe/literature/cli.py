"""Local, single-user literature workspace; credentials never reach its page."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import math
import os
from pathlib import Path
import re
import secrets
import threading
from urllib.parse import parse_qs, urlsplit
import uuid
import webbrowser

from .chatgpt import ChatGPTError, REASONING_EFFORTS


def account_error(error):
    """Expose only stable diagnostics; provider prose and credentials stay private."""
    code = error.code if isinstance(error.code, str) and re.fullmatch(r"[a-z0-9_]{1,100}", error.code) else "account_unavailable"
    if getattr(error, "usage_exhausted", False) is True or code == "usage_limit_exceeded":
        code = "usage_exhausted"
    elif code == "rate_limited":
        code = "rate_limit_exceeded"
    messages = {
        "usage_exhausted": "이 앱에서 사용할 수 있는 ChatGPT 사용량 한도에 도달해 AI 검토를 보류했습니다. 한도가 회복된 뒤 직접 다시 실행해 주세요.",
        "rate_limit_exceeded": "ChatGPT 요청이 일시적으로 제한되어 AI 검토를 보류했습니다. 잠시 후 직접 다시 실행해 주세요.",
        "sign_in_cancelled": "로그인을 취소했습니다. Continue with ChatGPT로 다시 연결할 수 있습니다.",
        "sign_in_timed_out": "로그인 대기가 끝났습니다. Continue with ChatGPT로 다시 연결해 주세요.",
        "browser_unavailable": "로그인 브라우저를 열지 못했습니다. 기본 브라우저 설정을 확인하고 다시 연결해 주세요.",
        "access_denied": "로그인 또는 플랜 사용 권한이 허용되지 않았습니다. 필요하면 다시 연결해 주세요.",
    }
    value = {"code": code, "error": messages.get(code, "계정 요청을 완료하지 못했습니다. 진단 코드: " + code)}
    if code in {"usage_exhausted", "rate_limit_exceeded"}:
        delay = getattr(error, "retry_after_seconds", None)
        if isinstance(delay, (int, float)) and not isinstance(delay, bool) and math.isfinite(delay) and delay >= 0:
            value["retry_after_seconds"] = delay
        reset = getattr(error, "reset_at", None)
        if isinstance(reset, str) and re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z", reset):
            value["reset_at"] = reset
    return (429 if code in {"usage_exhausted", "rate_limit_exceeded"} else 400), value


def default_state_dir():
    return Path.home() / ".local" / "share" / "skinpipe" / "chatgpt"


def _private_root(path):
    """Never put accounts or literature exports into a source checkout."""
    path = Path(path).expanduser().absolute()
    for part in (path, *path.parents):
        if part.is_symlink():
            raise ValueError("Private runtime paths must not contain symlinks")
        if (part / ".git").exists():
            raise ValueError("Use private runtime storage outside Git")
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    if path.stat().st_uid != os.getuid() or path.stat().st_mode & 0o077:
        raise ValueError("Private runtime directory must be owned by you with mode 700")
    return path


def _target(query, *, kind="protein", taxon=9606, action="any"):
    from ..platform.targets import target_request
    # Resolving an exact target is a separate verification. A typed gene name
    # only scopes retrieval; it is never a verified molecular identity.
    return target_request({"query": query, "kind": kind, "taxon_id": taxon,
                           "desired_action": action})


def collect_local(query, output_dir, *, state_dir=None, model="gpt-6-luna", reasoning_effort="max", review=False,
                  max_papers=8, kind="protein", taxon=9606, action="any", review_timeout=300):
    from .collect import collect_evidence
    if reasoning_effort is not None and (not isinstance(reasoning_effort, str) or reasoning_effort not in REASONING_EFFORTS):
        raise ValueError("Unsupported literature reasoning_effort")
    if type(review_timeout) is not int or not 1 <= review_timeout <= 300:
        raise ValueError("Literature review_timeout must be within 1..300 seconds")
    reviewer = None
    client = None
    if review:
        from .chatgpt import ChatGPTClient
        if not model:
            raise ValueError("Choose an AI model before requesting literature review")
        client = ChatGPTClient(state_dir or default_state_dir(), model=model, reasoning_effort=reasoning_effort)
        def reviewer(packet):
            models = client.list_models()
            choices = models.get('models') if isinstance(models, dict) else None
            if not isinstance(choices, list) or not any(row.get('id') == model for row in choices if isinstance(row, dict)):
                raise ChatGPTError('selected_model_unavailable')
            return client.review(packet, timeout=review_timeout)
    try:
        return collect_evidence(_target(query, kind=kind, taxon=taxon, action=action),
                                Path(output_dir), reviewer=reviewer, max_papers=max_papers)
    finally:
        if client is not None:
            client.close()


class LiteratureServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address, *, state_dir, output_root, client_factory=None):
        if address[0] != "127.0.0.1":
            raise ValueError("This ChatGPT connection is a local single-user feature")
        state_path, output_path = (Path(p).expanduser().resolve() for p in (state_dir, output_root))
        if state_path.is_relative_to(output_path) or output_path.is_relative_to(state_path):
            raise ValueError("Account storage and literature reports must be disjoint directories")
        self.state_dir, self.output_root = _private_root(state_dir), _private_root(output_root)
        self.csrf = secrets.token_urlsafe(32)
        self.csp_nonce = secrets.token_urlsafe(32)
        self.operation = threading.Lock()
        self.login_guard = threading.Lock()
        self.login_cancel = None
        self.login_done = None
        if client_factory is None:
            from .chatgpt import ChatGPTClient
            client_factory = ChatGPTClient
        self.client_factory = client_factory
        super().__init__(address, LiteratureHandler)
        self.origin = "http://127.0.0.1:" + str(self.server_port)

    @contextmanager
    def client(self, model=None, reasoning_effort=None):
        options = {"model": model}
        if reasoning_effort is not None:
            options["reasoning_effort"] = reasoning_effort
        client = self.client_factory(self.state_dir, **options)
        try:
            yield client
        finally:
            close = getattr(client, "close", None)
            if close:
                close()


class LiteratureHandler(BaseHTTPRequestHandler):
    server_version = "SkinPipeLocalLiterature/1"

    def log_message(self, *args):
        pass  # No OAuth callbacks, query text or identities in access logs.

    def _reply(self, code, value, content_type="application/json; charset=utf-8"):
        raw = value if isinstance(value, bytes) else json.dumps(value, ensure_ascii=False, allow_nan=False).encode()
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(raw)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'none'; script-src 'nonce-" + self.server.csp_nonce + "'; style-src 'unsafe-inline'; connect-src 'self'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'")
        try:
            self.end_headers()
            self.wfile.write(raw)
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError):
            # Closing the page does not turn a cancelled login into a server
            # exception or keep the operation lock reserved.
            return

    def _origin_ok(self, *, mutation=False):
        expected = urlsplit(self.server.origin).netloc
        if self.headers.get("Host") != expected:
            return False
        origin = self.headers.get("Origin")
        if origin is not None and origin != self.server.origin:
            return False
        if self.headers.get("Sec-Fetch-Site") not in {None, "same-origin", "none"}:
            return False
        if mutation:
            return origin == self.server.origin and secrets.compare_digest(
                self.headers.get("X-SkinPipe-CSRF", ""), self.server.csrf)
        return True

    def _run_dir(self, query):
        values = parse_qs(query, strict_parsing=True)
        run_id = values.get("id", [""])[0]
        if not re.fullmatch(r"[0-9a-f]{32}", run_id):
            raise ValueError("Invalid literature run identifier")
        path = self.server.output_root / run_id
        if path.is_symlink() or not path.is_dir():
            raise ValueError("Unknown literature run")
        return path, values

    def do_GET(self):
        if not self._origin_ok():
            return self._reply(403, {"error": "Local same-origin access required"})
        parsed = urlsplit(self.path)
        try:
            if parsed.path == "/":
                from .local_page import PAGE
                return self._reply(200, PAGE.replace("__NONCE__", self.server.csrf).replace("__CSP_NONCE__", self.server.csp_nonce).encode(), "text/html; charset=utf-8")
            if parsed.path == "/api/status":
                with self.server.client() as client:
                    state = client.status()
                with self.server.login_guard:
                    state["login_in_progress"] = self.server.login_cancel is not None
                return self._reply(200, state)
            if parsed.path == "/api/usage":
                with self.server.client() as client:
                    state = client.usage()
                return self._reply(200, state)
            if parsed.path == "/api/models":
                with self.server.client() as client:
                    models = client.list_models()
                return self._reply(200, models if isinstance(models, dict) else {
                    "models": models, "note": "Catalog only; actual request completion establishes access"})
            if parsed.path == "/api/reports":
                reports = []
                paths = [p for p in self.server.output_root.iterdir()
                         if re.fullmatch(r"[0-9a-f]{32}", p.name) and not p.is_symlink() and p.is_dir()]
                for path in sorted(paths, key=lambda p: (p / "dossier.json").stat().st_mtime
                                   if (p / "dossier.json").is_file() else 0, reverse=True):
                    if re.fullmatch(r"[0-9a-f]{32}", path.name) and not path.is_symlink() and (path / "dossier.json").is_file():
                        reports.append({"id": path.name})
                return self._reply(200, {"reports": reports[:100]})
            if parsed.path in {"/api/report", "/api/source"}:
                directory, values = self._run_dir(parsed.query)
                relative = values.get("file", ["dossier.json"])[0] if parsed.path == "/api/source" else "dossier.json"
                if len(relative) > 256 or "\\" in relative or any(p in {"", ".", ".."} for p in relative.split("/")):
                    raise ValueError("Invalid source path")
                path = directory / relative
                if not path.resolve().is_relative_to(directory.resolve()) or any(p.is_symlink() for p in (path, *path.parents)):
                    raise ValueError("Source must remain inside the saved run")
                if not path.is_file() or path.stat().st_size > 16 * 1024 * 1024:
                    raise ValueError("Source missing or exceeds 16 MiB")
                if path.stat().st_nlink != 1:
                    raise ValueError("Linked source files cannot be served")
                if parsed.path == "/api/source":
                    dossier_path = directory / "dossier.json"
                    if dossier_path.is_symlink() or dossier_path.stat().st_size > 16 * 1024 * 1024:
                        raise ValueError("Invalid dossier source closure")
                    dossier = json.loads(dossier_path.read_bytes())
                    if not isinstance(dossier, dict) or dossier.get("schema") != "skinpipe-literature-evidence-v1":
                        raise ValueError("Unsupported dossier source closure")
                    if not isinstance(dossier.get("artifact_paths"), list):
                        raise ValueError("No bound source paths")
                    allowed = {record.get("path"): record for record in dossier.get("artifact_paths", [])
                               if isinstance(record, dict)}
                    if relative not in allowed:
                        raise ValueError("Source is not bound to this literature dossier")
                raw = path.read_bytes()
                if parsed.path == "/api/source" and (
                    hashlib.sha256(raw).hexdigest() != allowed[relative].get("sha256")
                    or len(raw) != allowed[relative].get("size_bytes")
                ):
                    raise ValueError("Source differs from the preserved dossier")
                return self._reply(200, raw, "application/json; charset=utf-8" if path.name == "dossier.json" else "text/plain; charset=utf-8")
            return self._reply(404, {"error": "Unknown local route"})
        except ChatGPTError as error:
            return self._reply(*account_error(error))
        except (ValueError, OSError, ImportError):
            return self._reply(400, {"error": "Account or report unavailable; inspect local CLI diagnostics"})

    def do_POST(self):
        if not self._origin_ok(mutation=True):
            return self._reply(403, {"error": "Local same-origin consent required"})
        try:
            if self.headers.get("Transfer-Encoding") or self.headers.get("Content-Type", "").split(";")[0] != "application/json":
                raise ValueError("JSON content required")
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 4096:
                raise ValueError("Request exceeds local limit")
            body = json.loads(self.rfile.read(length))
            if not isinstance(body, dict):
                raise ValueError("Expected an object")
            if self.path == "/api/login-cancel":
                if body:
                    raise ValueError("Cancel accepts no supplied credentials")
                with self.server.login_guard:
                    cancel, done = self.server.login_cancel, self.server.login_done
                    if cancel is not None:
                        cancel.set()
                ended = done.wait(timeout=2) if done is not None else True
                return self._reply(200, {"cancelled": cancel is not None and ended,
                                         "login_in_progress": not ended})
            if not self.server.operation.acquire(blocking=False):
                return self._reply(409, {"error": "Another account or collection operation is active"})
            login_done = None
            try:
                if self.path == "/api/login":
                    if set(body) - {"add_account"} or ("add_account" in body and not isinstance(body["add_account"], bool)):
                        raise ValueError("Only an optional add_account boolean is accepted")
                    cancel, login_done = threading.Event(), threading.Event()
                    with self.server.login_guard:
                        self.server.login_cancel, self.server.login_done = cancel, login_done
                    with self.server.client() as client:
                        options = {"timeout": 180, "cancel_event": cancel}
                        if body.get("add_account") is True:
                            options["add_account"] = True
                        return self._reply(200, client.sign_in(**options))
                if self.path == "/api/logout":
                    if body:
                        raise ValueError("No supplied credentials accepted")
                    with self.server.client() as client:
                        return self._reply(200, client.disconnect())
                if self.path == "/api/collect":
                    allowed = {"query", "kind", "taxon", "action", "model", "reasoning_effort", "review", "max_papers"}
                    if set(body) - allowed or not isinstance(body.get("review", False), bool):
                        raise ValueError("Unsupported literature request")
                    run_id = uuid.uuid4().hex
                    result = collect_local(**body, state_dir=self.server.state_dir,
                                           output_dir=self.server.output_root / run_id)
                    return self._reply(200, {"id": run_id, "status": result.get("status"),
                                             "dossier_sha256": hashlib.sha256((self.server.output_root / run_id / "dossier.json").read_bytes()).hexdigest()})
                return self._reply(404, {"error": "Unknown local route"})
            finally:
                self.server.operation.release()
                if login_done is not None:
                    with self.server.login_guard:
                        if self.server.login_done is login_done:
                            self.server.login_cancel = self.server.login_done = None
                    login_done.set()
        except ChatGPTError as error:
            return self._reply(*account_error(error))
        except (ValueError, OSError, ImportError, TypeError):
            return self._reply(400, {"error": "Operation did not complete. Check login, model, permissions and CLI diagnostics; no scientific validation was granted."})


def main(argv=None):
    parser = argparse.ArgumentParser(description="SkinPipe local literature evidence and ChatGPT account connection")
    parser.add_argument("--state-dir", type=Path, default=default_state_dir())
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("status")
    login = commands.add_parser("login")
    login.add_argument("--add-account", action="store_true",
                       help="Connect another ChatGPT account; keep the current connection until the new sign-in succeeds")
    commands.add_parser("logout")
    commands.add_parser("models")
    collect = commands.add_parser("collect")
    collect.add_argument("query")
    collect.add_argument("--out", type=Path, required=True, help="New directory; existing directories are never overwritten")
    collect.add_argument("--kind", choices=["protein", "complex", "endpoint"], default="protein")
    collect.add_argument("--taxon", type=int, default=9606)
    collect.add_argument("--action", choices=["any", "inhibit", "activate"], default="any")
    collect.add_argument("--max-papers", type=int, default=8)
    collect.add_argument("--review", action="store_true")
    collect.add_argument("--model", default="gpt-6-luna")
    collect.add_argument("--reasoning-effort", choices=["none", "low", "medium", "high", "xhigh", "max"], default="max")
    serve = commands.add_parser("serve")
    serve.add_argument("--port", type=int, default=8771)
    serve.add_argument("--out-root", type=Path, default=Path.home() / ".local" / "share" / "skinpipe" / "literature")
    serve.add_argument("--open-browser", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "serve":
            with LiteratureServer(("127.0.0.1", args.port), state_dir=args.state_dir, output_root=args.out_root) as server:
                print(server.origin, flush=True)
                if args.open_browser:
                    webbrowser.open(server.origin)
                try:
                    server.serve_forever()
                except KeyboardInterrupt:
                    pass
            return
        if args.command == "collect":
            result = collect_local(args.query, args.out, state_dir=args.state_dir, model=args.model,
                                   reasoning_effort=args.reasoning_effort,
                                   review=args.review, max_papers=args.max_papers, kind=args.kind,
                                   taxon=args.taxon, action=args.action)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return
        from .chatgpt import ChatGPTClient
        with ChatGPTClient(args.state_dir) as client:
            if args.command == "login":
                result = client.sign_in(add_account=True) if args.add_account else client.sign_in()
            else:
                result = {"status": client.status, "logout": client.disconnect,
                          "models": client.list_models}[args.command]()
        print(json.dumps(result, ensure_ascii=False, indent=2))
    except ChatGPTError as error:
        _, safe = account_error(error)
        parser.exit(2, json.dumps(safe, ensure_ascii=False) + "\n")
    except (ValueError, OSError, ImportError) as error:
        # Transport exceptions and token responses are never printed.
        parser.exit(2, "SkinPipe literature operation unavailable (" + type(error).__name__ + "). Check local storage, login, selected model and the literature guide.\n")


if __name__ == "__main__":
    main()
