"""Personal bearer credentials and project roles; external headers grant nothing."""

from __future__ import annotations

import hashlib
import hmac
import secrets
from psycopg.types.json import Jsonb

from .queue import identifier

ROLES = {"viewer", "researcher", "admin"}


class ProjectAccess:
    def __init__(self, queue, config=None):
        self.queue = queue
        self.credentials = (config or {}).get("credentials", [])
        for credential in self.credentials:
            if not isinstance(credential, dict) or credential.get("role") not in ROLES:
                raise ValueError("Personal credential requires a project role")
            identifier(credential.get("user_id"), "user")
            identifier(credential.get("project_id"), "project")
            identifier(credential.get("owner_id"), "owner")
            sha = credential.get("token_sha256")
            if (
                not isinstance(sha, str)
                or len(sha) != 64
                or any(c not in "0123456789abcdef" for c in sha)
            ):
                raise ValueError(
                    "Store only SHA256 of high-entropy personal bearer tokens"
                )
        namespaces = {}
        for credential in self.credentials:
            owner, project = credential["owner_id"], credential["project_id"]
            if owner in namespaces and namespaces[owner] != project:
                raise ValueError(
                    "Each project requires a distinct research owner namespace"
                )
            namespaces[owner] = project
        if queue is not None:
            with queue.connection() as conn:
                conn.execute("""CREATE TABLE IF NOT EXISTS platform_project_credentials (
                    user_id TEXT NOT NULL, project_id TEXT NOT NULL, owner_id TEXT NOT NULL,
                    role TEXT NOT NULL CHECK(role IN ('viewer','researcher','admin')),
                    token_sha256 TEXT PRIMARY KEY, revoked BOOLEAN NOT NULL DEFAULT FALSE,
                    created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp())""")
                conn.execute("""CREATE TABLE IF NOT EXISTS platform_access_events (
                    id BIGSERIAL PRIMARY KEY,project_id TEXT NOT NULL,actor_id TEXT NOT NULL,
                    action TEXT NOT NULL,data JSONB NOT NULL,created_at TIMESTAMPTZ NOT NULL DEFAULT clock_timestamp())""")
                for owner, project in namespaces.items():
                    conn.execute("SELECT pg_advisory_xact_lock(hashtextextended(%s,0))", ("access-owner:" + owner,))
                    conflict = conn.execute("SELECT 1 FROM platform_project_credentials WHERE owner_id=%s AND project_id<>%s LIMIT 1", (owner, project)).fetchone()
                    if conflict:
                        raise ValueError("Configured project namespace collides with durable credentials")

    def authenticate(self, token):
        if not isinstance(token, str) or not token:
            return None
        actual = hashlib.sha256(token.encode()).hexdigest()
        if self.queue is not None:
            with self.queue.connection() as conn:
                row = conn.execute(
                    "SELECT * FROM platform_project_credentials WHERE token_sha256=%s",
                    (actual,),
                ).fetchone()
            if row is not None:
                return (
                    None
                    if row["revoked"]
                    else {
                        k: row[k]
                        for k in (
                            "user_id",
                            "project_id",
                            "owner_id",
                            "role",
                            "token_sha256",
                        )
                    }
                )
        return next(
            (
                dict(c)
                for c in self.credentials
                if hmac.compare_digest(actual, c["token_sha256"])
            ),
            None,
        )

    def manage(self, principal, method, body=None):
        """Explicit administrator actions only, confined to their own project."""
        self.authorize(principal, method, "/v2/access", body)
        if self.queue is None:
            raise RuntimeError("Durable access ledger unavailable")
        project, owner = principal["project_id"], principal["owner_id"]
        if method == "GET":
            with self.queue.connection() as conn:
                rows = conn.execute(
                    "SELECT user_id,role,revoked,created_at FROM platform_project_credentials WHERE project_id=%s ORDER BY created_at",
                    (project,),
                ).fetchall()
            return {
                "project_id": project,
                "members": [
                    dict(r, created_at=r["created_at"].isoformat()) for r in rows
                ],
                "configured_members": [
                    {k: c[k] for k in ("user_id", "role")}
                    for c in self.credentials
                    if c["project_id"] == project
                ],
            }
        if not isinstance(body, dict) or set(body) - {
            "action",
            "user_id",
            "role",
            "project_id",
        }:
            raise ValueError("Unsupported access action")
        user = identifier(body.get("user_id"), "user")
        action = body.get("action")
        if action not in {"issue", "revoke", "role"} or (
            action != "revoke" and body.get("role") not in ROLES
        ):
            raise ValueError("Explicit issue, revoke or role action required")
        if user == principal["user_id"] and action in {"revoke", "role"}:
            raise ValueError("An administrator cannot remove their own active access")
        token = secrets.token_urlsafe(48) if action == "issue" else None
        with self.queue.connection() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(hashtextextended(%s,0))", ("access-owner:" + owner,))
            # The same chemical/results namespace is never assigned to a
            # different project, including credentials loaded from the config.
            conflict = conn.execute(
                "SELECT 1 FROM platform_project_credentials WHERE owner_id=%s AND project_id<>%s LIMIT 1",
                (owner, project),
            ).fetchone()
            if conflict:
                raise ValueError("Project namespace collision")
            if action == "issue":
                conn.execute(
                    "INSERT INTO platform_project_credentials(user_id,project_id,owner_id,role,token_sha256) VALUES(%s,%s,%s,%s,%s)",
                    (
                        user,
                        project,
                        owner,
                        body["role"],
                        hashlib.sha256(token.encode()).hexdigest(),
                    ),
                )
            else:
                for c in self.credentials:
                    if c["project_id"] == project and c["user_id"] == user:
                        conn.execute(
                            "INSERT INTO platform_project_credentials(user_id,project_id,owner_id,role,token_sha256) VALUES(%s,%s,%s,%s,%s) ON CONFLICT DO NOTHING",
                            (user, project, owner, c["role"], c["token_sha256"]),
                        )
                if action == "revoke":
                    conn.execute(
                        "UPDATE platform_project_credentials SET revoked=true WHERE project_id=%s AND user_id=%s",
                        (project, user),
                    )
                else:
                    conn.execute(
                        "UPDATE platform_project_credentials SET role=%s WHERE project_id=%s AND user_id=%s",
                        (body["role"], project, user),
                    )
            conn.execute(
                "INSERT INTO platform_access_events(project_id,actor_id,action,data) VALUES(%s,%s,%s,%s)",
                (
                    project,
                    principal["user_id"],
                    action,
                    Jsonb({"user_id": user, "role": body.get("role")}),
                ),
            )
        return {
            "project_id": project,
            "user_id": user,
            "action": action,
            **(
                {"token": token, "token_disclosure": "one-time; store privately"}
                if token
                else {}
            ),
        }

    def authorize(self, principal, method, url, body=None):
        from urllib.parse import urlparse, unquote

        parts = [unquote(p) for p in urlparse(url).path.strip("/").split("/")]
        if not principal:
            raise PermissionError("Personal authorization required")
        if parts[:2] != ["v2", "access"] and method == "GET":
            return principal["owner_id"]
        if principal["role"] == "viewer":
            # Search and resolution do not mutate research sources or run compute.
            if method == "POST" and parts in (
                ["v2", "candidate-queries"],
                ["v2", "catalog-queries"],
                ["v2", "targets", "resolve"],
                ["v2", "research-inputs", "resolve"],
                ["v2", "exports"],
            ):
                return principal["owner_id"]
            raise PermissionError("Viewer role cannot modify research data")
        if parts[:2] == ["v2", "access"] and principal["role"] != "admin":
            raise PermissionError("Administrator role required")
        if body and body.get("project_id") not in {None, principal["project_id"]}:
            raise PermissionError("Request belongs to another project")
        return principal["owner_id"]
