# Verify GitHub Actions for a pushed commit

Use the exact remote commit and include `push` runs. A pull-request-only query
can return no runs for a successful push; an empty response does not prove
success or failure.

With GitHub CLI installed and authenticated:

```sh
REPO=kangk1204/SkinPipe_public
COMMIT=$(git rev-parse HEAD)
gh run list --repo "$REPO" --commit "$COMMIT" --event push --limit 100 \
  --json databaseId,workflowName,headSha,event,status,conclusion,url
gh run view RUN_ID --repo "$REPO" --json headSha,event,status,conclusion,jobs,url
```

Replace `RUN_ID` with the returned run ID. Check every expected workflow for
the same full `headSha`, wait for completion, and inspect failed job annotations.
The repository's **Actions** page is an alternative when CLI access is unavailable.
The REST endpoint `GET /repos/{owner}/{repo}/actions/runs?head_sha={sha}` also
includes push runs; follow pagination and retain the event and SHA in the record.

Report these states separately:

- **Passed:** the expected workflow jobs actually ran and completed successfully.
- **Pending:** an expected run is queued or running.
- **Failed:** a job ran and failed; inspect its log before changing code.
- **Blocked before execution:** an annotation says the job could not start, for
  example because of an account billing or spending limit. Record the reason;
  do not describe this as a test failure or a successful CI check.
- **Unverified:** an expected run is absent or cannot be read with current access.

Independent Ubuntu checks can continue when hosted execution is blocked. Record
the OS, Python version, tested source revision, commands, results and skipped
checks separately. They do not replace a missing hosted run. Do not disable
required tests, change billing, expose private contents or register a new runner
solely to turn a blocked check green.
