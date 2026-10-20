"""Bound an external calculation and reap its process group on interruption."""
from __future__ import annotations

import os
import signal
import subprocess
import sys


def _darwin_group_retired(group_id):
    """Darwin can return EPERM for a just-retired process group.

    A read-only process inventory distinguishes it from a live group we cannot
    signal. A failed inventory or any live member is never a cleanup success.
    """
    try:
        result = subprocess.run(["ps", "-axo", "pgid=,stat="],
                                capture_output=True, text=True, timeout=3)
    except (OSError, subprocess.TimeoutExpired):
        return False
    if result.returncode:
        return False
    members = [fields for line in result.stdout.splitlines()
               if (fields := line.split()) and fields[0] == str(group_id)]
    return all(len(fields) >= 2 and fields[1].startswith("Z") for fields in members)


def _stop_group(process):
    if os.name == "posix":
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(process.pid, sig)
            except ProcessLookupError:
                break
            except PermissionError:
                if sys.platform == "darwin" and process.poll() is not None and _darwin_group_retired(process.pid):
                    break
                raise
            if sig == signal.SIGTERM:
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    pass
    elif process.poll() is None:
        process.kill()
    process.wait()


def run_calculation(args, *, timeout, **kwargs):
    """Like ``subprocess.run``, with descendant cleanup on POSIX hosts.

    Files and logs are retained. A timeout is still an exception, never a
    successful partial result. Calculations must not daemonize out of the
    session they inherit. Windows uses the direct-child fallback.
    """
    if "start_new_session" in kwargs:
        raise ValueError("calculation session is managed by run_calculation")
    command = args
    if sys.platform.startswith("linux"):
        # Inner calculations have their own timeout groups. Tie their lifetime
        # to the supervisor too, so an outer bridge timeout cannot orphan them.
        command = [sys.executable, "-m", "skinpipe.processes", str(os.getpid()), *args]
    process = subprocess.Popen(command, start_new_session=os.name == "posix", **kwargs)
    try:
        process.wait(timeout=timeout)
    except BaseException:
        _stop_group(process)
        raise
    _stop_group(process)
    return subprocess.CompletedProcess(args, process.returncode)


def _parent_bound_exec():
    """Linux guard cleans the inner group even if its supervisor disappears."""
    import ctypes
    signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})
    expected_parent = int(sys.argv[1])
    libc = ctypes.CDLL(None, use_errno=True)
    child = []

    def terminate(_signum, _frame):
        if child:
            # Do not call Popen.wait in a signal handler: the interrupted wait
            # can own its waitpid lock. Kill immediately; the finally block
            # reaps the child after that interrupted wait has released it.
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(child[0].pid, sig)
                except ProcessLookupError:
                    break
        raise SystemExit(143)

    signal.signal(signal.SIGTERM, terminate)
    if libc.prctl(1, int(signal.SIGTERM), 0, 0, 0) != 0:  # PR_SET_PDEATHSIG
        raise OSError(ctypes.get_errno(), "cannot bind calculation lifetime to parent")
    # Cover parent death between Popen and prctl. Never start the calculation
    # if the original supervisor has already disappeared.
    if os.getppid() != expected_parent:
        raise SystemExit("calculation supervisor exited before startup")
    # Defer parent-death notification until the child handle is owned, closing
    # the small spawn/assignment race without preexec_fn in a threaded parent.
    previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM})
    try:
        child.append(subprocess.Popen([sys.executable, "-m", "skinpipe.processes", "--exec", *sys.argv[2:]],
                                      start_new_session=True))
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous)
    if os.getppid() != expected_parent:
        terminate(signal.SIGTERM, None)
    try:
        status = child[0].wait()
    finally:
        _stop_group(child[0])
    raise SystemExit(status if status >= 0 else 128 - status)


if __name__ == "__main__":
    if sys.argv[1] == "--exec":
        # A fork inherits the blocked signal mask. Reset it before executing
        # the actual calculation, including nested preparation supervisors.
        signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})
        os.execvpe(sys.argv[2], sys.argv[2:], os.environ)
    else:
        _parent_bound_exec()
