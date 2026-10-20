"""Writable user cache locations, independent of checkout/package permissions."""
from __future__ import annotations
import os
from pathlib import Path


def user_cache_dir():
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")).expanduser() / "skinpipe"
