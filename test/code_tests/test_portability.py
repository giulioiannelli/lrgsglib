"""Portability guards (Phase 0): the library must import without a GPU and
without the ``LRGSG_*`` environment being configured.

Both checks run in a *subprocess* so they exercise a clean interpreter:
- cupy is forced unavailable (``sys.modules['cupy'] = None``) to emulate a
  non-CUDA machine; ``import lrgsglib`` must still succeed.
- all ``LRGSG_*`` variables are stripped from the environment to emulate a fresh
  checkout / installed wheel; ``lrgsglib.config.lrgsg_env`` must still resolve
  every path to an absolute location derived from the package.
"""

import os
import subprocess
import sys

import pytest

# Block cupy, then import the full top-level package.
_NO_GPU_SNIPPET = (
    "import sys; sys.modules['cupy'] = None;"
    " import lrgsglib;"
    " from lrgsglib import shared;"
    " assert shared.cp is None;"
    " print('ok')"
)

# With no LRGSG_* env, every path constant must be absolute and rooted at the
# package location (relocatable derivation).
_NO_ENV_SNIPPET = (
    "import lrgsglib.config.lrgsg_env as e;"
    " names = [n for n in dir(e) if n.startswith('LRGSG_')];"
    " assert names, 'no LRGSG_* constants';"
    " vals = [getattr(e, n) for n in names];"
    " assert all(v.is_absolute() for v in vals), 'non-absolute path';"
    " assert str(e.LRGSG_LIB).endswith('src/lrgsglib'), e.LRGSG_LIB;"
    " assert e.LRGSG_DATA == e.LRGSG_LLIB / 'data', e.LRGSG_DATA;"
    " print('ok')"
)


def _run(snippet: str, env: dict) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-c", snippet],
        capture_output=True,
        text=True,
        env=env,
    )


@pytest.mark.code
def test_import_without_gpu():
    """`import lrgsglib` succeeds when cupy is unavailable (no CUDA box)."""
    res = _run(_NO_GPU_SNIPPET, dict(os.environ))
    assert res.returncode == 0, res.stderr
    assert res.stdout.strip().endswith("ok")


@pytest.mark.code
def test_env_paths_derive_without_environment():
    """`lrgsg_env` resolves all paths from the package when no LRGSG_* env set."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("LRGSG_")}
    res = _run(_NO_ENV_SNIPPET, env)
    assert res.returncode == 0, res.stderr
    assert res.stdout.strip().endswith("ok")
