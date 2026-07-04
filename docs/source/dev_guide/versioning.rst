Versioning & Releases
=====================

The single source of truth for the version is **git tags on this repository**.
There is no version string to bump by hand: the package version is derived from
the most recent tag by `setuptools_scm`_, and the documentation and version
switcher follow automatically.

.. _setuptools_scm: https://setuptools-scm.readthedocs.io/

How the version is derived
--------------------------

``pyproject.toml`` declares ``dynamic = ["version"]`` and configures
``scikit-build-core`` to obtain the version from ``setuptools_scm`` at build
time:

- **On a release tag** ``vX.Y.Z`` → clean version ``X.Y.Z``.
- **Between tags** → a development version such as ``0.2.1.dev5+g1a2b3c4``,
  where ``dev5`` counts commits since the last tag and ``g1a2b3c4`` is the
  current commit. This is the "version you are working on".

At runtime the version is importable::

    >>> import lrgsglib
    >>> lrgsglib.__version__        # doctest: +SKIP
    '0.2.1.dev5+g1a2b3c4'

``setuptools_scm`` writes ``src/lrgsglib/_version.py`` on every build (it is
git-ignored). If the repository has no tags at all, the build falls back to the
``fallback_version`` set in ``pyproject.toml``.

Cutting a release
-----------------

Releases are made purely by tagging. Semantic versioning (``MAJOR.MINOR.PATCH``)
is used.

.. code-block:: bash

   # 1. Make sure main is at the commit you want to release.
   git checkout main && git pull

   # 2. Create an annotated tag (the "v" prefix is required).
   git tag -a v0.2.0 -m "Release 0.2.0"

   # 3. Push the tag. This is what triggers the release.
   git push origin v0.2.0

Pushing the tag triggers the documentation deploy (below); publishing to PyPI,
if desired, is handled by the wheel-building workflow.

Documentation versions
-----------------------

The docs are deployed to GitHub Pages by ``.github/workflows/docs-deploy.yml``:

- **push to** ``main`` → builds and publishes the in-development docs at
  ``<site>/dev/``.
- **push a tag** ``vX.Y.Z`` → builds and publishes that release at
  ``<site>/X.Y.Z/``.

On every deploy, ``docs/gen_switcher.py`` regenerates ``switcher.json`` from the
git tags, so the **version switcher** in the navbar always lists ``dev`` plus
every released version, with the newest release marked ``(stable)``. The site
root redirects to the newest release.

Building the switcher locally
-----------------------------

The switcher dropdown is populated from a JSON file fetched at the URL given by
``DOCS_SWITCHER_URL`` (or ``DOCS_BASEURL``). For a local ``make livehtml`` the
default points at the local server, so the dropdown works offline against the
committed ``docs/source/_static/switcher.json``. To preview the tag-derived
list without deploying::

    python docs/gen_switcher.py --base-url http://127.0.0.1:8000
