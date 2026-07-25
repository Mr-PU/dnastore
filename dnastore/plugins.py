"""
Third-party codec plugin discovery.

DNAstore's built-in codecs (naive, rotating, fountain) live in
`dnastore.codec.REGISTRY`. This module lets *other, separately
installed packages* register additional codecs without ever touching
this repository -- the same mechanism pytest uses for plugins, mkdocs
uses for themes, and much of the scientific Python ecosystem uses for
extensibility in general: Python entry points.

A third-party package declares codecs in its own pyproject.toml:

    [project.entry-points."dnastore.codecs"]
    my_codec = "my_package.module:MyCodecClass"

Once that package is `pip install`-ed alongside dnastore, its codec
shows up automatically:

    import dnastore
    dnastore.available_codecs()   # includes "my_codec"
    dnastore.DNAStorage(codec="my_codec")

No PR against this repo, no waiting on a maintainer. See
`examples/plugin_example/` for a complete, working example package
(built and installed in this project's own test suite to prove the
mechanism actually works, not just described).
"""
from __future__ import annotations

import functools
from importlib.metadata import entry_points

PLUGIN_GROUP = "dnastore.codecs"


class PluginLoadError(Exception):
    """Raised when a registered entry point exists but fails to import --
    kept distinct from a plain ImportError so callers can tell "this
    codec name doesn't exist" apart from "this codec exists but its
    package is broken"."""


@functools.lru_cache(maxsize=1)
def discover_plugin_codecs() -> dict[str, "type"]:
    """Scan installed packages for dnastore.codecs entry points. Cached
    after first call -- installing a new plugin package requires a
    fresh Python process to pick it up, same as any other Python
    entry-point-based plugin system (pytest, mkdocs, etc.)."""
    discovered: dict[str, type] = {}
    try:
        eps = entry_points(group=PLUGIN_GROUP)
    except TypeError:
        # Python < 3.10 compatibility path (entry_points() took no kwargs)
        eps = entry_points().get(PLUGIN_GROUP, [])

    for ep in eps:
        try:
            discovered[ep.name] = ep.load()
        except Exception as e:
            raise PluginLoadError(
                f"failed to load codec plugin '{ep.name}' from {ep.value}: {e}"
            ) from e
    return discovered


def clear_plugin_cache() -> None:
    """Mainly for tests: force re-discovery within the same process
    (e.g. after installing a plugin package mid-test-run)."""
    discover_plugin_codecs.cache_clear()
