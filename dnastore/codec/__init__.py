from .base import Codec, Strand, BASES, MAX_STRAND_LENGTH
from .naive import NaiveCodec
from .rotating import RotatingCodec
from .fountain import FountainCodec, FountainDecodeError
from ..plugins import discover_plugin_codecs

# Built-in codecs. Third-party codecs are *not* added here -- they're
# discovered dynamically via entry points (see dnastore.plugins) so
# installing a plugin package never requires editing this file.
REGISTRY: dict[str, type[Codec]] = {
    "naive": NaiveCodec,
    "rotating": RotatingCodec,
    "fountain": FountainCodec,
}


def available_codecs() -> list[str]:
    """Built-in codec names plus any discovered via the dnastore.codecs
    entry-point group (third-party plugin packages)."""
    return sorted(set(REGISTRY) | set(discover_plugin_codecs()))


def get_codec(name: str, **kwargs) -> Codec:
    if name in REGISTRY:
        return REGISTRY[name](**kwargs)
    plugins = discover_plugin_codecs()
    if name in plugins:
        return plugins[name](**kwargs)
    raise ValueError(f"unknown codec '{name}'. Available: {available_codecs()}")


__all__ = [
    "Codec", "Strand", "BASES", "MAX_STRAND_LENGTH",
    "NaiveCodec", "RotatingCodec", "FountainCodec", "FountainDecodeError",
    "REGISTRY", "get_codec", "available_codecs",
]
