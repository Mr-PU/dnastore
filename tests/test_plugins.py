"""Tests for entry-point-based codec plugin discovery.

Note: these tests assume the demo plugin package in
examples/plugin_example/ has been installed into the test environment
(see CONTRIBUTING.md / README "Plugin architecture" section for the
install command). If it isn't installed, the plugin-specific tests are
skipped rather than failed -- absence of an optional demo package
shouldn't break the core test suite.
"""
import os
import pytest

import dnastore
from dnastore.codec import available_codecs, get_codec
from dnastore.plugins import discover_plugin_codecs, clear_plugin_cache, PLUGIN_GROUP

PLUGIN_INSTALLED = "xor_demo" in available_codecs()


def test_builtin_codecs_always_present():
    codecs = dnastore.available_codecs()
    assert {"naive", "rotating", "fountain"}.issubset(set(codecs))


def test_discover_plugin_codecs_returns_dict():
    result = discover_plugin_codecs()
    assert isinstance(result, dict)


@pytest.mark.skipif(not PLUGIN_INSTALLED, reason="demo plugin package not installed")
def test_demo_plugin_discovered():
    assert "xor_demo" in dnastore.available_codecs()


@pytest.mark.skipif(not PLUGIN_INSTALLED, reason="demo plugin package not installed")
def test_demo_plugin_codec_round_trip():
    codec = get_codec("xor_demo")
    data = os.urandom(500)
    strands = codec.encode(data)
    decoded = codec.decode(strands, len(data))
    assert decoded == data


@pytest.mark.skipif(not PLUGIN_INSTALLED, reason="demo plugin package not installed")
def test_demo_plugin_full_api_round_trip():
    from dnastore import DNAStorage
    dna = DNAStorage(codec="xor_demo")
    data = os.urandom(500)
    dna.store(data, name="f.bin")
    assert dna.retrieve("f.bin") == data


def test_unknown_codec_raises_with_helpful_message():
    with pytest.raises(ValueError, match="unknown codec"):
        get_codec("definitely-not-a-real-codec")


def test_plugin_cache_can_be_cleared():
    discover_plugin_codecs()  # populate cache
    clear_plugin_cache()      # should not raise
    discover_plugin_codecs()  # repopulate
