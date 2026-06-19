"""Unit tests for the de-hardwired path configuration in ``landsymm.config``.

These verify that the LANDSYMM_* environment variables drive path/naming
resolution (so a new scenario set or data location needs no source edits), and
that the baseline-file, harmonized-dir-name, and scenario-discovery helpers
behave as expected. Portable: no external data required.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest  # noqa: F401  (pytest provides tmp_path / monkeypatch fixtures)

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from landsymm import config  # noqa: E402

_ENV_VARS = [
    "LANDSYMM_DATA_DIR",
    "LANDSYMM_GEODATA_DIR",
    "LANDSYMM_REMAP_DIRNAME",
    "LANDSYMM_REMAP_VER",
    "LANDSYMM_PLUM_DIRNAME",
    "LANDSYMM_MEMBER",
]


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Start each test from a clean slate so an ambient LANDSYMM_* in the calling
    shell cannot leak into assertions about defaults."""
    for _v in _ENV_VARS:
        monkeypatch.delenv(_v, raising=False)


def test_env_overrides(tmp_path, monkeypatch):
    monkeypatch.setenv("LANDSYMM_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("LANDSYMM_REMAP_DIRNAME", "myremap")
    monkeypatch.setenv("LANDSYMM_REMAP_VER", "9_test")
    monkeypatch.setenv("LANDSYMM_PLUM_DIRNAME", "myplum")
    monkeypatch.setenv("LANDSYMM_MEMBER", "s2")

    assert config.get_data_dir() == tmp_path
    assert config.get_remap_output_dir() == tmp_path / "myremap"
    assert config.get_remap_baseline_dir() == tmp_path / "myremap" / "remaps_v9_test"
    bf = config.get_remap_baseline_files()
    assert bf["lu"].name == "LU.remapv9_test.txt"
    assert bf["cropf"].name == "cropfracs.remapv9_test.txt"
    assert bf["nfert"].name == "nfert.remapv9_test.txt"
    assert config.get_plum_output_dir() == tmp_path / "myplum"
    assert config.get_member() == "s2"
    assert config.get_geodata_dir() == tmp_path / "geodata_py"


def test_defaults(tmp_path, monkeypatch):
    monkeypatch.setenv("LANDSYMM_DATA_DIR", str(tmp_path))
    for v in _ENV_VARS:
        if v != "LANDSYMM_DATA_DIR":
            monkeypatch.delenv(v, raising=False)
    assert config.get_remap_dirname() == "output_hildaplus_remap_10b"
    assert config.get_remap_ver() == "10_old_62892_gL"
    assert config.get_member() == "s1"
    assert config.get_plum_output_dir() == tmp_path / "PLUMv2_LU_default_output"


def test_harm_dirname(monkeypatch):
    for v in ["LANDSYMM_MEMBER", "LANDSYMM_REMAP_VER"]:
        monkeypatch.delenv(v, raising=False)
    assert config.harm_dirname() == "s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py"
    assert (
        config.harm_dirname(member="s5", remap_ver="x", allow_unveg=False)
        == "s5.HILDA+_remap_vx.harm_py"
    )


def test_geodata_override(tmp_path, monkeypatch):
    gd = tmp_path / "elsewhere_geodata"
    gd.mkdir()
    monkeypatch.setenv("LANDSYMM_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("LANDSYMM_GEODATA_DIR", str(gd))
    assert config.get_geodata_dir() == gd
    monkeypatch.delenv("LANDSYMM_GEODATA_DIR")
    assert config.get_geodata_dir() == tmp_path / "geodata_py"


def test_discover_scenarios(tmp_path, monkeypatch):
    monkeypatch.setenv("LANDSYMM_MEMBER", "s1")
    for s in ["BAU", "NfN_NfN", "NaC_NaC"]:
        (tmp_path / s / "s1").mkdir(parents=True)
    (tmp_path / "README.txt").write_text("x")          # a file -> ignored
    (tmp_path / "no_member_dir").mkdir()                 # dir without s1 -> ignored
    assert config.discover_scenarios(tmp_path) == ["BAU", "NaC_NaC", "NfN_NfN"]
    assert config.discover_scenarios(tmp_path / "does_not_exist") == []
