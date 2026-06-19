"""Unit tests for wetland_into_forLPJG path selection (`_resolve_landcover_files`).

Verifies the three input modes (explicit files / flat-parent scan / default
forLPJG scan) and the scenario filter, without needing any wetland data.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from landsymm.wetlands.wetland_into_forLPJG import _resolve_landcover_files  # noqa: E402


def _touch(p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("Lon Lat Year PASTURE CROPLAND NATURAL BARREN\n")


def test_explicit_landcover(tmp_path):
    a = tmp_path / "a" / "landcover.txt"
    b = tmp_path / "b" / "landcover.txt"
    _touch(a)
    _touch(b)
    out = _resolve_landcover_files(landcover=[str(a), str(b)])
    assert out == [str(a.resolve()), str(b.resolve())]


def test_explicit_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        _resolve_landcover_files(landcover=[str(tmp_path / "nope" / "landcover.txt")])


def test_flat_parent_scan(tmp_path):
    for s in ["av", "bio", "pv"]:
        _touch(tmp_path / s / "landcover.txt")
    out = _resolve_landcover_files(flat_parent=str(tmp_path))
    assert [Path(p).parent.name for p in out] == ["av", "bio", "pv"]


def test_forlpjg_scan_default_member(tmp_path, monkeypatch):
    monkeypatch.delenv("LANDSYMM_MEMBER", raising=False)  # default member s1
    for s in ["BAU", "NfN_NfN"]:
        _touch(tmp_path / s / "s1.HILDA+_remap_vX.harm.allow_unveg_py.forLPJG" / "landcover.txt")
    # a non-matching dir (wrong member) should be ignored
    _touch(tmp_path / "BAU" / "s2.something.forLPJG" / "landcover.txt")
    out = _resolve_landcover_files(parent_dir=str(tmp_path))
    assert len(out) == 2
    assert all("s1." in p and p.endswith("landcover.txt") for p in out)


def test_scenarios_filter_applies_to_scan_only(tmp_path):
    for s in ["av", "bio", "pv"]:
        _touch(tmp_path / s / "landcover.txt")
    out = _resolve_landcover_files(flat_parent=str(tmp_path), scenarios=["bio"])
    assert [Path(p).parent.name for p in out] == ["bio"]
