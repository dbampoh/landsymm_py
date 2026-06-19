"""Integration test for the PLUMharmFigs diagnostic-figure module.

Runs ``run_plumharm_figs`` on a small year window for one scenario and asserts
that time-series PDFs and map PNGs are produced. This exercises the figs code
path that previously failed (see the four bugs fixed in plumharm_figs.py /
plumharm_import_ref.py: the None baseline-year default, the bad_base_yx
broadcast, and the advanced-index axis reordering at the regional-stats and
management-time-series call sites).

The test needs harmonized scenario outputs on disk and is SKIPPED automatically
when they are not present (e.g., in CI). To run it, set the env vars below to a
data location that contains harmonized output for the scenario (defaults point
at the paper-3 NFF data on this workstation).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# Defaults target the paper-3 NFF harmonized data on this workstation; override
# via the environment to point elsewhere.
_DATA_DIR = os.environ.get(
    "LANDSYMM_TEST_DATA_DIR",
    "/home/bampoh-d/Desktop/landsymm_lpjg/landsymm_mat/data",
)
_REMAP_DIRNAME = os.environ.get("LANDSYMM_TEST_REMAP_DIRNAME", "output_hildaplus_remap_10b_3_py")
_PLUM_DIRNAME = os.environ.get("LANDSYMM_TEST_PLUM_DIRNAME", "PLUMv2_LU_NFF_output")
_SCENARIO = os.environ.get("LANDSYMM_TEST_SCENARIO", "BAU")


@pytest.fixture
def _env(monkeypatch):
    if not Path(_DATA_DIR).is_dir():
        pytest.skip(f"Data dir not present: {_DATA_DIR}")
    monkeypatch.setenv("LANDSYMM_DATA_DIR", _DATA_DIR)
    monkeypatch.setenv("LANDSYMM_REMAP_DIRNAME", _REMAP_DIRNAME)
    monkeypatch.setenv("LANDSYMM_PLUM_DIRNAME", _PLUM_DIRNAME)
    monkeypatch.delenv("LANDSYMM_MEMBER", raising=False)
    monkeypatch.delenv("LANDSYMM_REMAP_VER", raising=False)


def test_plumharm_figs_runs(_env):
    from landsymm import config
    from landsymm.harmonization.run_plumharm_figs import main as figs_main

    data_dir = config.get_plum_output_dir()
    harm_sub = config.harm_dirname()
    harm_dir = data_dir / _SCENARIO / harm_sub
    if not harm_dir.is_dir() or not config.get_remap_baseline_files()["lu"].exists():
        pytest.skip(
            f"Harmonized data not available for {_SCENARIO} at {harm_dir}; "
            "run the harmonization pipeline first."
        )

    # Small 2-year window keeps the test quick; reads the existing harmonized data.
    figs_main(scenarios=[_SCENARIO], year1=2021, yearN=2022)

    figs_dir = data_dir / _SCENARIO / f"{harm_sub}_harms_figs"
    pdfs = list(figs_dir.glob("timeSeries_*.pdf"))
    pngs = list(figs_dir.glob("maps_*.png"))
    xlsx = list(figs_dir.glob("harm_by_numbers.*.xlsx"))
    assert pdfs, f"no time-series PDFs produced in {figs_dir}"
    assert pngs, f"no map PNGs produced in {figs_dir}"
    assert xlsx, f"no regional-statistics spreadsheets produced in {figs_dir}"
