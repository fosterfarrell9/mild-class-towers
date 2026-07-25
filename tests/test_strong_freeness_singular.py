#!/usr/bin/env python3
"""Tests for the Singular/Letterplace strong-freeness driver."""

from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "tools" / "strong_freeness_singular.py"
SPEC = importlib.util.spec_from_file_location(
    "strong_freeness_singular", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
driver = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = driver
SPEC.loader.exec_module(driver)
gb = driver.gb

UT_PHI_TEXT = (
    "[0,2,0,1,4,0,0,2,1,2,2,3,4,0,0,2,0,0,0,3,3,0,0,0,1,0,0;"
    "0,3,4,4,2,0,2,1,0,3,1,4,2,0,1,1,3,0,4,4,0,0,1,0,0,0,0;"
    "0,1,2,3,2,1,1,4,0,1,1,0,2,0,0,4,0,0,2,0,0,1,0,0,0,0,0]")
UT_PHI = gb.parse_matrix_text(UT_PHI_TEXT)

# 1. Script generation is deterministic and uses the descending variable
#    list matching the ascending --order convention.
script = driver.build_script(UT_PHI, "xyz", 6)
assert "ring r = 5,(z,y,x),Dp;" in script
assert "freeAlgebra(r, 6)" in script
assert script.count("GB_LEADS") == 1 and script.count("GB_END") == 1

# 2. Lead parsing tolerates coefficients and blank lines.
parsed = driver.parse_leads(
    "noise\nGB_LEADS\nx*y*y\n4*z*z*x\n\nGB_END\ntail")
assert parsed == ["xyy", "zzx"], parsed

# 3. If Singular is available, the fixture must certify strong freeness.
if shutil.which("Singular"):
    code = driver.main([
        "--matrix", UT_PHI_TEXT, "--order", "xyz", "--degree-bound", "6"])
    assert code == 0
else:
    print("Singular not found; run-level check skipped")

print("STRONG_FREENESS_SINGULAR_TEST PASS")
