"""Required oxides, opx-garnet tests and warnings."""

import numpy as np
import pandas as pd
import pytest

from xenopt import opx_garnet_tests, pair_tests, solve_pt, warning_text
from xenopt.chemistry import C2K
from xenopt.thermobarometers import T_brey_kohler90_ca_in_opx

from .sample import sample


def test_ca_in_opx_correction():
    _, c = sample()
    T = float(T_brey_kohler90_ca_in_opx(c, 4.0)[0]) - C2K
    corrected = float(T_brey_kohler90_ca_in_opx(c, 4.0, correction="NG10")[0]) - C2K
    assert corrected == pytest.approx(-628.7 + 2.0690 * T - 4.530e-4 * T**2)


def test_missing_required_oxide():
    wt, c = sample(cpx=["Na2O"])
    pt = solve_pt(c, "T98", "NG85", wt_perc=wt)
    assert pt.status[0] == "missing oxides"
    assert pt.missing_oxides[0] == "cpx Na2O"
    assert np.isnan(pt.T_C[0])


def test_cr_is_required_only_where_it_matters():
    # Cr-in-cpx needs cpx Cr2O3; Al-in-opx and two-pyroxene T only warn without it
    wt, c = sample(cpx=["Cr2O3"], opx=["Cr2O3"], gt=["Cr2O3"])
    pt = solve_pt(c, "NT00", "SUD21", wt_perc=wt)
    assert pt.status[0] == "missing oxides" and pt.missing_oxides[0] == "cpx Cr2O3"
    pt = solve_pt(c, "T98", "NG85", wt_perc=wt)
    assert pt.status[0] == "converged"
    assert pt.assumed_zero[0] == "opx Cr2O3, cpx Cr2O3, gt Cr2O3"


def test_blank_minor_oxide_is_assumed_zero():
    wt, c = sample(opx=["TiO2", "MnO"])
    pt = solve_pt(c, "T98", "NG85", wt_perc=wt)
    assert pt.status[0] == "converged"
    assert pt.assumed_zero[0] == "opx TiO2, opx MnO"
    # The sample's cpx K2O is blank; SUD21 uses it, NG85 does not
    wt, c = sample()
    assert solve_pt(c, "NT00", "SUD21", wt_perc=wt).assumed_zero[0] == "cpx K2O"
    # Not checked without wt_perc
    assert "assumed_zero" not in solve_pt(c, "T98", "NG85")


def test_opx_garnet_tests_on_sample():
    wt, c = sample()
    pt = solve_pt(c, "T98", "NG85", wt_perc=wt)
    tests = opx_garnet_tests(wt, c, pt.P_GPa).iloc[0]
    # T98 = 908 C at 1.26 GPa: good analyses in equilibrium, but at too low P for garnet lherzolite.
    # Ca-in-opx (corrected) is 52 C and NG09 10 C from T98.
    assert tests.drop("ng10_P_min").all()
    assert not tests["ng10_P_min"]
    # Tighter tolerances flag disequilibrium
    tight = opx_garnet_tests(wt, c, pt.P_GPa, cpx_opx_tol=(50, 50, 50), gt_opx_tol=5).iloc[0]
    assert not tight["ng10_cpx_opx"] and not tight["ng10_gt_opx"]
    # The PTEXL garnet-cpx check is off unless a tolerance is given; Krogh88 is 105 C above T98
    assert "ptexl_gt_cpx" not in tests
    assert opx_garnet_tests(wt, c, pt.P_GPa, gt_cpx_tol=120).ptexl_gt_cpx[0]
    assert not opx_garnet_tests(wt, c, pt.P_GPa, gt_cpx_tol=90).ptexl_gt_cpx[0]


def test_pair_tests_depend_on_minerals():
    wt, c = sample()
    pt = solve_pt(c, "T98", "NG85", wt_perc=wt)
    cpx = pd.DataFrame({"cpx_present": [True], "sudholz21_na": [True]})
    assert list(pair_tests({"cpx"}, wt, c, pt, cpx).columns) == ["sudholz21_na", "ziberna16_T"]
    assert "ng10_cpx_opx" in pair_tests({"opx", "cpx", "gt"}, wt, c, pt, cpx)
    assert pair_tests({"opx", "gt"}, wt, c, pt, cpx).empty


def test_warning_text():
    pt = pd.DataFrame(
        {
            "status": ["converged", "converged", "missing oxides", "missing minerals"],
            "missing_oxides": ["", "", "cpx Na2O", ""],
            "assumed_zero": ["opx TiO2", "", "", ""],
        }
    )
    tests = pd.DataFrame({"a": [False, True, False, False], "b": [False, True, True, True]})
    # Tests without a label appear under their own name
    assert warning_text(tests, pt).tolist() == [
        "a; b; assumed zero: opx TiO2",
        "",
        "missing oxides: cpx Na2O",
        "",
    ]

    # Labelled tests are described in words
    tests = pd.DataFrame({"ziberna16_ca_mg": [False, True], "ptexl_na_k": [False, False]})
    pt = pd.DataFrame({"status": ["converged", "converged"], "assumed_zero": ["", "cpx NiO"]})
    assert warning_text(tests, pt).tolist() == ["high Ca/(Ca+Mg); low Na+K", "low Na+K; assumed zero: cpx NiO"]

    # The lherzolite pressure limit failing alone gets a caveat; failing both
    # pressure limits only reports the stricter one
    tests = pd.DataFrame({"ng10_P_min": [False, False], "ng10_P_min_pyroxenite": [True, False]})
    pt = pd.DataFrame({"status": ["converged", "converged"]})
    assert warning_text(tests, pt).tolist() == ["low P (ok if pyroxenite)", "very low P"]

    # Very low T makes the low-T accuracy caution redundant
    tests = pd.DataFrame({"ng10_T_min": [True, False], "ng10_T_accuracy": [False, False]})
    assert warning_text(tests, pt).tolist() == ["low T (accuracy not assured)", "very low T"]
