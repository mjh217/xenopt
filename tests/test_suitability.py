import numpy as np
import pandas as pd

from xenopt import cations_per_formula_unit, cpx_tests, temperature_test

from .sample import NAMES, PERIDOTITE_CPX

GOOD = PERIDOTITE_CPX


def run(*rows, **thresholds):
    wt = pd.DataFrame(list(rows), columns=NAMES, dtype=float)
    return cpx_tests(wt, cations_per_formula_unit(wt, "cpx"), **thresholds)


def test_good_cpx_passes_all():
    tests = run(GOOD)
    assert tests.iloc[0].all(), tests.iloc[0][~tests.iloc[0]]


def test_missing_cpx_fails_all():
    tests = run(dict.fromkeys(NAMES, np.nan))
    assert not tests.iloc[0].any()


def test_nimis98_cr_al_field():
    cases = [  # (Cr2O3, Al2O3, passes)
        (0.4, 2.0, False),  # Cr2O3 too low
        (0.5, 2.0, False),  # the boundary is Cr2O3 > 0.5 (Ziberna et al., 2016)
        (1.0, 4.0, True),
        (1.0, 4.1, False),
        (2.25, 5.0, True),
        (2.25, 5.1, False),
    ]
    tests = run(*[{**GOOD, "Cr2O3": cr, "Al2O3": al} for cr, al, _ in cases])
    assert tests["nimis98_cr_al"].tolist() == [passes for *_, passes in cases]
    # Optionally Cr2O3 >= 0.5
    inclusive = run(*[{**GOOD, "Cr2O3": cr, "Al2O3": 2.0} for cr in (0.4, 0.5)], cr2o3_min_inclusive=True)
    assert inclusive["nimis98_cr_al"].tolist() == [False, True]
    assert not run({**GOOD, "Cr2O3": 0.8}, cr2o3_min=1.0)["nimis98_cr_al"][0]


def test_thresholds():
    # GOOD has Ca/(Ca+Mg) ~ 0.48, a_CaCrTs/Cr# ~ 0.05 and ~ 3.99 cations
    tests = run(GOOD, ca_mg_max=0.45, cacrts_cr_min=0.1, cation_range=(4.01, 4.03))
    assert not tests.loc[0, ["ziberna16_ca_mg", "sudholz21_cacrts", "cpx_cation_total"]].any()


def test_temperature_test():
    assert temperature_test(pd.Series([599.9, 600.0, np.nan])).tolist() == [False, True, False]
    assert not temperature_test(pd.Series([650.0]), T_min=700)[0]
