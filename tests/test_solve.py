import numpy as np
import pandas as pd
import pytest

from xenopt import BAROMETERS, THERMOMETERS, solve_pt
from xenopt.chemistry import CATIONS
from xenopt.thermobarometers import Calibration


def make_cations(rows: dict) -> dict:
    """Cations for each mineral; minerals missing from `rows` are NaN for every sample."""
    n = len(next(iter(rows.values())))
    blank = pd.DataFrame(np.nan, index=range(n), columns=list(CATIONS))
    return {m: pd.DataFrame(rows[m], columns=list(CATIONS)) if m in rows else blank for m in ["opx", "cpx", "gt"]}


@pytest.fixture
def slow_pressure(monkeypatch):
    """T barely depends on P, while P converges slowly towards 0 GPa (ratio 0.9 per iteration)."""
    monkeypatch.setitem(THERMOMETERS, "test", Calibration(lambda c, P: 1000 + 1e-3 * P, "", {"cpx": ()}))
    monkeypatch.setitem(BAROMETERS, "test", Calibration(lambda c, T: 900 * (T - 1000), "", {"cpx": ()}))
    return make_cations({"cpx": np.ones((1, len(CATIONS)))})


def test_needs_both_T_and_P_to_converge(slow_pressure):
    # T changes by < 0.1 K from the first iteration; stopping then would leave P = 3.6
    pt = solve_pt(slow_pressure, "test", "test")
    assert pt.status[0] == "converged"
    assert pt.n_iter[0] > 30
    assert pt.P_GPa[0] < 0.1


def test_not_converged_gives_nan(slow_pressure):
    pt = solve_pt(slow_pressure, "test", "test", max_iter=5)
    assert pt.status[0] == "not converged"
    assert pt.n_iter[0] == 5
    assert np.isnan(pt.T_C[0]) and np.isnan(pt.P_GPa[0])


def test_missing_minerals():
    c = make_cations({"cpx": np.ones((2, len(CATIONS)))})
    pt = solve_pt(c, "T98", "NG85")
    assert (pt.status == "missing minerals").all()
    assert (pt.n_iter == 0).all()
    assert pt.T_C.isna().all()


def test_not_finite():
    # No Cr in cpx: log of zero CaCr-Tschermak activity in the SUD21 barometer
    cpx = pd.DataFrame([dict(zip(CATIONS, [2, 0, 0.1, 0, 0.05, 0, 0, 0.85, 0.8, 0.1, 0, 0, 0, 0]))])
    pt = solve_pt({"cpx": cpx}, "NT00", "SUD21")
    assert pt.status[0] == "not finite"
    assert np.isnan(pt.P_GPa[0])
