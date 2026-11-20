import inspect

import numpy as np
import pandas as pd
import pytest

from xenopt import BAROMETERS, THERMOMETERS
from xenopt.chemistry import CATIONS
from xenopt.thermobarometers import OPX_MGTS, P_nickel_green85

from .sample import C


def opx(**cations):
    row = dict.fromkeys(CATIONS, 0.0)
    row.update(cations)
    return pd.DataFrame([row])


# A low-Na opx (Na < Cr + Ti) and a sodic one (Na > Cr + 2Ti)
LOW_NA = opx(Al=0.05, Cr=0.01, Ti=0.005, Na=0.004)
SODIC = opx(Al=0.05, Cr=0.01, Ti=0.005, Na=0.03)


@pytest.mark.parametrize(
    "scheme, low_na, sodic",
    [
        # Nickel & Green (1985): (Al - Cr - 2Ti + Na) / 2 throughout
        ("NG85", (0.05 - 0.01 - 0.01 + 0.004) / 2, (0.05 - 0.01 - 0.01 + 0.03) / 2),
        # Carswell (1991) E9 when Na > Cr + Ti: (Al + Cr - Na) / 2
        ("CA91", (0.05 - 0.01 - 0.01 + 0.004) / 2, (0.05 + 0.01 - 0.03) / 2),
        # Carswell & Gibb (1987b) when Na > Cr + 2Ti: (Al + Cr + 2Ti - Na) / 2
        ("CG87", (0.05 - 0.01 - 0.01 + 0.004) / 2, (0.05 + 0.01 + 0.01 - 0.03) / 2),
        # Taylor (1998): (Al - Cr + Ti + Na) / 2 throughout
        ("T98", (0.05 - 0.01 + 0.005 + 0.004) / 2, (0.05 - 0.01 + 0.005 + 0.03) / 2),
    ],
)
def test_mgts_schemes(scheme, low_na, sodic):
    assert OPX_MGTS[scheme](LOW_NA)[0] == pytest.approx(low_na)
    assert OPX_MGTS[scheme](SODIC)[0] == pytest.approx(sodic)


@pytest.mark.parametrize("scheme, threshold", [("CA91", 0.01 + 0.005), ("CG87", 0.01 + 2 * 0.005)])
def test_jadeite_schemes_are_continuous(scheme, threshold):
    """Both branches agree at the Na threshold, so P does not jump there."""
    below, above = opx(Al=0.05, Cr=0.01, Ti=0.005, Na=threshold - 1e-9), opx(Al=0.05, Cr=0.01, Ti=0.005, Na=threshold)
    assert OPX_MGTS[scheme](below)[0] == pytest.approx(OPX_MGTS[scheme](above)[0], abs=1e-8)


def test_mgts_in_only_matters_for_modified_schemes():
    T = 1273.15
    ng85 = [float(P_nickel_green85(C, T, "NG85", scope)[0]) for scope in ("activity", "activity_dV", "all")]
    assert np.ptp(ng85) == 0
    ta98 = [float(P_nickel_green85(C, T, "T98", scope)[0]) for scope in ("activity", "activity_dV", "all")]
    assert len(set(ta98)) == 3


def test_options():
    assert BAROMETERS["NG85"].resolve() == {"mgts": "NG85", "mgts_in": "activity_dV"}
    assert BAROMETERS["T98"].resolve() == {"mgts": "T98", "mgts_in": "activity_dV"}
    assert BAROMETERS["T98"].resolve({"mgts_in": "all"})["mgts_in"] == "all"
    with pytest.raises(ValueError, match="must be one of"):
        BAROMETERS["NG85"].resolve({"mgts": "XX"})
    with pytest.raises(ValueError, match="unknown option"):
        BAROMETERS["SUD21"].resolve({"mgts": "CA91"})


@pytest.mark.parametrize("key, cal", [*THERMOMETERS.items(), *BAROMETERS.items()])
def test_calibration_cites_source(key, cal):
    # Every calibration names its paper and the equation it implements
    source = inspect.getdoc(cal.func).split("Source: ")[-1]
    assert (
        "Source: " in inspect.getdoc(cal.func)
        and "(" in source
        and any(word in source for word in ("eq", "Table", "p. "))
    )
