"""Calibrations against values published in the papers, by hand from printed formulas, or by construction."""

import numpy as np
import pandas as pd
import pytest

from xenopt import BAROMETERS, OXIDES, THERMOMETERS, cations_per_formula_unit, thermobarometers
from xenopt.chemistry import C2K
from xenopt.thermobarometers import (
    _CUBIC,
    _CUBIC_B,
    _CUBIC_C,
    P_beyer15,
    P_grutter06,
    P_macgregor74,
    T_canil99,
    T_nimis24_ni,
    T_ryan96_ni,
    T_sudholz21_ni,
    _holdaway01,
    spinel_sites,
)

from .sample import NAMES, PERIDOTITE_CPX, T_C, WT, C

NI_PER_NIO = 58.690 / 74.689


def nio_wt(ni_ppm):
    return pd.DataFrame({"NiO": np.atleast_1d(ni_ppm) / 1e4 / NI_PER_NIO})


# Grutter et al. (2006, Table 2): garnet Cr2O3, CaO [wt%] and P38 [kbar]
GRUTTER06 = [
    (13.29, 3.19, 60.0),
    (12.76, 2.75, 59.7),
    (13.0, 5.36, 52.5),
    (12.81, 5.18, 52.5),
    (10.7, 2.92, 52.5),
    (6.89, 0.25, 48.3),
    (10.14, 5.71, 42.4),
    (6.91, 2.2, 42.4),
    (9.16, 7.02, 36.8),
    (12.05, 2.19, 59.1),
    (9.72, 2.77, 49.8),
]


@pytest.mark.parametrize("cr2o3, cao, P38", GRUTTER06)
def test_grutter06_p38(cr2o3, cao, P38):
    wt = {"gt": pd.DataFrame({"Cr2O3": [cr2o3], "CaO": [cao]})}
    assert 10 * float(P_grutter06({}, 1300.0, wt)[0]) == pytest.approx(P38, abs=0.15)


@pytest.mark.parametrize(
    "cr2o3, cao, geotherm, expected",
    [(13.29, 3.19, "35", 61.1), (12.05, 2.19, "46", 56.0), (9.72, 2.77, "48", 45.9), (10.14, 5.71, "40", 41.7)],
)
def test_grutter06_geotherm_correction(cr2o3, cao, geotherm, expected):
    wt = {"gt": pd.DataFrame({"Cr2O3": [cr2o3], "CaO": [cao]})}
    assert 10 * float(P_grutter06({}, 1300.0, wt, geotherm=geotherm)[0]) == pytest.approx(expected, abs=0.15)


# Nimis et al. (2024, Table 1): garnet X_Ca, X_Cr, Ni [ppm], P [GPa] and T [C] from
# RY96, CN99, SU21 and their eq. 2, all with 2900 ppm Ni in olivine
NIMIS24 = [
    (0.176, 0.298, 64.71, 6.06, 1120, 1112, 1309, 1187),
    (0.127, 0.127, 83.90, 6.48, 1222, 1171, 1227, 1270),
    (0.138, 0.170, 41.45, 5.75, 974, 1021, 1085, 1077),
]


@pytest.mark.parametrize("x_ca, x_cr, ni, P, T_ry96, T_cn99, T_su21, T_eq2", NIMIS24)
def test_ni_in_garnet_against_nimis24(x_ca, x_cr, ni, P, T_ry96, T_cn99, T_su21, T_eq2):
    # Garnet cations with these X_Ca = Ca/(Ca+Fe+Mg) and X_Cr = Cr/(Cr+Al)
    gt = pd.DataFrame({"Ca": [3 * x_ca], "Fe": [0.3], "Mg": [3 * (1 - x_ca) - 0.3]})
    gt = gt.assign(Cr=2 * x_cr, Al=2 - 2 * x_cr)
    wt = {"gt": nio_wt(ni), "ol": nio_wt(2900)}
    c = {"gt": gt}
    assert float(T_ryan96_ni(c, P, wt)[0]) - C2K == pytest.approx(T_ry96, abs=1)
    assert float(T_canil99(c, P, wt)[0]) - C2K == pytest.approx(T_cn99, abs=1)
    assert float(T_sudholz21_ni(c, P, wt, ni_ol="measured")[0]) - C2K == pytest.approx(T_su21, abs=3)
    assert float(T_nimis24_ni(c, P, wt)[0]) - C2K == pytest.approx(T_eq2, abs=1)


def test_ni_assumed_olivine():
    wt = {"gt": nio_wt(64.71), "ol": nio_wt(2900)}
    assert float(T_nimis24_ni({}, 6.06, wt, ni_ol="measured")[0]) == pytest.approx(float(T_nimis24_ni({}, 6.06, wt)[0]))
    # Canil and Sudholz assume 3000 ppm, a few degrees cooler than 2900
    assert float(T_canil99({}, 6.06, wt, ni_ol="assumed")[0]) < float(T_canil99({}, 6.06, wt)[0])
    # NiO reported as zero (below detection) gives no temperature, not 0 K
    assert np.isnan(float(T_ryan96_ni({}, 4.0, {"gt": nio_wt(0.0)})[0]))


def test_spinel_fe3_from_stoichiometry():
    # (Mg0.7 Fe2+0.3)(Al1.0 Cr0.8 Fe3+0.2)O4 entered with all Fe as FeOt gets its Fe3+ back
    moles = {"MgO": 0.7, "FeOt": 0.3 + 0.2, "Al2O3": 0.5, "Cr2O3": 0.4}
    weights = {o.name: o.weight for o in OXIDES}
    wt = {name: 0.0 for name in NAMES}
    wt.update({name: n * weights[name] for name, n in moles.items()})
    sp = cations_per_formula_unit(pd.DataFrame([wt]), "sp")
    sites = spinel_sites(sp).iloc[0]
    assert sites.Fe3 == pytest.approx(0.2)
    assert sites.Fe2 == pytest.approx(0.3)
    assert sites.Cr == pytest.approx(0.8)


def test_olivine_spinel_al_thermometers():
    # Coogan et al. (2014) is a refit of Wan et al. (2008); they agree closely for typical spinel
    assert T_C("CO14") == pytest.approx(T_C("WA08"), abs=20)


def test_spinel_peridotite_opx_thermometer():
    # Witt-Eickschen & Seck (1991) eq. 7 by hand: 636.54 + 2088.21 Al(M1) + 14527.32 Cr
    opx = C["opx"]
    al_m1 = float((0.5 * (opx.Al - opx.Cr - 2 * opx.Ti + opx.Na))[0])
    expected = 636.54 + 2088.21 * al_m1 + 14527.32 * float(opx.Cr[0])
    assert T_C("WES91", version="Cr-Al-opx") == pytest.approx(expected)


# Beyer et al. (2015, Table 3): experiments on basaltic composition OC1/2, garnet and cpx
# oxides SiO2 TiO2 Al2O3 FeO MgO MnO CaO Na2O [wt%], at P [GPa] and T [C].
# H3483 is printed as 56 GPa; 6 GPa from Table 2.
BEYER15 = {
    "A571": (3, 1430, [40.0, 0.4, 22.3, 11.3, 12.8, 0.29, 11.0, 0.12], [49.0, 0.37, 13, 5.3, 10.8, 0.10, 17.3, 2.45]),
    "H3547": (4, 1485, [40.1, 0.35, 22.2, 13.3, 11.3, 0.31, 10.9, 0.19], [51.3, 0.32, 13.7, 5.5, 9.1, 0.09, 14.6, 3.9]),
    "H3543": (4, 1460, [39.9, 0.45, 21.8, 14.2, 10.7, 0.19, 10.8, 0.17], [51.8, 0.36, 13.5, 5.9, 8.8, 0.07, 14.4, 3.9]),
    "H3540": (6, 1480, [40.0, 0.4, 21.4, 12.2, 11.0, 0.25, 12.7, 0.28], [53.9, 0.16, 12.4, 4.6, 8.7, 0.05, 13.0, 5.6]),
    "H3639": (5, 1445, [41.1, 0.35, 22.9, 9.2, 12.8, 0.28, 13.0, 0.20], [53.8, 0.23, 12.8, 3.4, 9.9, 0.06, 15.0, 4.6]),
    "H3483": (6, 1549, [40.6, 0.34, 22.1, 11.6, 11.7, 0.24, 12.1, 0.29], [53.4, 0.19, 13.8, 4.4, 8.7, 0.06, 12.8, 5.4]),
    "H3602": (6, 1518, [40.5, 0.52, 21.5, 13.4, 10.4, 0.26, 11.8, 0.27], [53.7, 0.5, 12.7, 5.1, 8.2, 0.05, 13.2, 4.9]),
    "H3510": (6, 1380, [39.7, 0.6, 21.0, 14.3, 9.6, 0.25, 12.3, 0.28], [54.5, 0.35, 12.9, 5.0, 8.2, 0.06, 12.9, 5.2]),
    "H3605x": (7, 1480, [40.1, 0.52, 21.6, 14.1, 9.9, 0.25, 12.0, 0.29], [54.9, 0.30, 12.8, 5.0, 8.2, 0.06, 12.8, 5.0]),
}
BEYER15_OXIDES = ["SiO2", "TiO2", "Al2O3", "FeOt", "MgO", "MnO", "CaO", "Na2O"]


@pytest.mark.parametrize("run", BEYER15)
def test_beyer15_reproduces_experiments(run):
    # Their calibration reproduces experimental P to about 0.2 GPa on average
    P, T, gt, cpx = BEYER15[run]
    c = {
        m: cations_per_formula_unit(pd.DataFrame([dict(zip(BEYER15_OXIDES, v))]).reindex(columns=NAMES).fillna(0), m)
        for m, v in [("gt", gt), ("cpx", cpx)]
    }
    assert float(P_beyer15(c, T + C2K)[0]) == pytest.approx(P, abs=0.35)


# Holdaway (2001) grossular polynomials (Wu & Zhao 2007, Table S2, eqs. A1-A3), only to test
# the almandine and pyrope ones used by the thermometer
GROSSULAR = {
    "Caa": [
        (0.337, "Fe Fe"), (-18.983, "Mg Mg"), (-9.5725, "Fe Mg"), (-19.34, "Fe Ca"), (-4.6665, "Fe Mn"),
        (-25.746, "Mg Ca"), (-38.938, "Mg Mn"), *_CUBIC,
    ],
    "Cab": [
        (0.04, "Fe Fe"), (0.068, "Mg Mg"), (0.065, "Mn Mn"), (-0.0515, "Fe Mg"), (0.27, "Fe Ca"), (0.1765, "Fe Mn"),
        (0.28, "Mg Ca"), (0.207, "Mg Mn"), (0.13, "Ca Mn"), *_CUBIC_B,
    ],
    "Cac": [
        (-1304.0, "Fe Fe"), (66114.0, "Mg Mg"), (1425.0, "Mn Mn"), (31326.5, "Fe Mg"), (39864.0, "Fe Ca"),
        (12356.0, "Fe Mn"), (51518.0, "Mg Ca"), (88610.5, "Mg Mn"), (2850.0, "Ca Mn"), *_CUBIC_C,
    ],
}  # fmt: skip


@pytest.mark.parametrize("part", ["a", "b", "c"])
def test_holdaway01_gibbs_duhem(part, monkeypatch):
    # In a Ca-Fe-Mg garnet, sum of X_i d(RT ln gamma_i) = 0 for any change in composition
    monkeypatch.setattr(thermobarometers, "HOLDAWAY01", {**thermobarometers.HOLDAWAY01, **GROSSULAR})
    rng = np.random.default_rng(0)

    def rt_ln_gamma(cation, x):
        return float(_holdaway01(cation + part, {k: np.array([v]) for k, v in x.items()})[0])

    for _ in range(50):
        x = dict(zip(["Fe", "Mg", "Ca"], rng.dirichlet([2, 2, 1])), Mn=0.0)
        d = rng.normal(size=3)
        d = dict(zip(["Fe", "Mg", "Ca"], d - d.mean()), Mn=0.0)
        eps = 1e-6
        plus = {k: x[k] + eps * d[k] for k in x}
        minus = {k: x[k] - eps * d[k] for k in x}
        residual = sum(x[i] * (rt_ln_gamma(i, plus) - rt_ln_gamma(i, minus)) / (2 * eps) for i in ["Fe", "Mg", "Ca"])
        scale = max(abs(rt_ln_gamma(i, x)) for i in ["Fe", "Mg", "Ca"])
        assert abs(residual) < 1e-7 * scale


def test_wu_zhao07_close_to_oneill_wood79():
    # Same exchange; Nimis & Grutter (2010) find WZ07 higher than T98 below 1200 C
    assert T_C("WZ07") == pytest.approx(T_C("OW79"), abs=100)


def test_ashchepkov17_thermometer_correction():
    # Their correction worked by hand from the cpx (wt%) and NT00 at 4 GPa
    T = T_C("NT00") + C2K
    cpx = WT["cpx"].iloc[0]
    expected = (
        -1e-6 * T**2 + 0.95 * T + 33.75 * (cpx.Al2O3 + cpx.Cr2O3) - 55 - 17.9 * (cpx.FeOt - 2.5) - 10.9 * cpx.Na2O - C2K
    )
    assert T_C("AS17") == pytest.approx(expected)


def test_ashchepkov17_close_to_nimis_taylor00():
    # The paper reports precision comparable to Nimis & Taylor (2000), for peridotite cpx
    wt = {"cpx": pd.DataFrame([PERIDOTITE_CPX])}
    c = {"cpx": cations_per_formula_unit(wt["cpx"], "cpx")}
    T = 1000 + C2K
    nt00 = 10 * float(BAROMETERS["NT00"].func(c, T)[0])
    assert 10 * float(BAROMETERS["AS17"].func(c, T, wt)[0]) == pytest.approx(nt00, abs=3)
    # The printed MgO/CaO (wt%) instead gives ~8 kbar less here
    assert 10 * float(BAROMETERS["AS17"].func(c, T, wt, kd="wt")[0]) < nt00 - 5


def test_brey_kohler90_fe_mg_follows_table_1():
    ol, cpx, gt, opx = C["ol"], C["cpx"], C["gt"], C["opx"]
    ln_kd = float(np.log(ol.Fe * cpx.Mg / (ol.Mg * cpx.Fe))[0])
    assert T_C("BK90olcpx") == pytest.approx((1489 + 4.64 * 40) / (ln_kd + 1.07) - C2K)
    ln_kd = float(np.log(gt.Fe * opx.Mg / (gt.Mg * opx.Fe))[0])
    assert T_C("BK90gtopx") == pytest.approx((1456 + 9.86 * 40) / (ln_kd + 0.55) - C2K)


def test_ryan96_zn_published_coefficients():
    zn_ppm = 0.23e4 * 65.38 / 81.379
    assert T_C("RY96Zn") == pytest.approx(1000 / (-0.90 + 0.26 * np.log(zn_ppm)) - C2K)
    assert T_C("RY96Zn", coefficients="PTEXL") < T_C("RY96Zn")


def test_options_from_other_implementations_are_close():
    # Within a few tens of degrees of the published versions
    assert T_C("WE77", sites="PTEXL") == pytest.approx(T_C("WE77"), abs=15)
    assert T_C("BM85", opx_fe="PTEXL") == pytest.approx(T_C("BM85"), abs=5)
    # Creighton (2009) as printed is hotter
    assert T_C("CR09", cr_term="printed") > T_C("CR09") + 50


def test_oneill_wood79_rises_with_pressure():
    # The O'Neill (1980) correction gives T rising with P, as the reaction volume requires
    assert T_C("OW79", P=5.0) > T_C("OW79", P=3.0)


def test_iterative_thermometers_converge():
    for key in ["BM85", "OW79"]:
        cal = THERMOMETERS[key]
        T = cal.func(C, 4.0, **cal.resolve())
        assert np.isfinite(np.asarray(T)).all()


# MacGregor (1974, Table 4): opx Al2O3 [wt%] with garnet + forsterite, at T [C] and P [kbar]
MACGREGOR = [(1000, 30, 2.4), (1100, 30, 3.3), (1100, 40, 1.7), (1200, 30, 4.9), (1200, 35, 3.1), (1300, 35, 4.9)]


@pytest.mark.parametrize("T, P, al2o3", MACGREGOR)
def test_macgregor_fit_reproduces_experiments(T, P, al2o3):
    wt = {"opx": pd.DataFrame({"Al2O3": [al2o3]})}
    assert 10 * float(P_macgregor74({}, T + C2K, wt)[0]) == pytest.approx(P, abs=3)
