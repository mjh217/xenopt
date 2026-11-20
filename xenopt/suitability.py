"""Tests of whether analyses and P-T results are suitable for thermobarometry.

cpx_tests() screens cpx for single-cpx thermobarometry, and opx_garnet_tests()
checks opx-cpx-garnet thermobarometry. pair_tests() picks the tests that apply
to a thermometer-barometer pair and warning_text() lists those a sample failed.

Each test is True when the analysis passes. Any test involving a missing
analysis is False. Thresholds default to the published values and can be
changed with keyword arguments.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from .chemistry import C2K
from .thermobarometers import (
    Cations,
    T_brey_kohler90_ca_in_opx,
    T_krogh88,
    T_nimis_grutter10,
    T_taylor98,
    cacrts_activity,
    cpx_cr_number,
)


def cpx_tests(
    wt_perc: pd.DataFrame,
    cpx: pd.DataFrame,
    cation_range: tuple[float, float] = (3.98, 4.02),
    cacrts_cr_min: float = 0.005,
    ca_mg_max: float = 0.5,
    cr2o3_min: float = 0.5,
    cr2o3_min_inclusive: bool = False,
) -> pd.DataFrame:
    """Run the cpx tests given cpx oxide wt% and cations per formula unit.

    The temperature test depends on the P-T solution and is evaluated
    separately with temperature_test().

    cation_range: allowed total cations per 6 O. 3.98-4.02 from Nimis & Taylor
        (2000) and Ziberna et al. (2016); Grutter (2009) used 3.96-4.04 and
        Hoggard et al. (2026) 3.97-4.03.
    cacrts_cr_min: lower limit of a_CaCrTs / Cr#. 0.005 from Sudholz et al.
        (2021); Ziberna et al. (2016) discard values <= 0.011.
    ca_mg_max: upper limit of Ca / (Ca + Mg), above which equilibrium with opx
        is suspect (Ziberna et al., 2016).
    cr2o3_min, cr2o3_min_inclusive: lower limit of cpx Cr2O3 [wt%] in the
        Cr2O3 vs Al2O3 test, Cr2O3 > 0.5 (Ziberna et al., 2016);
        cr2o3_min_inclusive=True also accepts Cr2O3 equal to it.
    """
    wt_perc = wt_perc.fillna(0)
    total = wt_perc.sum(axis=1)
    cr2o3, al2o3, mgo = wt_perc.Cr2O3, wt_perc.Al2O3, wt_perc.MgO
    cr_no = cpx_cr_number(cpx)
    cation_total = cpx.sum(axis=1, skipna=False)

    tests = pd.DataFrame(index=wt_perc.index)
    tests["cpx_present"] = total > 0
    # Data quality (Ziberna et al., 2016, step 1)
    tests["cpx_oxide_total"] = (total >= 98) & (total <= 102)
    tests["cpx_cation_total"] = (cation_total >= cation_range[0]) & (cation_total <= cation_range[1])

    # Garnet-peridotite field of the Cr2O3 vs Al2O3 diagram of Ramsay (1995), as
    # used by Nimis (1998) and Ziberna et al. (2016, step 3), and the high-Al
    # field of Al2O3 vs MgO (Nimis, 1998, Fig. 3; Ziberna et al., 2016, step 4)
    above_min = cr2o3 >= cr2o3_min if cr2o3_min_inclusive else cr2o3 > cr2o3_min
    tests["nimis98_cr_al"] = above_min & (((cr2o3 < 2.25) & (al2o3 <= 4)) | ((cr2o3 >= 2.25) & (al2o3 <= 5)))
    tests["nimis98_al_mg"] = (al2o3 >= 0.7) & (al2o3 >= 12.175 - 0.6375 * mgo)

    # Sudholz et al. (2021): not eclogitic (Na <= 0.2), Cr# within 0.10-0.65,
    # and a_CaCrTs / Cr# high enough to be insensitive to analytical errors
    tests["sudholz21_na"] = cpx.Na <= 0.2
    tests["sudholz21_cr_no"] = (cr_no >= 0.1) & (cr_no <= 0.65)
    tests["sudholz21_cacrts"] = cacrts_activity(cpx) / cr_no >= cacrts_cr_min

    # Ziberna et al. (2016): Cr2O3 above 5 wt% is outside the tested range
    # (step 3), and high Ca / (Ca + Mg) suggests disequilibrium with opx (step 7)
    tests["ziberna16_cr2o3"] = (cr2o3 > 0) & (cr2o3 <= 5)
    tests["ziberna16_ca_mg"] = cpx.Ca / (cpx.Ca + cpx.Mg) <= ca_mg_max

    # PTEXL garnet-peridotite test, Na + K >= Al + Cr - 0.08, a form of the
    # Al + Cr - Na discriminant of Grutter (2009)
    tests["ptexl_na_k"] = cpx.Na + cpx.K >= cpx.Al + cpx.Cr - 0.08

    return tests


def temperature_test(T_C: pd.Series, T_min: float = 600) -> pd.Series:
    """Ziberna et al. (2016, step 7): T below ~600 C suggests disequilibrium with opx.

    Nimis & Grutter (2010) suggest 700 C as a safer cut-off.
    """
    return T_C >= T_min


def opx_garnet_tests(
    wt_perc: Mapping[str, pd.DataFrame],
    c: Cations,
    P_GPa: pd.Series,
    P_min: float = 1.5,
    P_min_pyroxenite: float = 1.0,
    T_ta98_min: float = 700,
    cpx_opx_tol: tuple[float, float, float] = (90, 70, 50),
    gt_opx_tol: float = 70,
    gt_cpx_tol: float | None = None,
) -> pd.DataFrame:
    """Checks of Nimis & Grutter (2010) for opx-cpx-garnet thermobarometry.

    Temperatures are T98 at the given pressure, as in their scheme (T98 with
    NG85). Thresholds:

    P_min, P_min_pyroxenite: minimum P [GPa] for garnet to be stable in
        lherzolite and in pyroxenite. Lithology is often not recorded, so both
        are tested for every sample; warning_text() notes when only the
        lherzolite limit failed.
    T_ta98_min: minimum T [C] for reliable pyroxene thermometry.
    cpx_opx_tol: maximum |T_T98 - T_Ca-in-opx| [C] at T_T98 < 900, 900-1200
        and > 1200 C. PTEXL widens these by 10 C.
    gt_opx_tol: maximum |T_T98 - T_NG09| [C].
    gt_cpx_tol: if given, also check |T_T98 - T_Krogh88| [C] against it, as
        PTEXL does with 90. Off by default: Nimis & Grutter (2010) found
        garnet-cpx thermometers too imprecise for this (> 200 C off T98 for
        well-equilibrated samples).
    """
    index = c["opx"].index
    P = np.asarray(P_GPa, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        T_ta98 = pd.Series(np.asarray(T_taylor98(c, P), dtype=float) - C2K, index=index)
        T_ca_in_opx = np.asarray(T_brey_kohler90_ca_in_opx(c, P, "NG10"), dtype=float) - C2K
        T_ng09 = np.asarray(T_nimis_grutter10(c, P), dtype=float) - C2K
        T_kr88 = np.asarray(T_krogh88(c, P), dtype=float) - C2K

    tests = pd.DataFrame(index=index)
    # Analysis quality: oxide totals for all minerals and cation sums for pyroxenes
    for mineral in ("opx", "cpx", "gt"):
        total = wt_perc[mineral].sum(axis=1)
        tests[f"ng10_{mineral}_quality"] = (total >= 98.5) & (total <= 101.5)
        if mineral != "gt":
            tests[f"ng10_{mineral}_quality"] &= c[mineral].sum(axis=1, skipna=False) >= 3.99
    # Peridotitic or pyroxenitic garnet: Cr2O3 > 1 wt%, or Mg/(Mg+Fe) > 0.6 (Grutter et al., 2004)
    gt = c["gt"]
    tests["ng10_gt_type"] = (wt_perc["gt"].Cr2O3.fillna(0) > 1) | (gt.Mg / (gt.Mg + gt.Fe) > 0.6)
    # Plausibility of the result: pyroxene thermometry is unreliable at very low
    # T, and garnet is not stable at low P
    tests["ng10_T_min"] = T_ta98 >= T_ta98_min
    tests["ng10_P_min"] = P >= P_min
    tests["ng10_P_min_pyroxenite"] = P >= P_min_pyroxenite
    # Pyroxene equilibrium: T98 (mainly cpx) against the corrected Ca-in-opx
    # thermometer (opx only); about twice the standard error in each T range
    tol = np.select([T_ta98 < 900, T_ta98 <= 1200], cpx_opx_tol[:2], cpx_opx_tol[2])
    tests["ng10_cpx_opx"] = np.abs(T_ta98 - T_ca_in_opx) <= tol
    # Fe-Mg vs Ca-Mg equilibrium: T98 against the opx-garnet Fe-Mg thermometer
    tests["ng10_gt_opx"] = np.abs(T_ta98 - T_ng09) <= gt_opx_tol
    if gt_cpx_tol is not None:
        tests["ptexl_gt_cpx"] = np.abs(T_ta98 - T_kr88) <= gt_cpx_tol
    # Accuracy is not assured below 900 C or above 60 kbar
    tests["ng10_T_accuracy"] = T_ta98 >= 900
    tests["ng10_P_accuracy"] = P <= 6.0
    return tests


def pair_tests(
    minerals: set[str],
    wt_perc: Mapping[str, pd.DataFrame],
    c: Cations,
    pt: pd.DataFrame,
    cpx: pd.DataFrame,
    T_min: float = 600,
    **opx_garnet_thresholds,
) -> pd.DataFrame:
    """The tests that apply to a P-T pair using these minerals.

    Single-cpx pairs get the cpx tests (`cpx`, from cpx_tests) and the minimum
    T test; pairs using opx, cpx and garnet get opx_garnet_tests. Other pairs
    get no tests.
    """
    if set(minerals) == {"cpx"}:
        tests = cpx.drop(columns="cpx_present")
        tests["ziberna16_T"] = temperature_test(pt["T_C"], T_min)
        return tests
    if {"opx", "cpx", "gt"} <= set(minerals):
        return opx_garnet_tests(wt_perc, c, pt["P_GPa"], **opx_garnet_thresholds)
    return pd.DataFrame(index=pt.index)


# Failed test -> what the warning says. Tests not listed appear under their own name.
LABELS = {
    # cpx tests
    "cpx_oxide_total": "cpx oxide total out of range",
    "cpx_cation_total": "cpx cation total out of range",
    "nimis98_cr_al": "outside garnet-peridotite field",
    "nimis98_al_mg": "low Al2O3 for MgO",
    "sudholz21_na": "high Na",
    "sudholz21_cr_no": "Cr# out of range",
    "sudholz21_cacrts": "low CaCrTs activity",
    "ziberna16_cr2o3": "Cr2O3 out of range",
    "ziberna16_ca_mg": "high Ca/(Ca+Mg)",
    "ptexl_na_k": "low Na+K",
    "ziberna16_T": "low T",
    # opx-garnet tests
    "ng10_opx_quality": "poor opx analysis",
    "ng10_cpx_quality": "poor cpx analysis",
    "ng10_gt_quality": "poor garnet analysis",
    "ng10_gt_type": "low-Cr, low-Mg garnet",
    "ng10_T_min": "very low T",
    "ng10_P_min": "low P",
    "ng10_P_min_pyroxenite": "very low P",
    "ng10_cpx_opx": "pyroxene temperatures disagree",
    "ng10_gt_opx": "garnet-opx temperature disagrees",
    "ptexl_gt_cpx": "garnet-cpx temperature disagrees",
    "ng10_T_accuracy": "low T (accuracy not assured)",
    "ng10_P_accuracy": "high P (accuracy not assured)",
}

# Failed test -> (companion test, note added to the warning when the companion passed)
CAVEATS = {"ng10_P_min": ("ng10_P_min_pyroxenite", "ok if pyroxenite")}

# Failed test -> milder test it makes redundant, left out of the warning when both fail
SUPERSEDES = {"ng10_P_min_pyroxenite": "ng10_P_min", "ng10_T_min": "ng10_T_accuracy"}


def warning_text(tests: pd.DataFrame, pt: pd.DataFrame) -> pd.Series:
    """Why each result might not be trustworthy: failed tests and blank oxides taken as zero.

    Failed tests are described with the short phrases in LABELS, separated by
    "; ". Samples without a result list any missing oxides instead.
    """
    converged = pt["status"] == "converged"
    columns = list(tests.columns)
    passed = tests.to_numpy(dtype=bool)
    shown = ~passed
    for name, milder in SUPERSEDES.items():
        if name in columns and milder in columns:
            shown[:, columns.index(milder)] &= passed[:, columns.index(name)]
    labels = np.tile(np.array([LABELS.get(n, n) for n in columns], dtype=object), (len(tests), 1))
    for name, (companion, note) in CAVEATS.items():
        if name in columns and companion in columns:
            labels[passed[:, columns.index(companion)], columns.index(name)] = f"{LABELS.get(name, name)} ({note})"
    text = pd.Series(["; ".join(row_labels[row]) for row_labels, row in zip(labels, shown)], index=pt.index)
    if "assumed_zero" in pt:
        assumed = pt["assumed_zero"].map(lambda s: f"assumed zero: {s}" if s else "")
        text = (text + "; " + assumed).str.strip("; ")
    text = text.where(converged, "")
    if "missing_oxides" in pt:
        text = text.mask(pt["missing_oxides"] != "", "missing oxides: " + pt["missing_oxides"])
    return text
