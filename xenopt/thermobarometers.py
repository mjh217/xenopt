"""Thermometers and barometers.

Every calibration takes `c`, a mapping of mineral ("opx", "cpx", "gt", ...) to a
DataFrame of cations per formula unit (see chemistry.cations_per_formula_unit),
plus the other intensive variable. Thermometers take pressure in GPa and return
temperature in K; barometers take temperature in K and return pressure in GPa.

Defaults follow the original publications. Later modifications are selected
with keyword options, listed for each calibration in the registry.

To add a calibration, write a function with that signature and register it in
THERMOMETERS or BAROMETERS at the bottom of this module.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import pandas as pd

from .chemistry import C2K, element_ppm

Cations = Mapping[str, pd.DataFrame]


# --------------------------------------------------------------------------
# Site occupancies and activities shared between calibrations and tests
# --------------------------------------------------------------------------


def cpx_cr_number(cpx: pd.DataFrame) -> pd.Series:
    """Cr / (Cr + Al) in clinopyroxene."""
    return cpx.Cr / (cpx.Cr + cpx.Al)


def cacrts_activity(cpx: pd.DataFrame) -> pd.Series:
    """Activity of the CaCr-Tschermak component in clinopyroxene (Nimis & Taylor, 2000)."""
    return cpx.Cr - 0.81 * cpx_cr_number(cpx) * (cpx.Na + cpx.K)


def _enstatite_activity_ta98(px: pd.DataFrame) -> pd.Series:
    """Enstatite activity in a pyroxene, as defined by Taylor (1998)."""
    al_vi = 0.5 * (px.Al - px.Cr - 2 * px.Ti + px.Na)
    al_iv = 0.5 * (px.Al + px.Cr + 2 * px.Ti - px.Na)
    return (1 - px.Ca - px.Na) * (1 - al_vi - px.Cr - px.Ti) * (1 - 0.5 * al_iv) ** 2


def _opx_al_m1(opx: pd.DataFrame) -> pd.Series:
    """Al on the opx M1 site from charge balance, (Al - Cr - 2Ti + Na) / 2 (Nickel & Green, 1985)."""
    return 0.5 * (opx.Al - opx.Cr - 2 * opx.Ti + opx.Na)


# Mole fraction of Mg-Tschermak (MgAl2SiO6) in opx, as used in the NG85
# barometer. Fe3+ is taken as zero throughout.


def _mgts_carswell91(opx: pd.DataFrame) -> pd.Series:
    """Carswell (1991) eqs. E6 and E9: Na in excess of Cr + Ti forms jadeite, giving (Al + Cr - Na) / 2."""
    return _opx_al_m1(opx).where(opx.Na <= opx.Cr + opx.Ti, 0.5 * (opx.Al + opx.Cr - opx.Na))


def _mgts_carswell_gibb87(opx: pd.DataFrame) -> pd.Series:
    """Carswell & Gibb (1987b), eq. E7 of Carswell (1991), as used in PTEXL.

    Na in excess of Cr + 2Ti forms jadeite, giving (Al + Cr + 2Ti - Na) / 2.
    """
    jadeite = opx.Na - opx.Cr - 2 * opx.Ti
    return _opx_al_m1(opx).where(jadeite <= 0, 0.5 * (opx.Al - jadeite))


def _mgts_taylor98(opx: pd.DataFrame) -> pd.Series:
    """Taylor (1998) empirical Ti correction, (Al - Cr + Ti + Na) / 2."""
    return 0.5 * (opx.Al - opx.Cr + opx.Ti + opx.Na)


OPX_MGTS = {
    "NG85": _opx_al_m1,
    "CA91": _mgts_carswell91,
    "CG87": _mgts_carswell_gibb87,
    "T98": _mgts_taylor98,
}


# --------------------------------------------------------------------------
# Thermometers: f(c, P [GPa]) -> T [K]
# --------------------------------------------------------------------------


def T_taylor98(c: Cations, P):
    """Taylor (1998) two-pyroxene thermometer (Ca-Mg exchange between cpx and opx), "TA97" in the paper.

    Source: Taylor (1998), N. Jb. Miner. Abh. 172, 381-408: unnumbered equation on p. 403, with a(En) on p. 402.
    """
    opx, cpx = c["opx"], c["cpx"]
    x_ts = cpx.Al + cpx.Cr - cpx.Na
    ln_kd = np.log(_enstatite_activity_ta98(cpx) / _enstatite_activity_ta98(opx))
    return (24787 + 678 * P) / (15.67 + 14.37 * cpx.Ti + 3.69 * cpx.Fe - 3.25 * x_ts + ln_kd**2)


def T_nimis_taylor00(c: Cations, P):
    """Nimis & Taylor (2000) enstatite-in-cpx thermometer.

    Source: Nimis & Taylor (2000), Contrib. Mineral. Petrol. 139, 541-554: eq. 17, with a(En) from eq. 15.
    """
    cpx = c["cpx"]
    a_en = (1 - cpx.Ca - cpx.Na - cpx.K) * (1 - 0.5 * (cpx.Al + cpx.Cr + cpx.Na + cpx.K))
    return (23166 + 392.8 * P) / (
        13.25 + 15.35 * cpx.Ti + 4.50 * cpx.Fe - 1.55 * (cpx.Al + cpx.Cr - cpx.Na - cpx.K) + np.log(a_en) ** 2
    )


def T_ashchepkov17(c: Cations, P, wt):
    """Nimis & Taylor (2000) enstatite-in-cpx thermometer with the correction of Ashchepkov et al. (2017, p. 783).

    The correction uses cpx oxide wt% (FeOt as FeO). Their further step "for Cr
    low Cpx", (T + 390) P / (P - 70) 0.95, is negative at any mantle pressure as
    printed and is not applied.

    Source: Ashchepkov et al. (2017), Geosci. Front. 8, 775-795: unnumbered equations on p. 783, applied to
    Nimis & Taylor (2000), Contrib. Mineral. Petrol. 139, 541-554, eq. 17.
    """
    T = T_nimis_taylor00(c, P)
    cpx = wt["cpx"]
    return -1e-6 * T**2 + 0.95 * T + 33.75 * (cpx.Al2O3 + cpx.Cr2O3) - 55 - 17.9 * (cpx.FeOt - 2.5) - 10.9 * cpx.Na2O


def T_brey_kohler90_ca_in_opx(c: Cations, P, correction: str = "none"):
    """Brey & Kohler (1990) Ca-in-opx thermometer (their eq. 10).

    correction="NG10" applies the empirical correction of Nimis & Grutter (2010,
    eq. 1), which makes it consistent with the T98 thermometer at all T.

    Source: Brey & Kohler (1990), J. Petrol. 31, 1353-1378: eq. 10. correction="NG10": Nimis & Grutter (2010),
    Contrib. Mineral. Petrol. 159, 411-427, eq. 1.
    """
    T = (6425 + 26.4 * 10 * P) / (-np.log(c["opx"].Ca) + 1.843)
    if correction == "NG10":
        T_C = T - C2K
        T = -628.7 + 2.0690 * T_C - 4.530e-4 * T_C**2 + C2K
    return T


def T_nimis_grutter10(c: Cations, P):
    """Nimis & Grutter (2010) opx-garnet Fe-Mg exchange thermometer (eq. 4, "NG09" in the paper).

    Source: Nimis & Grutter (2010), Contrib. Mineral. Petrol. 159, 411-427: eq. 4.
    """
    opx, gt = c["opx"], c["gt"]
    x_site = gt.Ca + gt.Fe + gt.Mg + gt.Mn
    ln_kd = np.log(gt.Fe * opx.Mg / (gt.Mg * opx.Fe))
    return (1215 + 17.4 * 10 * P + 1495 * (gt.Ca + gt.Mn) / x_site) / (ln_kd + 0.732)


def T_krogh88(c: Cations, P):
    """Krogh (1988) garnet-cpx Fe-Mg exchange thermometer, with all Fe as Fe2+.

    Source: Krogh (1988), Contrib. Mineral. Petrol. 99, 44-48: eq. 14.
    """
    cpx, gt = c["cpx"], c["gt"]
    x_ca = gt.Ca / (gt.Ca + gt.Fe + gt.Mg + gt.Mn)
    ln_kd = np.log(gt.Fe * cpx.Mg / (gt.Mg * cpx.Fe))
    return (-6173 * x_ca**2 + 6731 * x_ca + 1879 + 10 * 10 * P) / (ln_kd + 1.393)


def _x_gt(gt: pd.DataFrame, cation: str) -> pd.Series:
    """Fraction of a garnet X-site cation, e.g. X_Ca = Ca / (Ca + Mg + Fe + Mn)."""
    return gt[cation] / (gt.Ca + gt.Mg + gt.Fe + gt.Mn)


def _ln_kd_fe_mg(a: pd.DataFrame, b: pd.DataFrame) -> pd.Series:
    """ln K_D = ln[(Fe/Mg)_a / (Fe/Mg)_b], with all Fe as Fe2+."""
    return np.log(a.Fe * b.Mg / (a.Mg * b.Fe))


def _x_fe(m: pd.DataFrame) -> pd.Series:
    return m.Fe / (m.Fe + m.Mg)


def _px_site_fractions(px: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Fractions of the pyroxene M1 and M2 sites filled by Mg + Fe (Nickel & Green, 1985).

    M1 holds Al (from charge balance), Cr and Ti, and M2 holds Ca, Na and Mn.
    """
    return 1 - _opx_al_m1(px) - px.Cr - px.Ti, 1 - px.Ca - px.Na - px.Mn


def _solve_T(f: Callable, n: int, lo: float = 300.0, hi: float = 3000.0, n_iter: int = 60) -> np.ndarray:
    """Solve T = f(T) [K] by bisection, for thermometers with terms that depend on T.

    f must decrease with T, so there is one root; where f is undefined (NaN),
    T is taken as too high. Samples without a root in [lo, hi] give NaN.
    """
    lo_T, hi_T = np.full(n, lo), np.full(n, hi)
    with np.errstate(divide="ignore", invalid="ignore"):
        for _ in range(n_iter):
            mid = 0.5 * (lo_T + hi_T)
            below = ~(np.asarray(f(mid), dtype=float) > mid)
            hi_T, lo_T = np.where(below, mid, hi_T), np.where(below, lo_T, mid)
        T = 0.5 * (lo_T + hi_T)
        found = np.abs(np.asarray(f(T), dtype=float) - T) < 1e-3
    return np.where(found, T, np.nan)


def T_brey_kohler90(c: Cations, P):
    """Brey & Kohler (1990) two-pyroxene thermometer, "T_BKN" (their eq. 9).

    Source: Brey & Kohler (1990), J. Petrol. 31, 1353-1378: eq. 9.
    """
    opx, cpx = c["opx"], c["cpx"]

    def ca_star(px):  # Ca on M2 corrected for Na (their eq. 7)
        return px.Ca / (1 - px.Na)

    ln_kd = np.log((1 - ca_star(cpx)) / (1 - ca_star(opx)))
    return (23664 + (24.9 + 126.3 * _x_fe(cpx)) * 10 * P) / (13.38 + ln_kd**2 + 11.59 * _x_fe(opx))


def T_brey_kohler90_na(c: Cations, P):
    """Brey & Kohler (1990) Na partitioning between opx and cpx (their eq. 11).

    Source: Brey & Kohler (1990), J. Petrol. 31, 1353-1378: eq. 11.
    """
    ln_d = np.log(c["opx"].Na / c["cpx"].Na)
    return (35000 + 61.5 * 10 * P) / (ln_d**2 + 19.8)


def T_bertrand_mercier85(c: Cations, P, opx_fe: str = "none"):
    """Bertrand & Mercier (1985) two-pyroxene thermometer.

    Ca on M2 of both pyroxenes is corrected for Na, and in cpx also for Fe with
    a T-dependent term, so T is solved for numerically. The Fe correction is
    calibrated for Fe/(Fe+Mg) < 0.5, and can push cpx Ca above 1 (no result)
    for Fe-rich cpx. opx_fe="PTEXL" adds the
    0.075 Fe/(Fe+Mg) correction to opx Ca used in PTEXL; it is not in the
    paper, which found the Fe effect in opx negligible.

    Source: Bertrand & Mercier (1985/86), Earth Planet. Sci. Lett. 76, 109-122: eq. 5 with the coefficients of
    Table 3 (as rounded in the abstract), M2 Ca corrected for Na (eq. 7) and, in cpx, Fe (eq. 8).
    """
    opx, cpx = c["opx"], c["cpx"]
    ca_opx = opx.Ca / (1 - opx.Na)
    if opx_fe == "PTEXL":
        ca_opx = ca_opx + 0.075 * _x_fe(opx)

    def f(T):
        ca_cpx = cpx.Ca / (1 - cpx.Na) + (-0.77 + T / 1000) * _x_fe(cpx)
        ln_k = np.log((1 - ca_cpx) / (1 - ca_opx))
        return (36273 + 399 * P) / (19.31 - 8.314 * ln_k - 12.15 * ca_cpx**2)

    return _solve_T(f, len(cpx))


def T_wells77(c: Cations, P, sites: str = "WB73"):
    """Wells (1977) two-pyroxene thermometer (his eq. 5).

    Enstatite activities are X_Mg^M1 * X_Mg^M2 from the ideal two-site model of
    Wood & Banno (1973, eqs. 24-25): Ca, Na and Mn fill M2; Al, Cr and Ti fill
    M1, with half the Al on M1 (Al coupled between M1 and tetrahedral sites);
    Mg/(Mg+Fe) is the same on both sites. sites="PTEXL" puts
    (Al - Cr - 2Ti + Na)/2 on M1 instead (Nickel & Green, 1985), as PTEXL does.

    Source: Wells (1977), Contrib. Mineral. Petrol. 62, 129-139: eq. 5; site model of Wood & Banno (1973),
    Contrib. Mineral. Petrol. 42, 109-124, eqs. 24-25.
    """
    opx, cpx = c["opx"], c["cpx"]

    def a_en(px):
        m1, m2 = _px_site_fractions(px)
        if sites == "WB73":
            m1 = 1 - px.Al / 2 - px.Cr - px.Ti
        return m1 * m2 * (1 - _x_fe(px)) ** 2

    ln_k = np.log(a_en(cpx) / a_en(opx))
    return 7341 / (3.355 + 2.44 * _x_fe(opx) - ln_k)


def T_ellis_green79(c: Cations, P):
    """Ellis & Green (1979) garnet-cpx Fe-Mg thermometer (their eq. 9).

    Source: Ellis & Green (1979), Contrib. Mineral. Petrol. 71, 13-22: eq. 9.
    """
    gt = c["gt"]
    return (3104 * _x_gt(gt, "Ca") + 3030 + 10.86 * 10 * P) / (_ln_kd_fe_mg(gt, c["cpx"]) + 1.9034)


def T_powell85(c: Cations, P):
    """Powell (1985) garnet-cpx Fe-Mg thermometer (his eq. 17).

    Source: Powell (1985), J. Metamorph. Geol. 3, 231-243: eq. 17.
    """
    gt = c["gt"]
    return (2790 + 10 * 10 * P + 3140 * _x_gt(gt, "Ca")) / (1.735 + _ln_kd_fe_mg(gt, c["cpx"]))


def T_ai94(c: Cations, P):
    """Ai (1994) garnet-cpx Fe-Mg thermometer (his eq. 3); Mg no. is 100 Mg/(Mg+Fe) of garnet.

    Source: Ai (1994), Contrib. Mineral. Petrol. 115, 467-473: eq. 3.
    """
    gt = c["gt"]
    x_ca = _x_gt(gt, "Ca")
    mg_no = 100 * (1 - _x_fe(gt))
    numerator = -1629 * x_ca**2 + 3648.55 * x_ca - 6.59 * mg_no + 1987.98 + 17.66 * 10 * P
    return numerator / (_ln_kd_fe_mg(gt, c["cpx"]) + 1.076)


def T_krogh_ravna00(c: Cations, P):
    """Krogh Ravna (2000) garnet-cpx Fe-Mg thermometer, with Mn in garnet.

    Source: Krogh Ravna (2000), J. Metamorph. Geol. 18, 211-219: unnumbered equation on p. 213.
    """
    gt = c["gt"]
    x_ca, x_mn, x_mg = _x_gt(gt, "Ca"), _x_gt(gt, "Mn"), 1 - _x_fe(gt)
    numerator = (
        1939.9
        + 3270 * x_ca
        - 1396 * x_ca**2
        + 3319 * x_mn
        - 3535 * x_mn**2
        + 1105 * x_mg
        - 3561 * x_mg**2
        + 2324 * x_mg**3
        + 169.4 * P
    )
    return numerator / (_ln_kd_fe_mg(gt, c["cpx"]) + 1.223)


def T_harley84(c: Cations, P):
    """Harley (1984) garnet-opx Fe-Mg thermometer (his eq. 11).

    Source: Harley (1984), Contrib. Mineral. Petrol. 86, 359-373: eq. 11.
    """
    gt = c["gt"]
    return (3740 + 1400 * _x_gt(gt, "Ca") + 22.86 * 10 * P) / (1.987 * _ln_kd_fe_mg(gt, c["opx"]) + 1.96)


# Brey & Kohler (1990, Table 1) Fe-Mg exchange: -R T ln K_D = dH - T dS + P dV,
# K_D = (Fe/Mg)_A / (Fe/Mg)_B; values are (dH/R, dS/R, dV/R) with P in kbar
_BK90_FE_MG = {("ol", "cpx"): (-1489, -1.07, -4.64), ("gt", "opx"): (-1456, -0.55, -9.86)}


def _T_brey_kohler90_fe_mg(c: Cations, P, a: str, b: str):
    dH, dS, dV = _BK90_FE_MG[(a, b)]
    return (dH + dV * 10 * P) / (dS - _ln_kd_fe_mg(c[a], c[b]))


def T_brey_kohler90_ol_cpx(c: Cations, P):
    """Brey & Kohler (1990, Table 1) olivine-cpx Fe-Mg exchange; not recommended by them.

    Source: Brey & Kohler (1990), J. Petrol. 31, 1353-1378: Table 1.
    """
    return _T_brey_kohler90_fe_mg(c, P, "ol", "cpx")


def T_brey_kohler90_gt_opx(c: Cations, P):
    """Brey & Kohler (1990, Table 1) garnet-opx Fe-Mg exchange; not recommended by them.

    Source: Brey & Kohler (1990), J. Petrol. 31, 1353-1378: Table 1.
    """
    return _T_brey_kohler90_fe_mg(c, P, "gt", "opx")


def _kohler_brey90_ln_d(c: Cations) -> pd.Series:
    """ln D_Ca = ln(Ca_ol / Ca_cpx), cations per 4 and 6 O."""
    return np.log(c["ol"].Ca / c["cpx"].Ca)


def T_kohler_brey90(c: Cations, P):
    """Kohler & Brey (1990) Ca exchange between olivine and cpx, their eqs. 5 and 6 solved for T.

    The high-T branch is used where it gives T above T_split = 1275.25 + 2.827 P[kbar].

    Source: Kohler & Brey (1990), Geochim. Cosmochim. Acta 54, 2375-2388: eqs. 5 (high T) and 6 (low T).
    """
    ln_d, P_kbar = _kohler_brey90_ln_d(c), 10 * P
    T_high = (11982 + 56.2 * P_kbar) / (3.61 - ln_d)
    T_low = (5792 + 42.5 * P_kbar) / (-1.25 - ln_d)
    return T_high.where(T_high >= 1275.25 + 2.827 * P_kbar, T_low)


def T_oneill_wood79(c: Cations, P):
    """O'Neill & Wood (1979) olivine-garnet Fe-Mg thermometer (their eq. 26).

    The volume term DV is that of the correction by O'Neill (1980), which
    reverses the sign of eq. 27 as printed in 1979 and fixes its brackets. It
    depends on T, so T is solved for numerically.

    Source: O'Neill & Wood (1979), Contrib. Mineral. Petrol. 70, 59-70: eq. 26; DV from O'Neill (1980), Contrib.
    Mineral. Petrol. 72, 337.
    """
    ol, gt = c["ol"], c["gt"]
    P_kbar = 10 * P
    x_ol = 1 - 2 * _x_fe(ol)  # X_Mg - X_Fe in olivine
    composition = 902 + x_ol * (498 + 1.51 * (P_kbar - 30)) - 98 * (_x_gt(gt, "Mg") - _x_gt(gt, "Fe"))
    composition = composition + 1347 * _x_gt(gt, "Ca")
    ln_kd = _ln_kd_fe_mg(gt, ol)

    def f(T):
        dT = T - 1073
        dv_term = (
            -462.5 * (1.0191 + 2.87e-5 * dT) * (P_kbar - 2.63e-4 * P_kbar**2 - 29.76)
            - 262.4 * (1.0292 + 4.5e-5 * dT) * (P_kbar - 3.9e-4 * P_kbar**2 - 29.65)
            + 454 * (1.020 + 2.84e-5 * dT) * (P_kbar - 2.36e-4 * P_kbar**2 - 29.79)
            + 278.3 * (1.0234 + 2.3e-5 * dT) * (P_kbar - 4.5e-4 * P_kbar**2 - 29.6)
        )
        return (composition + dv_term) / (ln_kd + 0.357)

    return _solve_T(f, len(ol))


def _log_ppm(wt: pd.DataFrame, oxide: str) -> pd.Series:
    """ln of the cation concentration [ppm]; NaN where the oxide is zero (below detection)."""
    ppm = element_ppm(wt, oxide)
    return np.log(ppm.where(ppm > 0))


def _ol_cr_number(ol: pd.DataFrame) -> pd.Series:
    return ol.Cr / (ol.Cr + ol.Al)


def T_de_hoog10_al(c: Cations, P, wt, version: str = "eq12"):
    """De Hoog et al. (2010) Al-in-olivine thermometer for garnet peridotites.

    version="eq12" includes olivine Cr/(Cr+Al); "eq13" drops it, for olivines
    whose Cr cannot be measured or may be Cr2+ (e.g. inclusions in diamond).

    Source: De Hoog et al. (2010), Chem. Geol. 270, 196-215: eq. 12 (version="eq12") or eq. 13 ("eq13").
    """
    al = _log_ppm(wt["ol"], "Al2O3")
    P_kbar = 10 * P
    if version == "eq13":
        return (11959 + 55.6 * P_kbar) / (14.530 - al)
    return (9423 + 51.4 * P_kbar + 1860 * _ol_cr_number(c["ol"])) / (13.409 - al)


def T_de_hoog10_cr(c: Cations, P, wt):
    """De Hoog et al. (2010) Cr-in-olivine thermometer (their eq. 14).

    Source: De Hoog et al. (2010), Chem. Geol. 270, 196-215: eq. 14.
    """
    cr = _log_ppm(wt["ol"], "Cr2O3")
    return (13444 + 48.5 * 10 * P - 4678 * _ol_cr_number(c["ol"])) / (14.53 - cr)


def T_bussweiler17(c: Cations, P, wt, version: str = "eq3"):
    """Bussweiler et al. (2017) experimental Al-in-olivine thermometer.

    version="eq3" (default) uses P and Al only, which the authors prefer as
    Cr# correlates with T in their experiments; "eq2" includes olivine Cr#.

    Source: Bussweiler et al. (2017), Lithos 272-273, 301-314: eq. 3 (version="eq3") or eq. 2 ("eq2").
    """
    al = _log_ppm(wt["ol"], "Al2O3")
    P_kbar = 10 * P
    if version == "eq2":
        return (14750 + 45.1 * P_kbar - 2831 * _ol_cr_number(c["ol"])) / (14.97 - al)
    return (11245 + 46.0 * P_kbar) / (13.68 - al)


def T_ryan96_zn(c: Cations, P, wt, coefficients: str = "RY96"):
    """Ryan et al. (1996) Zn-in-chromite thermometer (their eq. 2), independent of P.

    coefficients="PTEXL" uses -0.904 and 0.264 as in PTEXL, instead of the
    published -0.90 and 0.26; the source of these is not given.

    Source: Ryan et al. (1996), J. Geophys. Res. 101(B3), 5611-5625: eq. 2.
    """
    a, b = (-0.904, 0.264) if coefficients == "PTEXL" else (-0.90, 0.26)
    return 1000 / (a + b * _log_ppm(wt["sp"], "ZnO"))


def T_creighton09(c: Cations, P, cr_term: str = "corrected"):
    """Creighton (2009) single-crystal Mn-in-garnet thermometer (his eq. 10), independent of P.

    The X_Cr coefficient of ln a_spessartine is printed as -31.5 (his eq. 6).
    That makes ln a_spessartine about 5 units below ideal mixing for typical
    mantle garnets and T near 1800 C, so it is taken as a misprint for -3.15,
    the value in PTEXL (cr_term="corrected", default). cr_term="printed" uses -31.5.

    Source: Creighton (2009), Lithos 112S, 177-182: eq. 10, with ln a(spessartine) from eq. 6.
    """
    gt = c["gt"]
    x_cr = gt.Cr / (gt.Cr + gt.Al)
    k_cr = 31.5 if cr_term == "printed" else 3.15
    ln_a_spess = 401.6 * _x_gt(gt, "Mn") - k_cr * x_cr - 17.85
    a_pyr = 1.163 * _x_gt(gt, "Mg") - 0.248 * x_cr - 0.46
    return -1635.48 - 184.8 * (ln_a_spess - np.log(a_pyr)) + C2K


def T_sudholz22_gt_cpx(c: Cations, P):
    """Sudholz et al. (2022) garnet-cpx Fe-Mg thermometer (their eq. 1), calibrated at 1-7 GPa.

    X_Ca and X_Mg are Ca and Mg / (Ca + Fe + Mg) in garnet (no Mn), and
    Jd = Na - Cr - 2Ti in cpx.

    Source: Sudholz et al. (2022), Contrib. Mineral. Petrol. 177, 77: eq. 1.
    """
    gt, cpx = c["gt"], c["cpx"]
    x = gt.Ca + gt.Fe + gt.Mg
    jd = cpx.Na - cpx.Cr - 2 * cpx.Ti
    denominator = -0.008 * 10 * P + 0.259 * gt.Ca / x + 0.914 * gt.Mg / x - 0.159 * jd
    return 3356.34 / (denominator + _ln_kd_fe_mg(gt, cpx) + 1.265)


def T_sudholz22_gt_opx(c: Cations, P):
    """Sudholz et al. (2022) garnet-opx Fe-Mg thermometer (their eq. 2); X_Ca = Ca / (Ca + Fe + Mg) in garnet.

    Source: Sudholz et al. (2022), Contrib. Mineral. Petrol. 177, 77: eq. 2.
    """
    gt = c["gt"]
    x_ca = gt.Ca / (gt.Ca + gt.Fe + gt.Mg)
    return 1851.85 / (-0.007 * 10 * P - 1.83 * x_ca + _ln_kd_fe_mg(gt, c["opx"]) + 1.08)


# Holdaway (2001) garnet activity polynomials, as given by Wu & Zhao (2007,
# Table S2, eqs. A4-A9): RT ln(gamma) of almandine (Fe) and pyrope (Mg) is
# a*T + b*P[bar] + c [J]. Each term is (coefficient, cations multiplied).
_CUBIC = [
    (11.006, "Fe Fe Mg"), (-0.674, "Fe Fe Ca"), (-11.986, "Mg Mg Fe"), (37.966, "Mg Mg Ca"), (46.02, "Mg Mg Mn"),
    (19.34, "Ca Ca Fe"), (25.746, "Ca Ca Mg"), (46.02, "Mn Mn Mg"), (19.145, "Fe Mg Ca"), (45.53, "Fe Mg Mn"),
    (9.333, "Fe Ca Mn"), (77.876, "Mg Ca Mn"),
]  # fmt: skip
_CUBIC_B = [
    (-0.1, "Fe Fe Mg"), (-0.08, "Fe Fe Ca"), (-0.048, "Fe Fe Mn"), (0.068, "Mg Mg Fe"), (-0.136, "Mg Mg Ca"),
    (-0.076, "Mg Mg Mn"), (-0.27, "Ca Ca Fe"), (-0.28, "Ca Ca Mg"), (-0.13, "Ca Ca Mn"), (-0.048, "Mn Mn Fe"),
    (-0.076, "Mn Mn Mg"), (-0.13, "Mn Mn Ca"), (0.103, "Fe Mg Ca"), (-0.14, "Fe Mg Mn"), (-0.353, "Fe Ca Mn"),
    (-0.414, "Mg Ca Mn"),
]  # fmt: skip
_CUBIC_C = [
    (-23244.0, "Fe Fe Mg"), (2608.0, "Fe Fe Ca"), (-3234.0, "Fe Fe Mn"), (11344.0, "Mg Mg Fe"), (-132228.0, "Mg Mg Ca"),
    (-82498.0, "Mg Mg Mn"), (-39864.0, "Ca Ca Fe"), (-51518.0, "Ca Ca Mg"), (-2850.0, "Ca Ca Mn"),
    (-3234.0, "Mn Mn Fe"), (-82498.0, "Mn Mn Mg"), (-2850.0, "Mn Mn Ca"), (-62653.0, "Fe Mg Ca"),
    (-91682.0, "Fe Mg Mn"), (-24712.0, "Fe Ca Mn"), (-177221.0, "Mg Ca Mn"),
]  # fmt: skip
HOLDAWAY01 = {
    "Fea": [
        (5.993, "Mg Mg"), (-9.67, "Ca Ca"), (-11.006, "Fe Mg"), (0.674, "Fe Ca"), (-9.5725, "Mg Ca"),
        (-22.765, "Mg Mn"), (-4.6665, "Ca Mn"), *_CUBIC,
    ],
    "Feb": [
        (-0.034, "Mg Mg"), (0.135, "Ca Ca"), (0.024, "Mn Mn"), (0.1, "Fe Mg"), (0.08, "Fe Ca"), (0.048, "Fe Mn"),
        (-0.0515, "Mg Ca"), (0.07, "Mg Mn"), (0.1765, "Ca Mn"), *_CUBIC_B,
    ],
    "Fec": [
        (-5672.0, "Mg Mg"), (19932.0, "Ca Ca"), (1617.0, "Mn Mn"), (23244.0, "Fe Mg"), (-2608.0, "Fe Ca"),
        (3234.0, "Fe Mn"), (31326.5, "Mg Ca"), (45841.0, "Mg Mn"), (12356.0, "Ca Mn"), *_CUBIC_C,
    ],
    "Mga": [
        (-5.503, "Fe Fe"), (-12.873, "Ca Ca"), (-23.01, "Mn Mn"), (11.986, "Fe Mg"), (-9.5725, "Fe Ca"),
        (-22.765, "Fe Mn"), (-37.966, "Mg Ca"), (-46.02, "Mg Mn"), (-38.938, "Ca Mn"), *_CUBIC,
    ],
    "Mgb": [
        (0.05, "Fe Fe"), (0.14, "Ca Ca"), (0.038, "Mn Mn"), (-0.068, "Fe Mg"), (-0.0515, "Fe Ca"), (0.07, "Fe Mn"),
        (0.136, "Mg Ca"), (0.076, "Mg Mn"), (0.207, "Ca Mn"), *_CUBIC_B,
    ],
    "Mgc": [
        (11622.0, "Fe Fe"), (25759.0, "Ca Ca"), (41249.0, "Mn Mn"), (-11344.0, "Fe Mg"), (31326.5, "Fe Ca"),
        (45841.0, "Fe Mn"), (132228.0, "Mg Ca"), (82498.0, "Mg Mn"), (88610.5, "Ca Mn"), *_CUBIC_C,
    ],
}  # fmt: skip


def _holdaway01(name: str, x: Mapping[str, pd.Series]) -> pd.Series:
    """Evaluate one of the HOLDAWAY01 polynomials for garnet X-site fractions x["Fe"], x["Mg"], ..."""
    return sum(coef * np.prod([x[cation] for cation in term.split()], axis=0) for coef, term in HOLDAWAY01[name])


def T_wu_zhao07(c: Cations, P):
    """Wu & Zhao (2007) olivine-garnet Fe-Mg thermometer (their eq. 7).

    Garnet activities follow Holdaway (2001), with the polynomials of their
    Table S2; olivine is a symmetric Fe-Mg solution. X-site fractions of
    garnet are over Fe + Mg + Ca + Mn; all Fe as Fe2+.

    Source: Wu & Zhao (2007), J. Metamorph. Geol. 25, 497-505: eq. 7; garnet activities from Holdaway (2001)
    with the polynomials of their supplementary Table S2 (eqs. A4-A9).
    """
    ol, gt = c["ol"], c["gt"]
    x = {cation: _x_gt(gt, cation) for cation in ("Fe", "Mg", "Ca", "Mn")}
    h = {name: _holdaway01(name, x) for name in HOLDAWAY01}
    P_bar = 1e4 * P
    x_ol = 2 * _x_fe(ol) - 1  # X_Fe - X_Mg in olivine
    numerator = 58136.3 + P_bar * (0.406 - 2 * h["Feb"] + 2 * h["Mgb"]) - 2 * h["Fec"] + 2 * h["Mgc"] - 23653.8 * x_ol
    denominator = 38.743 + 2 * h["Fea"] - 2 * h["Mga"] + 6 * 8.314 * _ln_kd_fe_mg(gt, ol)
    return numerator / denominator


# --------------------------------------------------------------------------
# Spinel peridotites: olivine-spinel and opx-spinel
# --------------------------------------------------------------------------


def spinel_sites(sp: pd.DataFrame) -> pd.DataFrame:
    """Spinel cations per 3 cations and 4 O, with Fe3+ from charge balance (Droop, 1987).

    Columns Mg, Fe2, Fe3, Al, Cr, Ti, from cations per 4 O with all Fe as Fe2+
    (chemistry.cations_per_formula_unit). Fe3+ is limited to 0..Fe.
    """
    total = sp.sum(axis=1, skipna=False)
    scaled = sp.mul(3 / total, axis=0)
    fe3 = (8 * (1 - 3 / total)).clip(lower=0).clip(upper=scaled.Fe)
    return pd.DataFrame(
        {"Mg": scaled.Mg, "Fe2": scaled.Fe - fe3, "Fe3": fe3, "Al": scaled.Al, "Cr": scaled.Cr, "Ti": scaled.Ti}
    )


def _ol_sp_terms(c: Cations):
    """Olivine-spinel Mg-Fe2+ exchange: ln K_D = ln[(Mg/Fe)_ol / (Mg/Fe2+)_sp], and spinel sites."""
    ol, sp = c["ol"], spinel_sites(c["sp"])
    return np.log(ol.Mg * sp.Fe2 / (ol.Fe * sp.Mg)), ol, sp


def T_oneill_wall87(c: Cations, P):
    """O'Neill & Wall (1987) olivine-spinel Mg-Fe2+ thermometer (their eq. 20).

    Spinel mole fractions are per two octahedral sites (X_Cr = Cr/2, X_Fe3+ =
    Fe3+/2), except X_Ti = Ti per 4 O; X_Mg and X_Fe2+ are Mg and Fe2+ /
    (Mg + Fe2+). Fe3+ is from spinel stoichiometry (spinel_sites).

    Source: O'Neill & Wall (1987), J. Petrol. 28, 1169-1191: eq. 20.
    """
    ln_kd, ol, sp = _ol_sp_terms(c)
    P_kbar = 10 * P
    x_ol = 1 - 2 * _x_fe(ol)  # X_Mg - X_Fe in olivine
    x_sp = (sp.Mg - sp.Fe2) / (sp.Mg + sp.Fe2)
    numerator = (
        6530
        + 28 * P_kbar
        + (5000 + 10.8 * P_kbar) * x_ol
        - 1960 * (1 + sp.Ti) * x_sp
        + 18620 * sp.Cr / 2
        + 25150 * (sp.Fe3 / 2 + sp.Ti)
    )
    return numerator / (8.31441 * ln_kd + 4.705)


def T_ballhaus91(c: Cations, P):
    """Ballhaus et al. (1991) olivine-spinel Mg-Fe2+ thermometer, a refit of O'Neill & Wall (1987).

    As printed, a bracket is misplaced: 6530 + 280P + (7000 + 108P)(1 - 2X_Fe^ol),
    with P in GPa, is read by analogy with O'Neill & Wall's eq. 20. Spinel
    mole fractions are as in T_oneill_wall87.

    Source: Ballhaus et al. (1991), Contrib. Mineral. Petrol. 107, 27-40: unnumbered equation on p. 35.
    """
    ln_kd, ol, sp = _ol_sp_terms(c)
    x_sp = (sp.Mg - sp.Fe2) / (sp.Mg + sp.Fe2)
    numerator = (
        6530
        + 280 * P
        + (7000 + 108 * P) * (1 - 2 * _x_fe(ol))
        - 1960 * x_sp
        + 16150 * sp.Cr / 2
        + 25150 * (sp.Fe3 / 2 + sp.Ti)
    )
    return numerator / (8.31441 * ln_kd + 4.705)


def T_li95(c: Cations, P):
    """Li et al. (1995) olivine-spinel Mg-Fe2+ thermometer (their eqs. 9-10), independent of P.

    Y_i = i / (Cr + Al + Fe3+) in spinel; K_D is corrected to a Fe3+-free basis.

    Source: Li, Kornprobst, Vielzeuf & Fabries (1995), Chinese J. Geochem. 14, 68-77: eqs. 9-10.
    """
    ln_kd, ol, sp = _ol_sp_terms(c)
    y = sp.Cr + sp.Al + sp.Fe3
    y_cr = sp.Cr / y
    ln_kd0 = ln_kd - 2.0 * sp.Fe3 / y
    return (4299 * y_cr + 1283) / (ln_kd0 + 1.469 * y_cr + 0.363)


def _ln_kd_al_ol_sp(wt) -> pd.Series:
    """ln(Al2O3 in olivine / Al2O3 in spinel), wt%; NaN where olivine Al2O3 is zero."""
    al_ol = wt["ol"].Al2O3
    return np.log(al_ol.where(al_ol > 0) / wt["sp"].Al2O3)


def T_wan08(c: Cations, P, wt):
    """Wan et al. (2008) Al exchange between olivine and spinel, at 100 kPa; little P effect.

    Source: Wan et al. (2008), Am. Mineral. 93, 1142-1147: eq. 2.
    """
    sp = c["sp"]
    return 10000 / (0.512 + 0.873 * sp.Cr / (sp.Cr + sp.Al) - 0.91 * _ln_kd_al_ol_sp(wt))


def T_coogan14(c: Cations, P, wt):
    """Coogan et al. (2014) Al exchange between olivine and spinel (their eqs. 6-7), a refit of Wan et al. (2008).

    D'Souza et al. (2020) found no P effect at spinel-peridotite conditions,
    and pair it with the Kohler & Brey (1990) Ca-in-olivine barometer.

    Source: Coogan et al. (2014), Chem. Geol. 368, 1-10: eqs. 6-7.
    """
    sp = c["sp"]
    return 10000 / (0.575 + 0.884 * sp.Cr / (sp.Cr + sp.Al) - 0.897 * _ln_kd_al_ol_sp(wt))


def T_witt_eickschen_seck91(c: Cations, P, version: str = "opx-sp"):
    """Witt-Eickschen & Seck (1991) Al solubility in opx from spinel peridotites, independent of P.

    version="opx-sp" (their eq. 6) uses
    K_D = X_Fo^2 X_Al^M1(opx) / [X_Al^2(sp) X_Mg(sp) X_Mg^M1(opx)] and spinel
    Y_Cr = Cr/(Cr + Al + Fe3+). X_Mg^M1 of opx is taken as
    (1 - X_Al^M1 - Cr - Ti) Mg/(Mg + Fe), which the paper does not specify.
    version="Cr-Al-opx" (their eq. 7) uses opx alone: Al and Cr on M1. Al on
    M1 is 0.5 (Al - Cr - 2Ti + Na), after Carswell (1989).

    Source: Witt-Eickschen & Seck (1991), Contrib. Mineral. Petrol. 106, 431-439: eq. 6 (version="opx-sp") or
    eq. 7 ("Cr-Al-opx").
    """
    opx = c["opx"]
    al_m1 = _opx_al_m1(opx)
    if version == "Cr-Al-opx":
        return 636.54 + 2088.21 * al_m1 + 14527.32 * opx.Cr + C2K
    ol, sp = c["ol"], spinel_sites(c["sp"])
    y = sp.Cr + sp.Al + sp.Fe3
    y_cr = sp.Cr / y
    mg_m1 = (1 - al_m1 - opx.Cr - opx.Ti) * (1 - _x_fe(opx))
    kd = (1 - _x_fe(ol)) ** 2 * al_m1 / ((sp.Al / y) ** 2 * sp.Mg / (sp.Mg + sp.Fe2) * mg_m1)
    ln_kd = np.log(kd)
    return 2248.25 + 991.58 * ln_kd + 153.32 * ln_kd**2 + 539.05 * y_cr - 2005.74 * y_cr**2 + C2K


# --------------------------------------------------------------------------
# Ni in garnet: Ni partitioning between garnet and olivine. Most calibrations
# assume the Ni content of mantle olivine (about 2900-3000 ppm), so they need
# garnet only; ni_ol="measured" uses analysed olivine NiO instead.
# --------------------------------------------------------------------------


def _ln_d_ni(wt, ni_ol: str, assumed: float) -> pd.Series:
    """ln(Ni_gt / Ni_ol), ppm, with Ni_ol measured or assumed."""
    ln_ol = _log_ppm(wt["ol"], "NiO") if ni_ol == "measured" else np.log(assumed)
    return _log_ppm(wt["gt"], "NiO") - ln_ol


def T_griffin89_ni(c: Cations, P, wt, ni_ol: str = "assumed"):
    """Griffin et al. (1989) Ni-in-garnet thermometer, independent of P.

    ni_ol="assumed" is their eq. IIa, for olivine with 3000 ppm Ni;
    ni_ol="measured" is their eq. I, with K_D = Ni_gt / Ni_ol.

    Source: Griffin et al. (1989), Contrib. Mineral. Petrol. 103, 199-202: eq. IIa (ni_ol="assumed") or eq. I
    ("measured").
    """
    if ni_ol == "measured":
        log_100kd = 2 + _ln_d_ni(wt, "measured", np.nan) / np.log(10)
        return 1000 / (-0.428 * log_100kd + 0.84)
    log_ni = (_log_ppm(wt["gt"], "NiO") - np.log(30)) / np.log(10)
    return 1000 / (-0.435 * log_ni + 0.83)


def T_ryan96_ni(c: Cations, P, wt):
    """Ryan et al. (1996) Ni-in-garnet thermometer (their eq. 1), for olivine with 2900 ppm Ni; independent of P.

    Source: Ryan et al. (1996), J. Geophys. Res. 101(B3), 5611-5625: eq. 1.
    """
    return 1000 / (1.506 - 0.189 * _log_ppm(wt["gt"], "NiO"))


def T_canil99(c: Cations, P, wt, ni_ol: str = "measured"):
    """Canil (1999) experimental Ni-in-garnet thermometer (his eq. 3), independent of P.

    ni_ol="assumed" takes olivine Ni as 3000 ppm, as Sudholz et al. (2021) do.

    Source: Canil (1999), Contrib. Mineral. Petrol. 136, 240-246: eq. 3.
    """
    return 8772 / (2.53 - _ln_d_ni(wt, ni_ol, 3000))


def T_sudholz21_ni(c: Cations, P, wt, ni_ol: str = "assumed"):
    """Sudholz et al. (2021) experimental Ni-in-garnet thermometer (their eq. 4), independent of P.

    X_Ca = Ca / (Ca + Fe + Mg) and X_Cr = Cr / (Cr + Al) in garnet. They
    recommend assuming 3000 ppm Ni in olivine for single garnets (default).

    Source: Sudholz et al. (2021), Contrib. Mineral. Petrol. 176, 32: eq. 4.
    """
    gt = c["gt"]
    x_ca = gt.Ca / (gt.Ca + gt.Fe + gt.Mg)
    x_cr = gt.Cr / (gt.Cr + gt.Al)
    return -8254.568 / (3.023 * x_ca + 2.307 * x_cr + _ln_d_ni(wt, ni_ol, 3000) - 2.639)


# Nimis et al. (2024) eqs. 3-7: P-independent versions along conductive
# geotherms of Hasterok & Chapman (2011), surface heat flow in mW/m2
_NIMIS24_GEOTHERMS = {
    "T35": (6452, 0.451),
    "T40": (6735, 0.955),
    "T45": (6889, 1.25),
    "T50": (7197, 1.52),
    "T60": (7099, 1.67),
}


def T_nimis24_ni(c: Cations, P, wt, version: str = "eq2", ni_ol: str = "assumed"):
    """Nimis et al. (2024) Ni-in-garnet thermometer, with a P term (their eq. 2).

    Olivine Ni is assumed to be 2900 ppm, as in their calibration, unless
    ni_ol="measured". version="T35" ... "T60" are their eqs. 3-7, single-garnet
    versions along 35-60 mW/m2 conductive geotherms, independent of P (and
    always with 2900 ppm in olivine).

    Source: Nimis et al. (2024), Lithos 468-469, 107513: eq. 2 (version="eq2") or eqs. 3-7 ("T35" to "T60").
    """
    if version in _NIMIS24_GEOTHERMS:
        a, b = _NIMIS24_GEOTHERMS[version]
        return a / (b - _ln_d_ni(wt, "assumed", 2900))
    return (7451 + 243.6 * P) / (2.31 - _ln_d_ni(wt, ni_ol, 2900))


# --------------------------------------------------------------------------
# Barometers: f(c, T [K]) -> P [GPa]
# --------------------------------------------------------------------------


def P_nickel_green85(c: Cations, T, mgts: str = "NG85", mgts_in: str = "activity_dV"):
    """Nickel & Green (1985) Al-in-opx barometer for opx-garnet pairs (their eq. 5).

    NG85 use X_Al^M1 = (Al - Cr - 2Ti + Na) / 2 in the equilibrium constant, in
    dV and in the M1 Mg+Fe occupancy. `mgts` replaces it with a later Mg-Tschermak
    scheme (a key of OPX_MGTS), and `mgts_in` sets where the replacement is used:
    "activity" (equilibrium constant only), "activity_dV" (equilibrium constant
    and dV, as in PTEXL), or "all" (also the M1 Mg+Fe occupancy).

    Source: Nickel & Green (1985), Earth Planet. Sci. Lett. 73, 158-170: eq. 5. mgts="CA91": Carswell (1991),
    Mineral. Mag. 55, 19-31, eqs. E6 and E9; "CG87": eq. E7 there; "T98": Taylor (1998), N. Jb. Miner. Abh.
    172, 381-408, p. 402.
    """
    opx, gt = c["opx"], c["gt"]
    al_m1 = _opx_al_m1(opx)
    x_mgts = OPX_MGTS[mgts](opx)
    x_dv = x_mgts if mgts_in in ("activity_dV", "all") else al_m1
    x_occupancy = x_mgts if mgts_in == "all" else al_m1

    # Opx M1 and M2 occupancies
    mf_m1 = 1 - x_occupancy - opx.Cr - opx.Ti
    mf_m2 = 1 - opx.Ca - opx.Na - opx.Mn
    mg_no = opx.Mg / (opx.Mg + opx.Fe)
    mg_m1 = mf_m1 * mg_no
    fe_m1 = mf_m1 * (1 - mg_no)

    # Garnet X and Y site fractions
    x_site = gt.Ca + gt.Mg + gt.Fe + gt.Mn
    x_ca = gt.Ca / x_site
    x_fe = gt.Fe / x_site
    x_al = gt.Al / (gt.Al + gt.Cr)
    x_cr = gt.Cr / (gt.Al + gt.Cr)

    ln_k = np.log((1 - x_ca) ** 3 * x_al**2 / (mf_m1 * mf_m2**2 * x_mgts))
    dV = -1 / (183.3 + 178.98 * x_dv * (1 - x_dv))
    # X_Cr^M1 is not defined by NG85; it is taken as Cr per formula unit, as in PTEXL
    nonideal = (
        -9000 * x_ca**2
        - 6800 * x_cr**2
        + 3400 * mg_m1 * opx.Cr
        - 90853 * x_ca * x_cr
        - 7590 * x_fe * x_ca
        + 5157 * mg_m1 * fe_m1
        + 6047
    )
    P_kbar = dV * (nonideal + T * (-1.9872 * ln_k + 52.1 * x_ca * x_cr - 3.23))
    return 0.1 * P_kbar


def P_sudholz21(c: Cations, T):
    """Sudholz et al. (2021) Cr-in-cpx barometer.

    Source: Sudholz et al. (2021), Contrib. Mineral. Petrol. 176, 11: eq. 5.
    """
    cpx = c["cpx"]
    return 11.03 - 0.001088 * T * np.log(cacrts_activity(cpx)) + 1.526 * np.log(cpx_cr_number(cpx) / T)


def _al_in_opx_terms(c: Cations):
    """Opx and garnet terms shared by the Brey & Kohler (1990) and Brey et al. (2008) barometers.

    Returns ln K_D and the products of site fractions used in their non-ideal terms.
    """
    opx, gt = c["opx"], c["gt"]
    mf_m1, mf_m2 = _px_site_fractions(opx)
    x_mg = 1 - _x_fe(opx)
    mg_fe_m1 = mf_m1 * x_mg * mf_m1 * (1 - x_mg)
    x_ca, x_fe = _x_gt(gt, "Ca"), _x_gt(gt, "Fe")
    x_al = gt.Al / (gt.Al + gt.Cr)
    x_cr = 1 - x_al
    # Al on M1 in Tschermak's component, (Al - |Na - Cr - 2Ti|) / 2 (Carswell & Gibb, 1987)
    x_al_ts = _mgts_carswell_gibb87(opx)
    ln_kd = np.log((1 - x_ca) ** 3 * x_al**2 / (mf_m1 * mf_m2**2 * x_al_ts))
    return ln_kd, mg_fe_m1, x_ca, x_fe, x_cr


def P_brey_kohler90(c: Cations, T):
    """Brey & Kohler (1990) Al-in-opx barometer for opx-garnet pairs, "P_BKN" (their Table 3).

    Source: Brey & Kohler (1990), J. Petrol. 31, 1353-1378: Table 3.
    """
    ln_kd, mg_fe, x_ca, x_fe, x_cr = _al_in_opx_terms(c)
    c1 = (
        -8.3143 * T * ln_kd
        - 5510
        + 88.91 * T
        - 19 * T**1.2
        + 3 * x_ca**2 * 82458
        + mg_fe * (80942 - 46.7 * T)
        - 3 * x_fe * x_ca * 17793
        - x_ca * x_cr * (1.164e6 - 420.4 * T)
        - x_fe * x_cr * (-1.25e6 + 565 * T)
    )
    c2 = -0.832 - 8.78e-5 * (T - 298) + 3 * x_ca**2 * 3.305 - x_ca * x_cr * 13.45 + x_fe * x_cr * 10.5
    c3 = 16.6e-4
    P_kbar = (-c2 - np.sqrt(c2**2 + 4 * c3 * c1 / 1000)) / (2 * c3)
    return 0.1 * P_kbar


def P_brey08(c: Cations, T):
    """Brey et al. (2008) Al-in-opx barometer, a refit of Brey & Kohler (1990) to 6-10 GPa (their eqs. 1-4).

    Source: Brey et al. (2008), J. Petrol. 49, 3-24: eqs. 1-4.
    """
    ln_kd, mg_fe, x_ca, x_fe, x_cr = _al_in_opx_terms(c)
    c1 = (
        -8.3143 * T * ln_kd
        + 4970
        + 84.15 * T
        - 19 * T**1.2
        + 3 * x_ca**2 * 82456
        + mg_fe * (80941 - 36.3 * T)
        - 3 * x_fe * x_ca * 17795
        - x_ca * x_cr * (1.164e6 - 335 * T)
        - x_fe * x_cr * (-1.25e6 + 730 * T)
    )
    c2 = -0.533 - 1.62e-4 * (T - 298) + 3 * x_ca**2 * 3.305 - 18.25 * x_ca * x_cr + 3.5 * x_fe * x_cr
    c3 = -7.2e-4
    return (-c2 - np.sqrt(c2**2 + 0.004 * c3 * c1)) / (20 * c3)


def P_macgregor74(c: Cations, T, wt):
    """Al2O3 in opx coexisting with garnet in MgO-Al2O3-SiO2 (MacGregor, 1974).

    MacGregor gives the Al2O3 isopleths only graphically (his Fig. 2). This is
    the fit to his data by Finnerty & Boyd (1984, Appendix eq. A3), with Al2O3
    as a weight fraction: ln X = (-3736 - 97.1 P[kbar]) / T + 1.46. It ignores
    Ca, Cr and Fe.

    Source: MacGregor (1974), Am. Mineral. 59, 110-119 (data, Fig. 2); fit by Finnerty & Boyd (1984), Geochim.
    Cosmochim. Acta 48, 15-27: Appendix eq. A3.
    """
    x_al2o3 = wt["opx"].Al2O3 / 100
    return 0.1 * ((np.log(x_al2o3) - 1.46) * T + 3736) / -97.1


def P_nimis_taylor00(c: Cations, T, correction: str = "none"):
    """Nimis & Taylor (2000) Cr-in-cpx barometer.

    correction="N20" applies the empirical correction of Nimis et al. (2020,
    Appendix A.1), fitted to agree with opx-garnet barometry in xenoliths.

    Source: Nimis & Taylor (2000), Contrib. Mineral. Petrol. 139, 541-554: eq. 9. correction="N20": Nimis et al.
    (2020), Lithos 376-377, 105729: Appendix A.1.
    """
    cpx = c["cpx"]
    P_kbar = -T / 126.9 * np.log(cacrts_activity(cpx)) + 15.483 * np.log(cpx_cr_number(cpx) / T) + T / 71.38 + 107.8
    P = 0.1 * P_kbar
    if correction == "N20":
        P = 0.05024 * P**2 + 0.7633 * P + 0.1257
    return P


def P_ashchepkov17(c: Cations, T, wt, kd: str = "cations"):
    """Ashchepkov et al. (2011) Jd-Di cpx barometer for peridotites, corrected in Ashchepkov et al. (2017, Appendix).

    K_D = Na / AlCr x Mg / Ca, where AlCr approximates the activity of Al, Cr
    and Ti on M1. The paper prints MgO/CaO but every other term in formula
    units; with Mg/Ca in formula units the pressures agree with Nimis & Taylor
    (2000) to within ~0.3 GPa from 2 to 8 GPa, as the paper reports, whereas the
    wt% ratio (kd="wt") gives 0.5-1.3 GPa less (on the 2026 AusXenos database).

    Source: Ashchepkov et al. (2017), Geosci. Front. 8, 775-795: unnumbered equation in the Appendix (p. 792),
    correcting Ashchepkov et al. (2011).
    """
    cpx = c["cpx"]
    al_cr = (
        (cpx.Al - 0.01) * ((T - 600) / 650) ** 0.75
        + cpx.Cr * (T - 100) / 900
        + (4 * cpx.Ti - 0.0125) / (T - 801) * 650
        + (cpx.Fe - 0.23) * (T - 700) / 21000
    )
    mg_ca = wt["cpx"].MgO / wt["cpx"].CaO if kd == "wt" else cpx.Mg / cpx.Ca
    k_d = cpx.Na / al_cr * mg_ca
    P_kbar = (
        0.028 * k_d * T / (1 + 2.05 * cpx.Fe - 1.15 * cpx.Fe * (T - 600) / 1200)
        - np.log(1273 / T) * 40 * (-cpx.Al + 7 * cpx.Na - 15 * cpx.Ti + 10 * cpx.Cr + cpx.Mg / 4)
        + 7.5 * cpx.Si
        + 50 * (cpx.Na + 0.10 * cpx.Al - 2 * cpx.Ti + 0.05 * cpx.Mg - 0.2 * cpx.Ca - 0.75 * cpx.Na / cpx.Ca)
    )
    return 0.1 * P_kbar


def P_kohler_brey90(c: Cations, T):
    """Kohler & Brey (1990) Ca exchange between olivine and cpx as a barometer (their eqs. 5 and 6).

    The high-T branch is used where T >= T_split = 1275.25 + 2.827 P[kbar] at its pressure.

    Source: Kohler & Brey (1990), Geochim. Cosmochim. Acta 54, 2375-2388: eqs. 5 (high T) and 6 (low T).
    """
    ln_d = _kohler_brey90_ln_d(c)
    P_high = (-T * ln_d - 11982 + 3.61 * T) / 56.2
    P_low = (-T * ln_d - 5792 - 1.25 * T) / 42.5
    return 0.1 * P_high.where(T >= 1275.25 + 2.827 * P_high, P_low)


def _asymmetric_ternary(x1, x2, x3, w12, w21, w13, w31, w23, w32, w123):
    """RT ln(gamma) of component 1 in an asymmetric ternary Margules solution (one site).

    The form of Mukhopadhyay et al. (1993), as used by Ganguly et al. (1996)
    and Beyer et al. (2015, eqs. 14-15).
    """
    return (
        0.5 * w12 * x2 * (1 - x1 + x2 + 2 * x1 * (x1 - x2 - 1))
        + 0.5 * w21 * x2 * (1 - x1 - x2 - 2 * x1 * (x1 - x2 - 1))
        + 0.5 * w13 * x3 * (1 - x1 + x3 + 2 * x1 * (x1 - x3 - 1))
        + 0.5 * w31 * x3 * (1 - x1 - x3 - 2 * x1 * (x1 - x3 - 1))
        + w23 * x2 * x3 * (x2 - x3 - 0.5)
        + w32 * x2 * x3 * (x3 - x2 - 0.5)
        + w123 * x2 * x3 * (1 - 2 * x1)
    )


# Beyer et al. (2015, Table 4) clinopyroxene interaction parameters [J/mol]
_BEYER15_CPX = {
    "M1": {"MgFe": 0.0, "MgAl": 4900.0, "FeAl": 3600.0},
    "M2": {"CaMg": 24900.0, "CaFe": -49300.0, "MgFe": 0.0, "CaNa": 11600.0, "MgNa": -35900.0, "FeNa": 59300.0},
    "T": {"AlSi": 7000.0},
}


def P_beyer15(c: Cations, T):
    """Beyer et al. (2015) garnet-cpx barometer, Al in cpx with garnet (their eqs. 2-19).

    Calibrated on bimineralic eclogite (3-7 GPa, 900-1550 C); suits garnet
    pyroxenites without opx. Garnet is an asymmetric Ca-Mg-Fe solution
    (Ganguly et al., 1996) with X_i = i/(Ca + Mg + Fe); cpx sites are
    Al(IV) = 2 - Si, the rest of Al, Ti and Cr on M1 with Mg + Fe filling it,
    and Ca, Na and the remaining Mg + Fe on M2, with Fe/Mg the same on both
    sites. All Fe is Fe2+, and the small reciprocal terms (their eqs. 20-21)
    are neglected, as in their calibration. Cr is not in their experiments;
    it is put on M1 here.

    Source: Beyer et al. (2015), Contrib. Mineral. Petrol. 169, 18: eqs. 2-19, with the Margules parameters of
    Table 4.
    """
    gt, cpx = c["gt"], c["cpx"]
    R = 8.314

    # Garnet (their eqs. 12-15; Table 4 values are per site, so x3 for 3 sites)
    total = gt.Ca + gt.Mg + gt.Fe
    ca, mg, fe = gt.Ca / total, gt.Mg / total, gt.Fe / total
    w_camg, w_mgca = 21600 - 5.78 * T, 9800 - 5.78 * T
    w_cafe, w_feca = 870 - 1.69 * T, 680 - 1.69 * T
    w_mgfe, w_femg = 2100.0, 700.0
    w_camgfe = ((w_mgca - w_camg) + (w_cafe - w_feca) + (w_femg - w_mgfe)) / 2
    rt_ln_gr = 3 * _asymmetric_ternary(ca, mg, fe, w_camg, w_mgca, w_cafe, w_feca, w_mgfe, w_femg, w_camgfe)
    rt_ln_py = 3 * _asymmetric_ternary(mg, ca, fe, w_mgca, w_camg, w_mgfe, w_femg, w_cafe, w_feca, w_camgfe)

    # Clinopyroxene sites (their eqs. 3-6, 16-17)
    al_iv = 2 - cpx.Si
    al_m1 = cpx.Al - al_iv
    mg_fe = cpx.Mg + cpx.Fe
    mf_m1 = np.minimum(1 - al_m1 - cpx.Ti - cpx.Cr, mg_fe)
    x_fe = _x_fe(cpx)
    mg1, fe1 = mf_m1 * (1 - x_fe), mf_m1 * x_fe
    mg2, fe2 = (mg_fe - mf_m1) * (1 - x_fe), (mg_fe - mf_m1) * x_fe
    ca2, na2 = cpx.Ca, cpx.Na
    al_t, si_t = al_iv / 2, cpx.Si / 2
    x_di = ca2 * mg1 * si_t**2
    x_cats = 4 * ca2 * al_m1 * al_t * si_t

    # Clinopyroxene activity coefficients (their eqs. 18-19)
    w1, w2, wt = _BEYER15_CPX["M1"], _BEYER15_CPX["M2"], _BEYER15_CPX["T"]["AlSi"]
    m2 = (
        w2["CaMg"] * mg2**2
        + w2["CaNa"] * na2**2
        + w2["CaFe"] * fe2**2
        + mg2 * fe2 * (w2["CaMg"] + w2["CaFe"] - w2["MgFe"])
        + mg2 * na2 * (w2["CaMg"] + w2["CaNa"] - w2["MgNa"])
        + na2 * fe2 * (w2["CaNa"] + w2["CaFe"] - w2["FeNa"])
    )
    rt_ln_di = (
        w1["MgAl"] * al_m1**2
        + w1["MgFe"] * fe1**2
        + fe1 * al_m1 * (w1["MgFe"] + w1["MgAl"] - w1["FeAl"])
        + m2
        + wt * al_t**2
    )
    rt_ln_cats = (
        w1["MgAl"] * mg1**2
        + w1["FeAl"] * fe1**2
        + mg1 * fe1 * (w1["MgAl"] + w1["FeAl"] - w1["MgFe"])
        + m2
        + wt * si_t**2
    )

    # K of Di + CaTs = 2/3 Gr + 1/3 Py + ..., with X_Gr = X_Ca^3, X_Py = X_Mg^3 (their eq. 10)
    ln_k = np.log(x_di * x_cats / (ca**2 * mg)) + (rt_ln_di + rt_ln_cats - 2 / 3 * rt_ln_gr - 1 / 3 * rt_ln_py) / (
        R * T
    )
    P_bar = -(16500 - 23.50 * T + R * T * ln_k) / 0.719
    return 1e-4 * P_bar


def P_grutter06(c: Cations, T, wt, geotherm: str = "38"):
    """Grutter et al. (2006) Cr-in-garnet barometer "P38", independent of T.

    Pressure for a garnet saturated in chromite, from garnet Cr2O3 and CaO
    [wt%]; without chromite it is a minimum. Not valid for wehrlitic garnets.
    geotherm (mW/m2) other than "38" applies their correction for hotter
    geotherms, P_CG = 0.997 P38 - 0.388 CG + 14.9 [kbar].

    Source: Grutter et al. (2006), J. Petrol. 47, 801-820: unnumbered equations for P38 on p. 812 and the
    geotherm correction on p. 813.
    """
    cr2o3, cao = wt["gt"].Cr2O3, wt["gt"].CaO
    high_cr = 26.9 + 3.22 * cr2o3 - 3.03 * cao
    low_cr = 9.2 + 36 * (cr2o3 + 1.6) / (cao + 7.02)
    P_kbar = high_cr.where(cr2o3 >= 0.94 * cao + 5, low_cr)
    if geotherm != "38":
        P_kbar = 0.997 * P_kbar - 0.388 * float(geotherm) + 14.9
    return 0.1 * P_kbar


# --------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Calibration:
    func: Callable
    reference: str
    # Mineral -> oxides that must have been analysed (not blank) for a result
    oxides: Mapping[str, tuple[str, ...]]
    # Mineral -> oxides that are taken as zero if blank, with a warning
    minor_oxides: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    # Option name -> allowed values; the first value is the published version
    options: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    # Options that differ from the published defaults for this entry
    preset: Mapping[str, str] = field(default_factory=dict)
    # (option, value) -> {mineral: oxides} also required when that option is chosen
    option_oxides: Mapping[tuple[str, str], Mapping[str, tuple[str, ...]]] = field(default_factory=dict)
    # Whether `func` also takes `wt`, the oxide wt% of each mineral (for trace-element calibrations)
    needs_wt: bool = False

    @property
    def minerals(self) -> tuple[str, ...]:
        return tuple(self.oxides)

    def required_oxides(self, options: Mapping[str, str] | None = None) -> dict[str, tuple[str, ...]]:
        """Required oxides with the given options (resolved with `resolve`)."""
        options = self.resolve(options)
        required = dict(self.oxides)
        for (name, value), extra in self.option_oxides.items():
            if options.get(name) == value:
                for mineral, oxides in extra.items():
                    required[mineral] = required.get(mineral, ()) + oxides
        return required

    def resolve(self, given: Mapping[str, str] | None = None) -> dict[str, str]:
        """Options to call `func` with: defaults, then this entry's preset, then `given`."""
        options = {name: values[0] for name, values in self.options.items()}
        options.update(self.preset)
        for name, value in (given or {}).items():
            if name not in self.options:
                raise ValueError(f"unknown option '{name}' (options: {', '.join(self.options) or 'none'})")
            if value not in self.options[name]:
                raise ValueError(f"option {name} must be one of {', '.join(self.options[name])}, not '{value}'")
            options[name] = value
        return options


# Required oxides are the main constituents, needed for a meaningful formula
# unit, plus any minor oxide whose blank would typically move the result by more
# than about half the calibration's uncertainty (judged by zeroing each oxide in
# turn on the 2026 AusXenos database). Other oxides the calibration uses are
# minor: often blank when below detection, and taken as zero with a warning.
# E.g. a blank cpx Cr2O3 gives no Cr-in-cpx pressure, but a blank opx Cr2O3
# moves NG85 pressures by only ~0.02 GPa (median; 0.2 GPa at the 90th percentile).
_BASE = ("SiO2", "Al2O3", "FeOt", "MgO", "CaO")
_OL = ("SiO2", "FeOt", "MgO")
_SP = ("Al2O3", "Cr2O3", "FeOt", "MgO")

_NG85_OPTIONS = {"mgts": tuple(OPX_MGTS), "mgts_in": ("activity_dV", "activity", "all")}
_AL_IN_OPX_OXIDES = {"opx": _BASE, "gt": _BASE}
_AL_IN_OPX_MINOR = {"opx": ("TiO2", "Cr2O3", "MnO", "Na2O"), "gt": ("Cr2O3", "MnO")}
# Two-pyroxene thermometers with Ca on M2 corrected for Na: cpx Na matters, opx Na hardly
_TWO_PX_OXIDES = {"opx": _BASE, "cpx": _BASE + ("Na2O",)}
_TWO_PX_MINOR = {"opx": ("Na2O",), "cpx": ("TiO2", "Cr2O3")}
_GT_CPX_OXIDES = {"cpx": _BASE, "gt": _BASE}
_GT_MINOR = {"gt": ("MnO",)}

THERMOMETERS = {
    # Two pyroxenes or opx alone
    "T98": Calibration(
        T_taylor98,
        "Taylor (1998) two-pyroxene",
        oxides={"opx": _BASE, "cpx": _BASE + ("Na2O",)},
        minor_oxides={"opx": ("TiO2", "Cr2O3", "Na2O"), "cpx": ("TiO2", "Cr2O3")},
    ),
    "BKN": Calibration(
        T_brey_kohler90, "Brey & Kohler (1990) two-pyroxene", oxides=_TWO_PX_OXIDES, minor_oxides=_TWO_PX_MINOR
    ),
    "BM85": Calibration(
        T_bertrand_mercier85,
        "Bertrand & Mercier (1985) two-pyroxene",
        oxides=_TWO_PX_OXIDES,
        minor_oxides=_TWO_PX_MINOR,
        options={"opx_fe": ("none", "PTEXL")},
    ),
    "WE77": Calibration(
        T_wells77,
        "Wells (1977) two-pyroxene",
        oxides={"opx": _BASE, "cpx": _BASE + ("Na2O",)},
        minor_oxides={"opx": ("TiO2", "Cr2O3", "MnO", "Na2O"), "cpx": ("TiO2", "Cr2O3", "MnO")},
        options={"sites": ("WB73", "PTEXL")},
    ),
    "BKNa": Calibration(
        T_brey_kohler90_na,
        "Brey & Kohler (1990) Na partitioning between opx and cpx",
        oxides={"opx": _BASE + ("Na2O",), "cpx": _BASE + ("Na2O",)},
    ),
    "BKCa": Calibration(
        T_brey_kohler90_ca_in_opx,
        "Brey & Kohler (1990) Ca-in-opx",
        oxides={"opx": _BASE},
        options={"correction": ("none", "NG10")},
    ),
    # cpx alone
    "NT00": Calibration(
        T_nimis_taylor00,
        "Nimis & Taylor (2000) enstatite-in-cpx",
        oxides={"cpx": _BASE + ("Cr2O3", "Na2O")},
        minor_oxides={"cpx": ("TiO2", "K2O")},
    ),
    "AS17": Calibration(
        T_ashchepkov17,
        "Nimis & Taylor (2000) enstatite-in-cpx, corrected by Ashchepkov et al. (2017)",
        oxides={"cpx": _BASE + ("Cr2O3", "Na2O")},
        minor_oxides={"cpx": ("TiO2", "K2O")},
        needs_wt=True,
    ),
    # Garnet-cpx Fe-Mg exchange
    "KR88": Calibration(T_krogh88, "Krogh (1988) garnet-cpx Fe-Mg", oxides=_GT_CPX_OXIDES, minor_oxides=_GT_MINOR),
    "EG79": Calibration(
        T_ellis_green79, "Ellis & Green (1979) garnet-cpx Fe-Mg", oxides=_GT_CPX_OXIDES, minor_oxides=_GT_MINOR
    ),
    "PO85": Calibration(T_powell85, "Powell (1985) garnet-cpx Fe-Mg", oxides=_GT_CPX_OXIDES, minor_oxides=_GT_MINOR),
    "AI94": Calibration(T_ai94, "Ai (1994) garnet-cpx Fe-Mg", oxides=_GT_CPX_OXIDES, minor_oxides=_GT_MINOR),
    "KR00": Calibration(
        T_krogh_ravna00, "Krogh Ravna (2000) garnet-cpx Fe-Mg", oxides=_GT_CPX_OXIDES, minor_oxides=_GT_MINOR
    ),
    "SU22gtcpx": Calibration(
        T_sudholz22_gt_cpx,
        "Sudholz et al. (2022) garnet-cpx Fe-Mg",
        oxides=_GT_CPX_OXIDES,
        minor_oxides={"cpx": ("TiO2", "Cr2O3", "Na2O")},
    ),
    # Garnet-opx Fe-Mg exchange
    "NG09": Calibration(
        T_nimis_grutter10,
        "Nimis & Grutter (2010) opx-garnet Fe-Mg",
        oxides={"opx": _BASE, "gt": _BASE},
        minor_oxides=_GT_MINOR,
    ),
    "HA84": Calibration(
        T_harley84, "Harley (1984) garnet-opx Fe-Mg", oxides={"opx": _BASE, "gt": _BASE}, minor_oxides=_GT_MINOR
    ),
    "SU22gtopx": Calibration(
        T_sudholz22_gt_opx, "Sudholz et al. (2022) garnet-opx Fe-Mg", oxides={"opx": _BASE, "gt": _BASE}
    ),
    "BK90gtopx": Calibration(
        T_brey_kohler90_gt_opx, "Brey & Kohler (1990) garnet-opx Fe-Mg", oxides={"opx": _BASE, "gt": _BASE}
    ),
    # Olivine
    "BK90olcpx": Calibration(
        T_brey_kohler90_ol_cpx, "Brey & Kohler (1990) olivine-cpx Fe-Mg", oxides={"ol": _OL, "cpx": _BASE}
    ),
    "KB90": Calibration(
        T_kohler_brey90, "Kohler & Brey (1990) Ca exchange, olivine-cpx", oxides={"ol": _OL + ("CaO",), "cpx": _BASE}
    ),
    "OW79": Calibration(
        T_oneill_wood79,
        "O'Neill & Wood (1979) olivine-garnet Fe-Mg",
        oxides={"ol": _OL, "gt": _BASE},
        minor_oxides=_GT_MINOR,
    ),
    "WZ07": Calibration(
        T_wu_zhao07,
        "Wu & Zhao (2007) olivine-garnet Fe-Mg",
        oxides={"ol": _OL, "gt": _BASE},
        minor_oxides=_GT_MINOR,
    ),
    "DH10Al": Calibration(
        T_de_hoog10_al,
        "De Hoog et al. (2010) Al-in-olivine",
        oxides={"ol": _OL + ("Al2O3",)},
        options={"version": ("eq12", "eq13")},
        option_oxides={("version", "eq12"): {"ol": ("Cr2O3",)}},
        needs_wt=True,
    ),
    "DH10Cr": Calibration(
        T_de_hoog10_cr, "De Hoog et al. (2010) Cr-in-olivine", oxides={"ol": _OL + ("Al2O3", "Cr2O3")}, needs_wt=True
    ),
    "BU17": Calibration(
        T_bussweiler17,
        "Bussweiler et al. (2017) Al-in-olivine",
        oxides={"ol": _OL + ("Al2O3",)},
        options={"version": ("eq3", "eq2")},
        option_oxides={("version", "eq2"): {"ol": ("Cr2O3",)}},
        needs_wt=True,
    ),
    # Spinel peridotites
    "OW87": Calibration(
        T_oneill_wall87,
        "O'Neill & Wall (1987) olivine-spinel Mg-Fe",
        oxides={"ol": _OL, "sp": _SP},
        minor_oxides={"sp": ("TiO2",)},
    ),
    "BA91": Calibration(
        T_ballhaus91,
        "Ballhaus et al. (1991) olivine-spinel Mg-Fe",
        oxides={"ol": _OL, "sp": _SP},
        minor_oxides={"sp": ("TiO2",)},
    ),
    "LI95": Calibration(T_li95, "Li et al. (1995) olivine-spinel Mg-Fe", oxides={"ol": _OL, "sp": _SP}),
    "WA08": Calibration(
        T_wan08,
        "Wan et al. (2008) Al exchange between olivine and spinel",
        oxides={"ol": _OL + ("Al2O3",), "sp": _SP},
        needs_wt=True,
    ),
    "CO14": Calibration(
        T_coogan14,
        "Coogan et al. (2014) Al exchange between olivine and spinel",
        oxides={"ol": _OL + ("Al2O3",), "sp": _SP},
        needs_wt=True,
    ),
    "WES91": Calibration(
        T_witt_eickschen_seck91,
        "Witt-Eickschen & Seck (1991) Al in opx with spinel",
        oxides={"opx": _BASE + ("Cr2O3",)},
        minor_oxides={"opx": ("TiO2", "Na2O")},
        options={"version": ("opx-sp", "Cr-Al-opx")},
        option_oxides={("version", "opx-sp"): {"ol": _OL, "sp": _SP}},
    ),
    # Ni in garnet
    "GR89Ni": Calibration(
        T_griffin89_ni,
        "Griffin et al. (1989) Ni-in-garnet",
        oxides={"gt": ("NiO",)},
        options={"ni_ol": ("assumed", "measured")},
        option_oxides={("ni_ol", "measured"): {"ol": ("NiO",)}},
        needs_wt=True,
    ),
    "RY96Ni": Calibration(T_ryan96_ni, "Ryan et al. (1996) Ni-in-garnet", oxides={"gt": ("NiO",)}, needs_wt=True),
    "CA99Ni": Calibration(
        T_canil99,
        "Canil (1999) Ni-in-garnet",
        oxides={"gt": ("NiO",)},
        options={"ni_ol": ("measured", "assumed")},
        option_oxides={("ni_ol", "measured"): {"ol": ("NiO",)}},
        needs_wt=True,
    ),
    "SU21Ni": Calibration(
        T_sudholz21_ni,
        "Sudholz et al. (2021) Ni-in-garnet",
        oxides={"gt": _BASE + ("Cr2O3", "NiO")},
        options={"ni_ol": ("assumed", "measured")},
        option_oxides={("ni_ol", "measured"): {"ol": ("NiO",)}},
        needs_wt=True,
    ),
    "NI24Ni": Calibration(
        T_nimis24_ni,
        "Nimis et al. (2024) Ni-in-garnet",
        oxides={"gt": ("NiO",)},
        options={"version": ("eq2", *_NIMIS24_GEOTHERMS), "ni_ol": ("assumed", "measured")},
        option_oxides={("ni_ol", "measured"): {"ol": ("NiO",)}},
        needs_wt=True,
    ),
    # Single minerals, independent of P
    "RY96Zn": Calibration(
        T_ryan96_zn,
        "Ryan et al. (1996) Zn-in-chromite",
        oxides={"sp": ("ZnO",)},
        options={"coefficients": ("RY96", "PTEXL")},
        needs_wt=True,
    ),
    "CR09": Calibration(
        T_creighton09,
        "Creighton (2009) single-crystal Mn-in-garnet",
        oxides={"gt": _BASE + ("Cr2O3", "MnO")},
        options={"cr_term": ("corrected", "printed")},
    ),
}

BAROMETERS = {
    # Al in opx coexisting with garnet
    "NG85": Calibration(
        P_nickel_green85,
        "Nickel & Green (1985) Al-in-opx",
        oxides=_AL_IN_OPX_OXIDES,
        minor_oxides=_AL_IN_OPX_MINOR,
        options=_NG85_OPTIONS,
    ),
    "T98": Calibration(
        P_nickel_green85,
        "Nickel & Green (1985) with Taylor (1998) Mg-Tschermak",
        oxides=_AL_IN_OPX_OXIDES,
        minor_oxides=_AL_IN_OPX_MINOR,
        options=_NG85_OPTIONS,
        preset={"mgts": "T98"},
    ),
    "BKN": Calibration(
        P_brey_kohler90,
        "Brey & Kohler (1990) Al-in-opx",
        oxides=_AL_IN_OPX_OXIDES,
        minor_oxides=_AL_IN_OPX_MINOR,
    ),
    "BBG08": Calibration(
        P_brey08, "Brey et al. (2008) Al-in-opx", oxides=_AL_IN_OPX_OXIDES, minor_oxides=_AL_IN_OPX_MINOR
    ),
    "MC74": Calibration(
        P_macgregor74,
        "MacGregor (1974) Al-in-opx (MAS), as fitted by Finnerty & Boyd (1984)",
        oxides={"opx": ("Al2O3",)},
        needs_wt=True,
    ),
    # Cr in cpx coexisting with garnet
    "SUD21": Calibration(
        P_sudholz21,
        "Sudholz et al. (2021) Cr-in-cpx",
        oxides={"cpx": _BASE + ("Cr2O3", "Na2O")},
        minor_oxides={"cpx": ("K2O",)},
    ),
    "NT00": Calibration(
        P_nimis_taylor00,
        "Nimis & Taylor (2000) Cr-in-cpx",
        oxides={"cpx": _BASE + ("Cr2O3", "Na2O")},
        minor_oxides={"cpx": ("K2O",)},
        options={"correction": ("none", "N20")},
    ),
    # Jd-Di exchange in cpx
    "AS17": Calibration(
        P_ashchepkov17,
        "Ashchepkov et al. (2011, corrected 2017) Jd-Di cpx",
        oxides={"cpx": _BASE + ("TiO2", "Cr2O3", "Na2O")},
        options={"kd": ("cations", "wt")},
        needs_wt=True,
    ),
    # Ca in olivine coexisting with cpx
    "KB90": Calibration(
        P_kohler_brey90, "Kohler & Brey (1990) Ca-in-olivine", oxides={"ol": _OL + ("CaO",), "cpx": _BASE}
    ),
    # Al in cpx coexisting with garnet (bimineralic eclogite and garnet pyroxenite)
    "BE15": Calibration(
        P_beyer15,
        "Beyer et al. (2015) garnet-cpx",
        oxides={"cpx": _BASE + ("Na2O",), "gt": _BASE},
        minor_oxides={"cpx": ("TiO2", "Cr2O3")},
    ),
    # Cr in garnet saturated with chromite
    "GR06": Calibration(
        P_grutter06,
        "Grutter et al. (2006) Cr-in-garnet (P38)",
        oxides={"gt": ("Cr2O3", "CaO")},
        options={"geotherm": tuple(str(q) for q in [38, *range(35, 38), *range(39, 49)])},
        needs_wt=True,
    ),
}
