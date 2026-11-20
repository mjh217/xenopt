"""Simultaneous solution of a thermometer and a barometer."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd

from .chemistry import C2K
from .thermobarometers import BAROMETERS, THERMOMETERS, Cations

CONVERGED = "converged"
MISSING_MINERALS = "missing minerals"
MISSING_OXIDES = "missing oxides"
NOT_FINITE = "not finite"
NOT_CONVERGED = "not converged"


def blank_oxides(wt_perc: Mapping[str, pd.DataFrame], oxides: Mapping[str, Iterable[str]]) -> pd.Series:
    """For each sample, the listed oxides left blank in an analysed mineral, e.g. "opx Na2O, opx TiO2"."""
    found = []
    for mineral, names in oxides.items():
        wt = wt_perc[mineral]
        analysed = wt.notna().any(axis=1)
        for name in names:
            found.append(pd.Series(np.where(analysed & wt[name].isna(), f"{mineral} {name}", ""), index=wt.index))
    if not found:
        return pd.Series("", index=next(iter(wt_perc.values())).index)
    return pd.concat(found, axis=1).apply(lambda row: ", ".join(s for s in row if s), axis=1)


def _blanks_as_zero(wt_perc: Mapping[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Blank oxides as zero in analysed minerals; minerals not analysed stay NaN."""
    return {m: wt.fillna(0).where(wt.notna().any(axis=1), np.nan) for m, wt in wt_perc.items()}


def _merge(*mappings: Mapping[str, Iterable[str]]) -> dict[str, list[str]]:
    merged: dict[str, list[str]] = {}
    for mapping in mappings:
        for mineral, names in mapping.items():
            merged.setdefault(mineral, [])
            merged[mineral] += [n for n in names if n not in merged[mineral]]
    return merged


def solve_pt(
    c: Cations,
    thermometer: str,
    barometer: str,
    thermometer_options: Mapping[str, str] | None = None,
    barometer_options: Mapping[str, str] | None = None,
    wt_perc: Mapping[str, pd.DataFrame] | None = None,
    check_blanks: bool = True,
    P0: float = 4.0,
    T0: float = 1200.0,
    T_tol: float = 0.1,
    P_tol: float = 0.01,
    max_iter: int = 100,
) -> pd.DataFrame:
    """Iterate a thermometer and barometer to a consistent P-T for every sample.

    Each iteration updates T at the current P, then P at the new T. A sample
    has converged once both |dT| < T_tol [K] and |dP| < P_tol [GPa]. P0 [GPa]
    and T0 [K] are the starting values. The options select variants of each
    calibration (see Calibration.options); by default the published versions.

    `wt_perc` is oxide wt% as returned by io.read_geochem, with blanks as NaN.
    Calibrations based on trace-element concentrations need it. If it is given
    and `check_blanks` is True, samples with a blank required oxide
    (Calibration.required_oxides) are not solved, and blank minor oxides
    (Calibration.minor_oxides), which are taken as zero, are listed in an
    "assumed_zero" column. With `check_blanks` False every blank in an analysed
    mineral is taken as zero.

    Returns a DataFrame indexed like the input with columns T_C, P_GPa, n_iter
    and status (and missing_oxides and assumed_zero). T_C and P_GPa are NaN
    unless status is "converged".
    """
    thermo, baro = THERMOMETERS[thermometer], BAROMETERS[barometer]
    thermo_kwargs, baro_kwargs = thermo.resolve(thermometer_options), baro.resolve(barometer_options)
    for cal, kwargs in [(thermo, thermo_kwargs), (baro, baro_kwargs)]:
        if cal.needs_wt:
            if wt_perc is None:
                raise ValueError(f"{cal.reference} needs oxide wt% (wt_perc)")
            kwargs["wt"] = wt_perc if check_blanks else _blanks_as_zero(wt_perc)
    check = wt_perc is not None and check_blanks
    required = sorted(set(thermo.required_oxides(thermometer_options)) | set(baro.required_oxides(barometer_options)))
    index = c[required[0]].index

    has_minerals = np.logical_and.reduce([c[m].notna().all(axis=1).to_numpy() for m in required])
    status = np.where(has_minerals, NOT_CONVERGED, MISSING_MINERALS).astype(object)
    if check:
        oxides = _merge(thermo.required_oxides(thermometer_options), baro.required_oxides(barometer_options))
        missing = blank_oxides(wt_perc, oxides)
        status[has_minerals & (missing.to_numpy() != "")] = MISSING_OXIDES
    active = status == NOT_CONVERGED
    n_iter = np.zeros(len(index), dtype=int)
    T = np.full(len(index), T0, dtype=float)
    P = np.full(len(index), P0, dtype=float)

    # Every sample is computed each iteration but only active ones are updated.
    # Logs of non-positive activities give NaN, which is caught as not finite.
    with np.errstate(divide="ignore", invalid="ignore"):
        for iteration in range(1, max_iter + 1):
            if not active.any():
                break
            T_new = np.asarray(thermo.func(c, P, **thermo_kwargs), dtype=float)
            P_new = np.asarray(baro.func(c, T_new, **baro_kwargs), dtype=float)

            finite = np.isfinite(T_new) & np.isfinite(P_new)
            converged = (np.abs(T_new - T) < T_tol) & (np.abs(P_new - P) < P_tol)

            T = np.where(active, T_new, T)
            P = np.where(active, P_new, P)
            n_iter[active] = iteration
            status[active & ~finite] = NOT_FINITE
            status[active & finite & converged] = CONVERGED
            active &= finite & ~converged

    ok = status == CONVERGED
    result = pd.DataFrame(
        {
            "T_C": np.where(ok, T - C2K, np.nan),
            "P_GPa": np.where(ok, P, np.nan),
            "n_iter": n_iter,
            "status": status,
        },
        index=index,
    )
    if check:
        assumed = blank_oxides(wt_perc, _merge(thermo.minor_oxides, baro.minor_oxides))
        result["missing_oxides"] = missing.where(status == MISSING_OXIDES, "")
        result["assumed_zero"] = assumed.where(ok, "")
    return result
