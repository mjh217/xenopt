"""Oxide and mineral constants, and conversion of oxide wt% to cations per formula unit."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

C2K = 273.15  # [K] conversion between Celsius and Kelvin


@dataclass(frozen=True)
class Oxide:
    name: str  # column label in the geochemistry spreadsheet, e.g. "Al2O3"
    cation: str  # label of the corresponding cation, e.g. "Al"
    weight: float  # molecular weight [g/mol]
    n_cations: int  # cations per oxide formula
    n_oxygens: int  # oxygens per oxide formula


# All iron is reported as FeOt and treated as Fe2+.
OXIDES = (
    Oxide("SiO2", "Si", 60.084, 1, 2),
    Oxide("TiO2", "Ti", 79.865, 1, 2),
    Oxide("Al2O3", "Al", 101.961, 2, 3),
    Oxide("Cr2O3", "Cr", 151.990, 2, 3),
    Oxide("FeOt", "Fe", 71.846, 1, 1),
    Oxide("MnO", "Mn", 70.937, 1, 1),
    Oxide("NiO", "Ni", 74.689, 1, 1),
    Oxide("MgO", "Mg", 40.304, 1, 1),
    Oxide("CaO", "Ca", 56.079, 1, 1),
    Oxide("Na2O", "Na", 61.979, 2, 1),
    Oxide("K2O", "K", 94.196, 2, 1),
    Oxide("ZnO", "Zn", 81.379, 1, 1),
    Oxide("V2O3", "V", 149.881, 2, 3),
    Oxide("P2O5", "P", 141.944, 2, 5),
)
OXIDE_NAMES = tuple(o.name for o in OXIDES)
CATIONS = tuple(o.cation for o in OXIDES)
O_WEIGHT = 15.999  # [g/mol]


def element_ppm(wt_perc: pd.DataFrame, oxide: str) -> pd.Series:
    """Concentration [ppm by weight] of the cation of `oxide`, e.g. Al from Al2O3 wt%."""
    ox = OXIDES[OXIDE_NAMES.index(oxide)]
    cation_weight = (ox.weight - ox.n_oxygens * O_WEIGHT) / ox.n_cations
    return wt_perc[oxide] * 1e4 * ox.n_cations * cation_weight / ox.weight


@dataclass(frozen=True)
class Mineral:
    name: str
    prefix: str  # column prefix in the geochemistry spreadsheet, e.g. "CP" in "CPAl2O3"
    n_oxygens: int  # oxygens per formula unit
    n_cations: int  # ideal number of cations per formula unit


MINERALS = {
    "ol": Mineral("olivine", "OL", 4, 3),
    "opx": Mineral("orthopyroxene", "OP", 6, 4),
    "cpx": Mineral("clinopyroxene", "CP", 6, 4),
    "gt": Mineral("garnet", "GT", 12, 8),
    "sp": Mineral("spinel", "SP", 4, 3),
}


def cations_per_formula_unit(wt_perc: pd.DataFrame, mineral: str) -> pd.DataFrame:
    """Convert oxide wt% of one mineral into cations per formula unit.

    `wt_perc` has one row per sample and one column per oxide in OXIDE_NAMES.
    Blank (NaN) oxides are taken as zero. Rows for samples without an analysis
    of this mineral should be all NaN and stay NaN in the result. Columns of the
    result are named by cation (Si, Ti, ...).
    """
    weight = np.array([o.weight for o in OXIDES])
    n_cations = np.array([o.n_cations for o in OXIDES])
    n_oxygens = np.array([o.n_oxygens for o in OXIDES])

    values = wt_perc[list(OXIDE_NAMES)]
    moles = values.fillna(0).to_numpy(dtype=float) / weight
    oxygens = (moles * n_oxygens).sum(axis=1, keepdims=True)  # oxygens per 100 g of mineral
    with np.errstate(divide="ignore", invalid="ignore"):
        cations = moles * n_cations * MINERALS[mineral].n_oxygens / oxygens
    cations[values.isna().all(axis=1).to_numpy()] = np.nan

    return pd.DataFrame(cations, index=wt_perc.index, columns=list(CATIONS))
