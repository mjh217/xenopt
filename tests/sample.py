"""Test sample: garnet-spinel pyroxenite BM96-2MC from Bullen Merri (Lu et al., 2018), in the 2026 AusXenos database.

It is the one sample there with olivine, opx, cpx, garnet and spinel all analysed. Blank
oxides are NaN, as read_geochem gives them.
"""

import numpy as np
import pandas as pd

from xenopt import OXIDES, THERMOMETERS, cations_per_formula_unit
from xenopt.chemistry import C2K

NAMES = [o.name for o in OXIDES]
nan = np.nan
#           SiO2  TiO2  Al2O3 Cr2O3 FeOt  MnO   NiO   MgO   CaO   Na2O  K2O  ZnO   V2O3 P2O5
WT_PERC = {
    "opx": [53.3, 0.08, 4.63, 0.14, 11.6, 0.16, 0.12, 29.7, 0.75, 0.09, nan, nan, nan, nan],
    "cpx": [51.2, 0.37, 6.62, 0.28, 5.12, 0.09, 0.08, 14.5, 20.3, 1.55, nan, nan, nan, nan],
    "gt": [41.2, 0.08, 23.2, 0.45, 13.3, 0.41, nan, 17.4, 5.3, 0.0, nan, nan, nan, nan],
    "ol": [39.5, 0.0, 0.03, 0.0, 18.0, 0.15, 0.48, 41.6, 0.07, nan, nan, nan, nan, nan],
    "sp": [0.06, 0.1, 59.7, 3.74, 17.8, 0.08, 0.59, 17.3, 0.0, nan, nan, 0.23, nan, nan],
}

# A garnet-peridotite cpx that passes every cpx test: xenocryst FS03-C6R3-CPX-1 from Blanchetown
# (Sudholz et al., 2022), in the 2026 AusXenos database
PERIDOTITE_CPX = dict(zip(NAMES, [55.5, 0.14, 2.46, 1.31, 2.48, 0.08, 0.04, 16.03, 20.77, 1.69, 0.0, 0, 0, 0]))


def sample(**blank):
    """Oxide wt% and cations of the sample, with the given {mineral: [oxides]} also left blank."""
    wt = {m: pd.DataFrame([v], columns=NAMES) for m, v in WT_PERC.items()}
    for mineral, oxides in blank.items():
        wt[mineral].loc[0, oxides] = np.nan
    # Blank oxides count as zero in the formula unit
    return wt, {m: cations_per_formula_unit(w.fillna(0), m) for m, w in wt.items()}


WT, C = sample()


def calculate(registry, key, x, **options):
    cal = registry[key]
    kwargs = cal.resolve(options)
    if cal.needs_wt:
        kwargs["wt"] = {m: w.fillna(0) for m, w in WT.items()}
    return float(np.asarray(cal.func(C, x, **kwargs))[0])


def T_C(key, P=4.0, **options):
    """Temperature [C] of the sample at P [GPa]."""
    return calculate(THERMOMETERS, key, P, **options) - C2K
