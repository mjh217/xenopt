"""Reading the geochemistry spreadsheet and writing results."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .chemistry import MINERALS, OXIDE_NAMES
from .suitability import temperature_test


def read_geochem(path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Read a geochemistry CSV with one row per sample.

    Oxide columns are found by header, e.g. "CPAl2O3" for cpx Al2O3 (prefixes
    are set in chemistry.MINERALS). Metadata are the columns before the first
    oxide column; total and check columns are ignored.

    Returns (metadata, wt_perc). wt_perc maps each mineral to a DataFrame of
    oxide wt%. Blank oxides stay NaN, so they can be told apart from reported
    zeros; rows of samples without an analysis of that mineral (all blank or
    zero) are entirely NaN.
    """
    table = pd.read_csv(path)
    oxide_columns = [m.prefix + ox for m in MINERALS.values() for ox in OXIDE_NAMES]
    found = [table.columns.get_loc(col) for col in oxide_columns if col in table.columns]
    if not found:
        raise ValueError(f"No oxide columns (e.g. 'CPAl2O3') found in {path}")
    metadata = table.iloc[:, : min(found)]

    wt_perc = {}
    for key, mineral in MINERALS.items():
        columns = [mineral.prefix + ox for ox in OXIDE_NAMES]
        missing = [col for col in columns if col not in table.columns]
        if missing:
            warnings.warn(f"Columns missing from {path}, taken as zero: {', '.join(missing)}", stacklevel=2)
        block = table.reindex(columns=columns)
        for col in block.columns:
            try:
                block[col] = pd.to_numeric(block[col])
            except ValueError as err:
                raise ValueError(f"Non-numeric value in column {col} of {path}: {err}") from None
        block.columns = list(OXIDE_NAMES)
        block.loc[block.sum(axis=1) <= 0] = np.nan
        wt_perc[key] = block

    return metadata, wt_perc


def read_localities(path) -> dict[str, str]:
    """Map locality name to its short form, from the localities CSV."""
    table = pd.read_csv(path)
    return dict(zip(table["LOCALITY"], table["SHORT FORM"]))


# --------------------------------------------------------------------------
# Fixed-width .dat output, e.g. for GMT plotting scripts
# --------------------------------------------------------------------------

# Shorter codes used in the .dat file names, where they differ
_DAT_CODES = {"SUD21": "S21"}


def dat_tag(thermometer: str, barometer: str, thermometer_suffix: str = "", barometer_suffix: str = "") -> str:
    """File-name tag such as "PNG85_TT98" for the .dat outputs.

    The suffixes name any non-default options, e.g. "-CA91" for "PNG85-CA91_TT98".
    """
    barometer = _DAT_CODES.get(barometer, barometer) + barometer_suffix
    thermometer = _DAT_CODES.get(thermometer, thermometer) + thermometer_suffix
    return f"P{barometer}_T{thermometer}"


def write_dat_pt(path, pt: pd.DataFrame) -> None:
    """PT_results*.dat: P [GPa] and T [C] as fixed-width columns (6.3f, 5.0f)."""
    with open(path, "w") as f:
        for P, T in zip(pt["P_GPa"], pt["T_C"]):
            if np.isfinite(P) and np.isfinite(T):
                f.write(f"{P:6.3f}{T:5.0f}.\n")
            else:
                f.write("   NaN   NaN\n")


def write_dat_tests(path, tests: pd.DataFrame) -> None:
    """suitability_tests*.dat: one 0/1 column per test, two characters wide."""
    with open(path, "w") as f:
        for row in tests.to_numpy(dtype=int):
            f.write("".join(f"{v:2d}" for v in row) + "\n")


def write_dat_ids(path, samples: pd.Series, localities: pd.Series) -> None:
    """sample_ids.txt: sample name and locality, space separated."""
    with open(path, "w") as f:
        for sample, locality in zip(samples, localities):
            f.write(f"{sample} {locality}\n")


def write_dat(
    directory, metadata: pd.DataFrame, results: list[tuple[str, pd.DataFrame]], tests: pd.DataFrame, T_min: float = 600
) -> list[Path]:
    """Write sample_ids.txt and, per combination, PT_results_<tag>.dat and suitability_tests_<tag>.dat.

    `results` holds (tag, output of solve_pt) for each combination; see dat_tag.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    locality = metadata["SHORT FORM"] if "SHORT FORM" in metadata else metadata["LOCALITY"]
    written = [directory / "sample_ids.txt"]
    write_dat_ids(written[0], metadata["SAMPLE"], locality)

    for tag, pt in results:
        combo_tests = tests.copy()
        # The temperature test goes before the last cpx test (ptexl_na_k)
        combo_tests.insert(len(tests.columns) - 1, "ziberna16_T", temperature_test(pt["T_C"], T_min))
        written += [directory / f"PT_results_{tag}.dat", directory / f"suitability_tests_{tag}.dat"]
        write_dat_pt(written[-2], pt)
        write_dat_tests(written[-1], combo_tests)
    return written
