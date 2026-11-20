import argparse

import numpy as np
import pandas as pd
import pytest

from xenopt.chemistry import MINERALS, OXIDE_NAMES
from xenopt.cli import main, parse_combo
from xenopt.io import dat_tag, read_geochem, write_dat_pt, write_dat_tests

from .sample import WT_PERC


def write_csv(path, rows):
    """Geochemistry CSV in the layout of the 2026 AusXenos database."""
    columns = ["SAMPLE", "LOCALITY", "SOURCE"]
    for m in MINERALS.values():
        columns += [m.prefix + ox for ox in OXIDE_NAMES] + [m.prefix + "Tot"]
    columns += ["", "OL_Check"]
    table = pd.DataFrame(rows).reindex(columns=columns)
    table.to_csv(path, index=False)
    return path


def sample_row(name):
    row = {"SAMPLE": name, "LOCALITY": "Bullen Merri", "SOURCE": "Lu et al. (2018)", "OL_Check": 0}
    for m, wt in WT_PERC.items():
        row.update({MINERALS[m].prefix + ox: v for ox, v in zip(OXIDE_NAMES, wt)})
    return row


def test_read_geochem(tmp_path):
    path = write_csv(tmp_path / "in.csv", [sample_row("A"), {"SAMPLE": "B", "CPSiO2": 54.0, "CPCaO": 20.0}])
    metadata, wt = read_geochem(path)
    assert list(metadata.columns) == ["SAMPLE", "LOCALITY", "SOURCE"]
    assert list(wt["cpx"].columns) == list(OXIDE_NAMES)
    # Blank oxides stay NaN, so they can be told apart from reported zeros
    assert np.isnan(wt["cpx"].loc[1, "Al2O3"]) and wt["cpx"].loc[1, "SiO2"] == 54.0
    assert wt["opx"].loc[1].isna().all()
    assert wt["ol"].loc[1].isna().all()


def test_read_geochem_rejects_text(tmp_path):
    path = write_csv(tmp_path / "in.csv", [{**sample_row("A"), "CPMgO": "15.6x"}])
    with pytest.raises(ValueError, match="CPMgO"):
        read_geochem(path)


def test_dat_formats(tmp_path):
    pt = pd.DataFrame({"P_GPa": [4.1496, np.nan], "T_C": [948.54, np.nan]})
    write_dat_pt(tmp_path / "pt.dat", pt)
    assert (tmp_path / "pt.dat").read_text() == " 4.150  949.\n   NaN   NaN\n"
    write_dat_tests(tmp_path / "t.dat", pd.DataFrame({"a": [True, False], "b": [False, True]}))
    assert (tmp_path / "t.dat").read_text() == " 1 0\n 0 1\n"
    assert dat_tag("T98", "NG85") == "PNG85_TT98"
    assert dat_tag("NT00", "SUD21") == "PS21_TNT00"


def test_parse_combo():
    combo = parse_combo("T98:NG85,mgts=T98,mgts_in=all")
    assert (combo.thermometer, combo.barometer) == ("T98", "NG85")
    assert combo.barometer_options == {"mgts": "T98", "mgts_in": "all"}
    assert combo.name == "T98_NG85-T98-all"
    assert combo.dat_tag == "PNG85-T98-all_TT98"
    assert parse_combo("NT00:SUD21").dat_tag == "PS21_TNT00"
    for bad in ["T98", "XX:NG85", "T98:NG85,mgts", "T98:NG85,mgts=XX", "T98:SUD21,mgts=CA91"]:
        with pytest.raises(argparse.ArgumentTypeError):
            parse_combo(bad)


def test_cli(tmp_path):
    path = write_csv(tmp_path / "in.csv", [sample_row("A"), {"SAMPLE": "B", "LOCALITY": "Elsewhere"}])
    pd.DataFrame({"LOCALITY": ["Bullen Merri"], "SHORT FORM": ["bullen"]}).to_csv(tmp_path / "locs.csv", index=False)
    out = tmp_path / "out" / "results.csv"
    assert (
        main(
            [
                str(path),
                "-o",
                str(out),
                "--localities",
                str(tmp_path / "locs.csv"),
                "--dat-dir",
                str(tmp_path / "dat"),
                "--cations",
                str(tmp_path / "cations.csv"),
            ]
        )
        == 0
    )

    results = pd.read_csv(out)
    assert results["SHORT FORM"].tolist() == ["bullen", "Elsewhere"]
    assert results["status_T98_NG85"].tolist() == ["converged", "missing minerals"]
    assert results["T_C_T98_NG85"][0] == pytest.approx(907.8, abs=0.1)
    assert results["has_gt"].tolist() == [1, 0]
    # A converged at 1.26 GPa, too low for garnet lherzolite; B has no result, so no warning
    assert "low P (ok if pyroxenite)" in results["warnings_T98_NG85"][0]
    assert pd.isna(results["warnings_T98_NG85"][1])
    assert "ng10_cpx_opx_T98_NG85" in results and "ziberna16_T_NT00_SUD21" in results

    assert (tmp_path / "dat" / "sample_ids.txt").read_text() == "A bullen\nB Elsewhere\n"
    assert (tmp_path / "dat" / "PT_results_PNG85_TT98.dat").read_text().splitlines()[1] == "   NaN   NaN"
    assert len((tmp_path / "dat" / "suitability_tests_PS21_TNT00.dat").read_text().split("\n")[0].split()) == 12
