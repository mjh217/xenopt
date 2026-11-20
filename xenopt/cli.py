"""Command line interface: python -m xenopt INPUT.csv [options]"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .chemistry import MINERALS, cations_per_formula_unit
from .io import dat_tag, read_geochem, read_localities, write_dat
from .solve import CONVERGED, solve_pt
from .suitability import cpx_tests, pair_tests, warning_text
from .thermobarometers import BAROMETERS, THERMOMETERS

DEFAULT_COMBOS = ["T98:NG85", "NT00:SUD21"]

EPILOG = """\
Options of a calibration follow its name, separated by commas, e.g.
  -c T98:NG85,mgts=CA91              NG85 with the Carswell (1991) Mg-Tschermak
  -c T98:NG85,mgts=T98,mgts_in=all  NG85 with Taylor (1998) applied everywhere

Hoggard et al. (2026) used
  -c T98:NG85 -c NT00:SUD21 --cpx-cation-range 3.97 4.03 --blanks-as-zero
and rejected cpx failing cpx_oxide_total, cpx_cation_total, nimis98_al_mg, ziberna16_cr2o3,
ziberna16_T, sudholz21_na or sudholz21_cacrts (see the README).
"""


@dataclass(frozen=True)
class Combo:
    thermometer: str
    barometer: str
    thermometer_options: dict = field(default_factory=dict)
    barometer_options: dict = field(default_factory=dict)

    @staticmethod
    def _suffix(options: dict) -> str:
        return "".join(f"-{value}" for value in options.values())

    @property
    def name(self) -> str:
        """Column suffix, e.g. "T98_NG85" or "T98_NG85-CA91"."""
        thermometer = self.thermometer + self._suffix(self.thermometer_options)
        return f"{thermometer}_{self.barometer}{self._suffix(self.barometer_options)}"

    @property
    def dat_tag(self) -> str:
        return dat_tag(
            self.thermometer,
            self.barometer,
            self._suffix(self.thermometer_options),
            self._suffix(self.barometer_options),
        )


def _parse_calibration(text: str, registry: dict, kind: str) -> tuple[str, dict]:
    name, *pairs = text.split(",")
    if name not in registry:
        raise argparse.ArgumentTypeError(f"unknown {kind} '{name}' (choose from {', '.join(registry)})")
    options = {}
    for pair in pairs:
        key, sep, value = pair.partition("=")
        if not sep:
            raise argparse.ArgumentTypeError(f"{kind} option '{pair}' should be NAME=VALUE")
        options[key] = value
    try:
        registry[name].resolve(options)
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"{kind} {name}: {err}") from None
    return name, options


def parse_combo(text: str) -> Combo:
    thermometer_text, sep, barometer_text = text.partition(":")
    if not sep:
        raise argparse.ArgumentTypeError(f"'{text}' is not THERMOMETER:BAROMETER")
    thermometer, thermometer_options = _parse_calibration(thermometer_text, THERMOMETERS, "thermometer")
    barometer, barometer_options = _parse_calibration(barometer_text, BAROMETERS, "barometer")
    return Combo(thermometer, barometer, thermometer_options, barometer_options)


def list_calibrations() -> str:
    lines = []
    for title, registry in [("Thermometers", THERMOMETERS), ("Barometers", BAROMETERS)]:
        lines.append(f"{title}:")
        for key, cal in registry.items():
            lines.append(f"  {key:6s} {cal.reference}")
            for mineral, oxides in cal.required_oxides().items():
                minor = " ".join(cal.minor_oxides.get(mineral, ()))
                line = f"           {mineral}: requires {' '.join(oxides)}"
                lines.append(line + (f"; minor {minor}" if minor else ""))
            defaults = cal.resolve()
            for option, values in cal.options.items():
                choices = [f"{v} (default)" if v == defaults[option] else v for v in values]
                lines.append(f"           {option}: {', '.join(choices)}")
            for (option, value), extra in cal.option_oxides.items():
                if defaults[option] != value:
                    needs = "; ".join(f"{m} {' '.join(oxides)}" for m, oxides in extra.items())
                    lines.append(f"           {option}={value} also requires {needs}")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="xenopt",
        description="P-T estimates for mantle xenoliths from mineral geochemistry.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", nargs="?", help="geochemistry CSV, one row per sample")
    parser.add_argument("-o", "--output", default="xenopt_results.csv", help="results CSV (default: %(default)s)")
    parser.add_argument(
        "-c",
        "--combo",
        action="append",
        type=parse_combo,
        metavar="THERMO:BARO",
        help=f"thermometer:barometer pair to solve; repeat for several (default: {' '.join(DEFAULT_COMBOS)})",
    )
    parser.add_argument("--localities", help="localities CSV, to add the short locality names")
    parser.add_argument("--cations", metavar="CSV", help="also write cations per formula unit to this file")
    parser.add_argument(
        "--dat-dir", metavar="DIR", help="also write sample_ids.txt and fixed-width .dat files of P-T and test results"
    )
    parser.add_argument("--max-iter", type=int, default=100, help="maximum P-T iterations (default: %(default)s)")
    parser.add_argument(
        "--blanks-as-zero",
        action="store_true",
        help="take every blank oxide as zero without checking which oxides each calibration needs",
    )
    tests = parser.add_argument_group("cpx suitability test thresholds (defaults are the published values)")
    tests.add_argument(
        "--cpx-cation-range",
        nargs=2,
        type=float,
        default=(3.98, 4.02),
        metavar=("MIN", "MAX"),
        help="allowed cations per 6 O (default: %(default)s)",
    )
    tests.add_argument(
        "--cpx-cacrts-min", type=float, default=0.005, help="minimum a_CaCrTs / Cr# (default: %(default)s)"
    )
    tests.add_argument("--cpx-ca-mg-max", type=float, default=0.5, help="maximum Ca / (Ca + Mg) (default: %(default)s)")
    tests.add_argument(
        "--cpx-cr2o3-min", type=float, default=0.5, help="Cr2O3 [wt%%] must exceed this (default: %(default)s)"
    )
    tests.add_argument("--cpx-cr2o3-min-inclusive", action="store_true", help="accept Cr2O3 equal to --cpx-cr2o3-min")
    tests.add_argument("--T-min", type=float, default=600, help="minimum temperature [C] (default: %(default)s)")
    ng10 = parser.add_argument_group("opx-garnet test thresholds (Nimis & Grutter 2010)")
    ng10.add_argument(
        "--ng10-P-min", type=float, default=1.5, help="minimum P [GPa] for lherzolite (default: %(default)s)"
    )
    ng10.add_argument(
        "--ng10-P-min-pyroxenite",
        type=float,
        default=1.0,
        help="minimum P [GPa] for pyroxenite; between the two the warning says 'ok if pyroxenite' "
        "(default: %(default)s)",
    )
    ng10.add_argument(
        "--ng10-T-min", type=float, default=700, help="minimum T98 temperature [C] (default: %(default)s)"
    )
    ng10.add_argument(
        "--ng10-cpx-opx-tol",
        nargs=3,
        type=float,
        default=(90, 70, 50),
        metavar=("LOW", "MID", "HIGH"),
        help="max |T_T98 - T_Ca-in-opx| [C] below 900, 900-1200 and above 1200 C; PTEXL uses 100 80 60 "
        "(default: %(default)s)",
    )
    ng10.add_argument(
        "--ng10-gt-opx-tol", type=float, default=70, help="max |T_T98 - T_NG09| [C] (default: %(default)s)"
    )
    ng10.add_argument(
        "--gt-cpx-tol",
        type=float,
        help="also flag |T_T98 - T_Krogh88| above this [C], as PTEXL does with 90 (default: off)",
    )
    parser.add_argument("--list", action="store_true", help="list thermometers, barometers and their options")
    args = parser.parse_args(argv)

    if args.list:
        print(list_calibrations())
        return 0
    if args.input is None:
        parser.error("the input CSV is required")
    combos = args.combo or [parse_combo(c) for c in DEFAULT_COMBOS]

    metadata, wt_perc = read_geochem(args.input)
    n = len(metadata)
    print(f"Found {n} samples in {args.input}")
    if "SAMPLE" in metadata:
        duplicates = metadata["SAMPLE"][metadata["SAMPLE"].duplicated(keep=False)]
        if len(duplicates):
            print(f"Warning: duplicate sample names: {', '.join(sorted(set(duplicates.astype(str))))}")
    if args.localities:
        short = metadata["LOCALITY"].map(read_localities(args.localities))
        unmatched = sorted(set(metadata["LOCALITY"][short.isna()].dropna()))
        if unmatched:
            print(f"Warning: localities without a short form (kept as is): {', '.join(unmatched)}")
        metadata = metadata.assign(**{"SHORT FORM": short.fillna(metadata["LOCALITY"])})

    cations = {m: cations_per_formula_unit(wt, m) for m, wt in wt_perc.items()}
    for m in MINERALS:
        print(f"  {cations[m].notna().all(axis=1).sum():6d} have {m} analyses")

    tests = cpx_tests(
        wt_perc["cpx"],
        cations["cpx"],
        cation_range=tuple(args.cpx_cation_range),
        cacrts_cr_min=args.cpx_cacrts_min,
        ca_mg_max=args.cpx_ca_mg_max,
        cr2o3_min=args.cpx_cr2o3_min,
        cr2o3_min_inclusive=args.cpx_cr2o3_min_inclusive,
    )
    print("cpx suitability tests passed:")
    for name, passed in tests.items():
        print(f"  {passed.sum():6d}  {name}")

    output = [metadata]
    output.append(pd.DataFrame({f"has_{m}": cations[m].notna().all(axis=1) for m in MINERALS}))
    results = []
    for combo in combos:
        pt = solve_pt(
            cations,
            combo.thermometer,
            combo.barometer,
            combo.thermometer_options,
            combo.barometer_options,
            wt_perc=wt_perc,
            check_blanks=not args.blanks_as_zero,
            max_iter=args.max_iter,
        )
        results.append((combo.dat_tag, pt))
        minerals = set(THERMOMETERS[combo.thermometer].minerals) | set(BAROMETERS[combo.barometer].minerals)
        combo_tests = pair_tests(
            minerals,
            wt_perc,
            cations,
            pt,
            tests,
            T_min=args.T_min,
            P_min=args.ng10_P_min,
            P_min_pyroxenite=args.ng10_P_min_pyroxenite,
            T_ta98_min=args.ng10_T_min,
            cpx_opx_tol=tuple(args.ng10_cpx_opx_tol),
            gt_opx_tol=args.ng10_gt_opx_tol,
            gt_cpx_tol=args.gt_cpx_tol,
        )
        combo_warnings = warning_text(combo_tests, pt)

        counts = pt["status"].value_counts()
        n_warned = ((pt["status"] == CONVERGED) & (combo_warnings != "")).sum()
        print(
            f"{combo.name}: "
            + ", ".join(f"{counts.get(s, 0)} {s}" for s in [CONVERGED] + [s for s in counts.index if s != CONVERGED])
            + f"; {n_warned} converged with warnings"
        )
        # Tests already in the output (the cpx tests) are not repeated for each pair
        columns = pd.concat([pt, combo_tests.drop(columns=tests.columns, errors="ignore")], axis=1)
        columns["warnings"] = combo_warnings
        output.append(columns.add_suffix(f"_{combo.name}"))
    output.append(tests)

    table = pd.concat(output, axis=1)
    bools = table.select_dtypes(bool).columns
    table[bools] = table[bools].astype(int)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False, float_format="%.5f")
    print(f"Wrote {args.output}")

    if args.cations:
        Path(args.cations).parent.mkdir(parents=True, exist_ok=True)
        frames = [metadata] + [cations[m].add_prefix(MINERALS[m].prefix) for m in MINERALS]
        pd.concat(frames, axis=1).to_csv(args.cations, index=False, float_format="%.6f")
        print(f"Wrote {args.cations}")

    if args.dat_dir:
        written = write_dat(args.dat_dir, metadata, results, tests, args.T_min)
        print(f"Wrote {len(written)} files to {args.dat_dir} ({', '.join(tag for tag, _ in results)})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
