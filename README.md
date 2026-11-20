# xenoPT

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24386935.svg)](https://doi.org/10.5281/zenodo.24386935)

Compute equilibration pressure-temperature (P-T) estimates using published thermobarometers for mantle xenoliths from mineral major-element geochemistry. This code was developed as part of Geoscience Australia's [*Resourcing Australia's Prosperity*](https://www.ga.gov.au/scientific-topics/resourcing-australias-prosperity) initiative (eCat record 151570).

## Usage

To estimate P-T for every sample in `geochem.csv` using the Nimis & Taylor (2000) enstatite-in-cpx thermometer
and the Sudholz et al. (2021) Cr-in-cpx barometer, run:

```bash
python -m xenopt geochem.csv -c NT00:SUD21
```

The results are written to `xenopt_results.csv`. Options:

- `-o FILE` names the results file.
- `-c THERMO:BARO` picks a thermometer/barometer pair; repeat it to solve several pairs in one
  run (default `T98:NG85` and `NT00:SUD21`). Options within a given thermobarometer follow its name after
  commas, e.g. `-c T98:NG85,mgts=CA91`. `python -m xenopt --list` shows everything that is available.
- `--localities FILE` adds the `SHORT FORM` locality name to outputs.
- `--cations FILE` writes a file of cations per formula unit.
- `--dat-dir DIR` writes fixed-width text files of outputs for easy use in other programs like GMT: `sample_ids.txt` (sample and
  locality), `PT_results_<tag>.dat` (P [GPa] and T [°C]) and `suitability_tests_<tag>.dat` (0/1
  for each [cpx test](#clinopyroxene-suitability-tests): `cpx_present`, `cpx_oxide_total`,
  `cpx_cation_total`, `nimis98_cr_al`, `nimis98_al_mg`, `sudholz21_na`, `sudholz21_cr_no`,
  `sudholz21_cacrts`, `ziberna16_cr2o3`, `ziberna16_ca_mg`, `ziberna16_T`, `ptexl_na_k`). Tags
  read `P<barometer>_T<thermometer>`, e.g. `PNG85_TT98`, `PS21_TNT00` (SUD21 shortened to S21)
  or `PNG85-CA91_TT98` with options.
- `--cpx-cation-range`, `--cpx-cacrts-min`, `--cpx-ca-mg-max`, `--cpx-cr2o3-min`,
  `--cpx-cr2o3-min-inclusive` and `--T-min` change the cpx suitability thresholds, and the
  `--ng10-*` options the opx-garnet ones. `--gt-cpx-tol 90` adds the garnet-cpx check from PTEXL (off by
  default).
- `--blanks-as-zero` interprets every blank oxide in the sample as zero, skipping the required-oxides check.

Or `pip install -e .` once, then run `xenopt ...` from anywhere.

### Input

A CSV with one sample per row in the layout of the [*AusXenos*](https://github.com/mjh217/ausxenos) database.
Oxide columns are found by header: mineral prefix (`OL`, `OP`, `CP`, `GT`, `SP`) plus oxide (`SiO2 TiO2 Al2O3 Cr2O3 FeOt MnO NiO MgO CaO Na2O K2O ZnO V2O3 P2O5`),
e.g. `CPAl2O3` is aluminium oxide wt% in clinopyroxene.
All Fe is treated as FeO.
Columns before the first oxide column are carried through as metadata.
A mineral with no oxides at all is treated as not analysed.
Blank oxides are kept apart from reported zeros (see [Required oxides](#required-oxides)).

### Output

The results file starts with the metadata columns copied from the input (e.g. `SAMPLE`, `LOCALITY`,
`SOURCE`, `SHORT FORM` if `--localities` is given). Next come `has_ol`, `has_opx`,
`has_cpx`, `has_gt` and `has_sp` which is 1 if the sample has an analysis of that mineral and 0 if not.
Then, for each selected thermobarometer pair, it lists:

- `T_C_<pair>`, `P_GPa_<pair>`, `n_iter_<pair>`: equilibration P-T conditions (NaN if not converged).
- `status_<pair>`: `converged`, `missing minerals`, `missing oxides`, `not finite` (e.g. log of a non-positive activity) or `not converged`.
- `missing_oxides_<pair>`: required oxides that were blank (if any).
- `assumed_zero_<pair>`: blank minor oxides that were taken as zero for the purposes of calculation.
- Results of any tests that rely on the P-T result (1/0 = pass/fail): `ziberna16_T_<pair>` for single-cpx thermobarometyr, and
  all the [opx-garnet tests](#orthopyroxene-garnet-tests) for thermobarometer pairs using opx, cpx and garnet.
- `warnings_<pair>`: indicator of why a result might not be trustworthy, e.g.
  `pyroxene temperatures disagree; low T (accuracy not assured); assumed zero: opx TiO2`. Blank
  means no concerns for this sample and thermobarometer pairing.

Last come the [cpx suitability tests](#clinopyroxene-suitability-tests). They use only the cpx composition,
so they are run on every sample, irrespective of which thermobarometers are chosen, and are written once rather than for
each thermobarometer pair (samples without cpx get 0). They only appear in the warnings of pairs that use cpx alone
(e.g. `NT00:SUD21`). A failed test adds a warning but does not automatically remove a P-T result.

### Warnings

Each failed test is described by a short phrase, separated by `; `. The test columns give the
pass/fail details, and the [cpx](#clinopyroxene-suitability-tests) and [opx-garnet](#orthopyroxene-garnet-tests)
sections explain each test.

| Warning | Test | | Warning | Test |
|---|---|---|---|---|
| `cpx oxide total out of range` | `cpx_oxide_total` | | `poor opx` / `cpx` / `garnet analysis` | `ng10_opx_quality`, `ng10_cpx_quality`, `ng10_gt_quality` |
| `cpx cation total out of range` | `cpx_cation_total` | | `low-Cr, low-Mg garnet` | `ng10_gt_type` |
| `outside garnet-peridotite field` | `nimis98_cr_al` | | `very low T` | `ng10_T_min` |
| `low Al2O3 for MgO` | `nimis98_al_mg` | | `low P (ok if pyroxenite)` | `ng10_P_min` |
| `high Na` | `sudholz21_na` | | `very low P` | `ng10_P_min_pyroxenite` |
| `Cr# out of range` | `sudholz21_cr_no` | | `pyroxene temperatures disagree` | `ng10_cpx_opx` |
| `low CaCrTs activity` | `sudholz21_cacrts` | | `garnet-opx temperature disagrees` | `ng10_gt_opx` |
| `Cr2O3 out of range` | `ziberna16_cr2o3` | | `garnet-cpx temperature disagrees` | `ptexl_gt_cpx` |
| `high Ca/(Ca+Mg)` | `ziberna16_ca_mg` | | `low T (accuracy not assured)` | `ng10_T_accuracy` |
| `low Na+K` | `ptexl_na_k` | | `high P (accuracy not assured)` | `ng10_P_accuracy` |
| `low T` | `ziberna16_T` | | | |

`very low P` and `very low T` replace `low P` and `low T (accuracy not assured)` when both apply.
Blank minor oxides taken as zero are listed as `assumed zero: ...`. Samples with no P-T result due to
a required oxide being blank are labelled `missing oxides: ...`.

## Available thermobarometers

Default thermobarometer equations are the original published versions, with obvious misprints corrected.
Later modifications by other authors are listed as options.
Every formula was checked against its original paper, cited with the equation number in the docstring of its function in
`xenopt/thermobarometers.py`. All Fe is Fe<sup>2+</sup>, except in spinel (see below).
`python -m xenopt --list` shows the keys, required oxides and options.

**Thermometers**

| Key | Calibration | Minerals | Options |
|---|---|---|---|
| T98 | Taylor (1998) two-pyroxene; Taylor calls it TA97 and Nimis & Grutter (2010) TA98 | opx, cpx | |
| BKN | Brey & Kohler (1990) two-pyroxene, eq. (9) | opx, cpx | |
| BM85 | Bertrand & Mercier (1985) two-pyroxene | opx, cpx | `opx_fe=PTEXL` |
| WE77 | Wells (1977) two-pyroxene, eq. (5) | opx, cpx | `sites=PTEXL` |
| BKNa | Brey & Kohler (1990) Na partitioning between opx and cpx, eq. (11) | opx, cpx | |
| BKCa | Brey & Kohler (1990) Ca-in-opx, eq. (10) | opx | `correction=NG10` |
| NT00 | Nimis & Taylor (2000) enstatite-in-cpx | cpx | |
| AS17 | NT00 with the correction of Ashchepkov et al. (2017, p. 783) | cpx | |
| KR88 | Krogh (1988) garnet-cpx Fe-Mg | cpx, gt | |
| EG79 | Ellis & Green (1979) garnet-cpx Fe-Mg, eq. (9) | cpx, gt | |
| PO85 | Powell (1985) garnet-cpx Fe-Mg, eq. (17) | cpx, gt | |
| AI94 | Ai (1994) garnet-cpx Fe-Mg, eq. (3) | cpx, gt | |
| KR00 | Krogh Ravna (2000) garnet-cpx Fe-Mg | cpx, gt | |
| SU22gtcpx | Sudholz et al. (2022) garnet-cpx Fe-Mg, eq. (1) | cpx, gt | |
| NG09 | Nimis & Grutter (2010) opx-garnet Fe-Mg, eq. (4) | opx, gt | |
| HA84 | Harley (1984) garnet-opx Fe-Mg, eq. (11) | opx, gt | |
| SU22gtopx | Sudholz et al. (2022) garnet-opx Fe-Mg, eq. (2) | opx, gt | |
| BK90gtopx | Brey & Kohler (1990, Table 1) garnet-opx Fe-Mg | opx, gt | |
| BK90olcpx | Brey & Kohler (1990, Table 1) olivine-cpx Fe-Mg | ol, cpx | |
| KB90 | Kohler & Brey (1990) Ca exchange between olivine and cpx, eqs. (5)-(6) | ol, cpx | |
| OW79 | O'Neill & Wood (1979) olivine-garnet Fe-Mg, eq. (26), with O'Neill (1980) | ol, gt | |
| WZ07 | Wu & Zhao (2007) olivine-garnet Fe-Mg, eq. (7) | ol, gt | |
| DH10Al | De Hoog et al. (2010) Al-in-olivine, eq. (12) | ol | `version=eq13` |
| DH10Cr | De Hoog et al. (2010) Cr-in-olivine, eq. (14) | ol | |
| BU17 | Bussweiler et al. (2017) Al-in-olivine, eq. (3) | ol | `version=eq2` |
| OW87 | O'Neill & Wall (1987) olivine-spinel Mg-Fe, eq. (20) | ol, sp | |
| BA91 | Ballhaus et al. (1991) olivine-spinel Mg-Fe | ol, sp | |
| LI95 | Li et al. (1995) olivine-spinel Mg-Fe, eqs. (9)-(10) | ol, sp | |
| WA08 | Wan et al. (2008) Al exchange between olivine and spinel | ol, sp | |
| CO14 | Coogan et al. (2014) Al exchange between olivine and spinel, eqs. (6)-(7) | ol, sp | |
| WES91 | Witt-Eickschen & Seck (1991) Al in opx with spinel, eq. (6) | opx, sp, ol | `version=Cr-Al-opx`: eq. (7), opx only |
| RY96Zn | Ryan et al. (1996) Zn-in-chromite, eq. (2) | sp | `coefficients=PTEXL` |
| CR09 | Creighton (2009) single-crystal Mn-in-garnet, eq. (10) | gt | `cr_term=printed` |
| GR89Ni | Griffin et al. (1989) Ni-in-garnet, eq. (IIa) | gt | `ni_ol=measured`: eq. (I), needs ol |
| RY96Ni | Ryan et al. (1996) Ni-in-garnet, eq. (1) | gt | |
| CA99Ni | Canil (1999) Ni-in-garnet, eq. (3) | gt, ol | `ni_ol=assumed`: 3000 ppm |
| SU21Ni | Sudholz et al. (2021) Ni-in-garnet, eq. (4) | gt | `ni_ol=measured` |
| NI24Ni | Nimis et al. (2024) Ni-in-garnet with P, eq. (2) | gt | `ni_ol=measured`; `version=T35`...`T60`: eqs. (3)-(7) |

**Barometers**

| Key | Calibration | Minerals | Options |
|---|---|---|---|
| NG85 | Nickel & Green (1985) Al-in-opx, eq. (5) | opx, gt | `mgts`, `mgts_in` (below) |
| T98 | NG85 with the Taylor (1998) Mg-Tschermak, i.e. `NG85,mgts=T98` | opx, gt | as NG85 |
| BKN | Brey & Kohler (1990) Al-in-opx, Table 3 | opx, gt | |
| BBG08 | Brey et al. (2008) Al-in-opx, eqs. (1)-(4), a refit of BKN up to 10 GPa | opx, gt | |
| MC74 | MacGregor (1974) Al<sub>2</sub>O<sub>3</sub> in opx with garnet in MAS, as fitted by Finnerty & Boyd (1984) | opx | |
| SUD21 | Sudholz et al. (2021) Cr-in-cpx | cpx | |
| NT00 | Nimis & Taylor (2000) Cr-in-cpx | cpx | `correction=N20` |
| AS17 | Ashchepkov et al. (2011) Jd-Di exchange in cpx for peridotites, as corrected in Ashchepkov et al. (2017, Appendix) | cpx | `kd=wt` |
| KB90 | Kohler & Brey (1990) Ca exchange between olivine and cpx, eqs. (5)-(6) | ol, cpx | |
| BE15 | Beyer et al. (2015) Al in cpx with garnet, for eclogite and garnet pyroxenite | cpx, gt | |
| GR06 | Grutter et al. (2006) Cr-in-garnet "P38", for chromite-saturated garnet | gt | `geotherm=35`...`48` |

Notes on thermobarometers:

- BKCa, NG09 and KR88 are used by the opx-garnet tests. Nimis & Grutter (2010) recommend NG09
  with NG85 for cpx-free garnet harzburgites and orthopyroxenites. They consider the
  garnet-cpx (KR88, EG79, PO85, AI94, KR00) and olivine-garnet (OW79) thermometers too
  imprecise to recommend, and Brey & Kohler (1990) do not recommend their Fe-Mg exchange
  thermometers (BK90gtopx, BK90olcpx). Several of the others have since been superseded (WE77, BM85, EG79,
  PO85, MC74), but are kept for comparison.
- Some of the single-mineral calibrations (NT00, SUD21, AS17, MC74, GR06, the olivine, chromite and garnet
  thermometers) assume that the mineral equilibrated with garnet, chromite or olivine
  (e.g. GR06 gives the true P only for chromite-saturated garnet, otherwise a minimum;
  the Ni-in-garnet thermometers need garnet that equilibrated with mantle olivine).
  Whether this criterion is met for a given sample is not checked.
- WE77, LI95, WES91, RY96Zn, CR09, the olivine-spinel Al thermometers (WA08, CO14) and the
  Ni-in-garnet thermometers (other than the default NI24Ni) do not depend on P.
- The GR06 barometer does not depend on T.
- Trace-element calibrations (Al and Cr in olivine, Zn in chromite, Ni in garnet) need the oxide
  measured to ppm-level precision (long-count EPMA, proton probe or LA-ICP-MS). For example, garnet NiO is often reported
  at or below 0.03 wt%, yielding Ni-in-garnet temperatures (e.g. RY96Ni) that
  are not meaningful. A zero (below detection) gives no result.
- Spinel Fe<sup>3+</sup> is estimated from stoichiometry (3 cations per 4 O; Droop, 1987), as done by Ballhaus et
  al. (1991), and is used by the olivine-spinel Mg-Fe and WES91 thermometers.
- D'Souza et al. (2020) combine the Al-in-olivine thermometer (CO14) with the Ca-in-olivine
  barometer (KB90) for spinel peridotites: `-c CO14:KB90`.
- KB90 as both thermometer and barometer uses the same reaction twice, which gives a meaningless result.
- BKCa `correction=NG10` refers to the Nimis & Grutter (2010) eq. (1) correction, applied at all T.
- NT00 barometer `correction=N20` is the empirical correction of Nimis et al. (2020, Appendix A.1),
  fitted to agree with opx-garnet barometry in xenoliths, which they call a band-aid measure.
- BM85: Bertrand & Mercier (1985) correct Ca-in-cpx for Fe with a term that depends on T, so T is
  solved for numerically. `opx_fe=PTEXL` adds the 0.075 Fe/(Fe+Mg) term to opx Ca that is used in
  PTEXL (the original paper neglects Fe in opx and the source of this term is not given). It raises T by a
  median of 2 °C on the AusXenos database.
- WE77: Wells uses the ideal two-site activity model of Wood & Banno (1973), eqs. (24)-(25), with
  half the Al placed on the M1 site. `sites=PTEXL` uses Nickel & Green (1985) charge balance instead, as implemented in PTEXL.
  It lowers T by a median of 6 °C on the AusXenos database.
- OW79: the volume term is that of the correction by O'Neill (1980), which reverses the sign of
  eq. (27) as printed in 1979 and fixes its brackets.
- BA91: printed as (6530 + 280P + 7000 + 108P)(1 - 2X<sub>Fe</sub><sup>ol</sup>), this has been read as
  6530 + 280P + (7000 + 108P)(1 - 2X<sub>Fe</sub><sup>ol</sup>), as in O'Neill & Wall (1987) eq. (20).
- WES91 eq. (6): X<sub>Mg</sub> on the opx M1 site is taken as (1 - Al<sup>M1</sup> - Cr - Ti) Mg/(Mg+Fe) because the paper does not
  define it. Originally calibrated on Eifel and Vogelsberg spinel peridotites, it becomes very sensitive to spinel Cr/(Cr+Al+Fe<sup>3+</sup>) at values above about 0.5.
- WZ07: garnet activities are from Holdaway (2001), with the polynomials of Wu & Zhao (2007)
  supplementary Table S2, checked for Gibbs-Duhem consistency in Ca-Fe-Mg garnet (N.B., Mn terms
  cannot be checked this way). Nimis & Grutter (2010) find it precise but that it gives systematically higher T than T98 below 1200 °C.
- BE15: a thermodynamic barometer for bimineralic garnet and cpx, calibrated on eclogite
  experiments at 3-7 GPa, which is useful for garnet pyroxenites without opx. It needs tetrahedral Al
  in cpx (Si < 2 per 6 O), so it gives no result for most peridotitic cpx. The printed garnet
  activity expressions, eqs. (14)-(15), have typos, so we use the standard asymmetric ternary form they cite from
  Ganguly et al. (1996). Cr, which is absent from their experiments, is placed on the cpx M1 site. It
  reproduces their experiments (Table 3) within 0.3 GPa. On the opx-cpx-garnet samples in the
  AusXenos database, T98:BE15 is a median 0.1 GPa lower than T98:NG85.
- AS17 barometer: an empirical Jd-Di exchange barometer, meant for any peridotitic cpx (garnet or
  spinel facies). The paper prints K<sub>D</sub> = Na/AlCr × MgO/CaO, with every other term in formula
  units. When interpreted as Mg/Ca in formula units (default), it agrees with NT00 to within ~0.3 GPa over the 2-8 GPa range,
  as the paper reports. Interpreting as the wt% oxide ratio (`kd=wt`) produces pressures 0.5-1.3 GPa lower.
  TiO<sub>2</sub> and Cr<sub>2</sub>O<sub>3</sub> are taken as required since leaving each blank moves P by a median 0.3 and 0.7 GPa, respectively, for AusXenos.
- AS17 thermometer: Ashchepkov et al. pair their barometer with NT00, corrected using cpx wt%
  of Al<sub>2</sub>O<sub>3</sub> + Cr<sub>2</sub>O<sub>3</sub>, FeO and Na<sub>2</sub>O. Their further step "for Cr low Cpx", (T + 390) P/(P - 70) 0.95,
  is negative at any mantle P as printed and is not applied.
- SU22gtcpx and SU22gtopx use X<sub>Ca</sub> = Ca/(Ca+Fe+Mg) in garnet, without Mn.
- Ni-in-garnet: Nimis et al. (2024, Table 1) computed CA99Ni and SU21Ni assuming 2900 ppm Ni-in-olivine,
  but the defaults here follow the original papers (3000 ppm where assumed).
- CR09: the Cr coefficient of ln a(spessartine) is printed as -31.5 in eq. (6). That puts ln a
  about 5 below ideal mixing for typical mantle garnets, so it was interpreted as a misprint for -3.15.
  `cr_term=printed` uses the original value of -31.5.
- RY96Zn: published coefficients are -0.90 and 0.26, but PTEXL uses -0.904 and 0.264. This
  results in a median T that is 31 °C cooler on the AusXenos database and is available as `coefficients=PTEXL`.
- MC74: MacGregor gives the Al<sub>2</sub>O<sub>3</sub> isopleths only in a figure. The equation we use is the fit
  in the appendix of Finnerty & Boyd (1984), eq. (A3), which reproduces MacGregor's garnet-field experiments
  (Table 4) within 3 kbar.
- BK90gtopx and BK90olcpx follow Brey & Kohler (1990) Table 1.

NG85: Our Nickel & Green (1985) barometer assumes Fe<sup>3+</sup> is zero and has a number of options:

- `mgts` sets the opx Mg-Tschermak mole fraction:
  - `NG85` (default): (Al - Cr - 2Ti + Na)/2, as in their original paper.
  - `CA91`: Carswell (1991) eqs. (E6) and (E9). When Na > Cr + Ti, the excess Na forms jadeite, so (Al + Cr - Na)/2.
  - `CG87`: Carswell & Gibb (1987b), eq. (E7) in Carswell (1991). If Na > Cr + 2Ti, then (Al + Cr + 2Ti - Na)/2.
    Carswell (1991) says this version is only valid if Ti sits in a NaTi<sub>0.5</sub>Mg<sub>0.5</sub>Si<sub>2</sub>O<sub>6</sub> component.
  - `T98`: Taylor (1998) uses (Al - Cr + Ti + Na)/2, which is an empirical Ti correction.
- `mgts_in` sets which NG85 terms use the modified Mg-Tschermak value. NG85 use X<sub>Al</sub><sup>M1</sup> in the equilibrium
  constant, in ΔV, and in the M1 Mg+Fe occupancy. The later papers do not say which of these to replace.
  - `activity_dV` (default): updates X<sub>Al</sub><sup>M1</sup> in the equilibrium constant and ΔV, while the M1 Mg+Fe occupancy keeps the stoichiometric NG85 value.
  - `activity`: updates the equilibrium constant only.
  - `all`: updates X<sub>Al</sub><sup>M1</sup> in all three.

We compared these options to the original NG85 for the 97 samples with a T98:NG85 result in the 2026 version of the
AusXenos database. Each option only differs from NG85 for some compositions, so it changes P in only some of the samples:

- `CA91` differs when opx Na > Cr + Ti, which is the case in 32 samples. It raises P in those by a median of 0.05 GPa, up to 0.25 GPa.
- `CG87` differs when Na > Cr + 2Ti (14 samples), raising P by a median of 0.05 GPa, up to 0.21 GPa.
- `T98` differs whenever opx contains Ti (the 78 samples with opx TiO<sub>2</sub> reported), lowering P by a median of 0.11 GPa, up to 0.55 GPa.
- The three `mgts_in` choices change P by at most 0.05 GPa.

## Required oxides

A blank cell may mean an oxide was not analysed or that it was below the limit of detection, so blanks are
treated according to how heavily each thermobarometer depends on that oxide:

- **Required**: SiO<sub>2</sub>, Al<sub>2</sub>O<sub>3</sub>, FeO<sub>t</sub>, MgO and CaO for a sound formula unit, plus minor oxides whose
  blank would typically move the result by more than about half of the calibration's stated uncertainty
  (e.g. cpx Na<sub>2</sub>O for the two-pyroxene thermometers, cpx Cr<sub>2</sub>O<sub>3</sub> and Na<sub>2</sub>O for NT00 and SUD21, olivine
  CaO for KB90, garnet Cr<sub>2</sub>O<sub>3</sub> and MnO for CR09, and the trace element itself for the olivine and
  chromite thermometers). Some calibrations also need olivine Cr<sub>2</sub>O<sub>3</sub>: DH10Al, and BU17 with `version=eq2`.
  If any of these is blank, the thermobarometer pair gives no result, with `status = missing oxides` and the
  blanks listed in `missing_oxides`.
- **Minor**: the other oxides the thermobarometer uses. A blank is taken as zero and listed under `assumed zero` in the warnings.

A sample with blank cpx Cr<sub>2</sub>O<sub>3</sub> therefore yields no Cr-in-cpx pressure, but only raises a warning with T98:NG85.
`python -m xenopt --list` shows the lists for each calibration. `--blanks-as-zero` turns this
check off and interprets every blank as zero.

To decide whether an oxide is required or minor, we set each oxide to zero in turn for samples in the 2026 AusXenos
database, recomputed P-T, and assessed the median (90th percentile) changes. Using T98:NG85 as an example, 97 samples yield a P-T
result and, for each oxide in only those samples where it was reported, we obtain:

| Oxide | Median (90th percentile) change | Conclusion |
|---|---|---|
| SiO<sub>2</sub>, Al<sub>2</sub>O<sub>3</sub>, FeO<sub>t</sub>, MgO, CaO in opx, cpx and garnet | not tested | Required, for a sound formula unit |
| cpx Na<sub>2</sub>O | 236 °C (322 °C), 0.88 (1.41) GPa | Required |
| cpx Cr<sub>2</sub>O<sub>3</sub> | 5 °C (37 °C), 0.02 (0.18) GPa | Minor |
| cpx TiO<sub>2</sub> | 6 °C (31 °C), 0.02 (0.11) GPa | Minor |
| opx Cr<sub>2</sub>O<sub>3</sub> | 3 °C (11 °C), 0.04 (0.28) GPa | Minor |
| opx TiO<sub>2</sub> | 3 °C (8 °C), 0.08 (0.17) GPa | Minor |
| opx Na<sub>2</sub>O | 1 °C (4 °C), 0.02 (0.09) GPa | Minor |
| opx MnO | 0 °C (1 °C), 0.01 (0.02) GPa | Minor |
| garnet Cr<sub>2</sub>O<sub>3</sub> | 0 °C (2 °C), 0.01 (0.06) GPa | Minor |
| garnet MnO | 0 °C (0 °C), 0.00 (0.00) GPa | Minor |

A blank minor oxide can still shift P by up to 0.5 GPa in a handful of samples (opx Cr<sub>2</sub>O<sub>3</sub>, TiO<sub>2</sub> or Na<sub>2</sub>O),
which is why they still appear as warnings.

## Clinopyroxene suitability tests

Single-clinopyroxene thermobarometry has largely been developed for use in heavy mineral concentrates where,
unlike with xenoliths, it is no longer possible to visually confirm coexistance with other minerals like orthopyroxene and garnet. 
The single-cpx NT00 thermometer and SUD21 barometer pairing, neverthless, still requires that the cpx grew in a
garnet peridotite in equilibrium with both opx and garnet. Enstatite-in-cpx (opx saturation) sets T,
while Cr exchange with garnet sets P. There have therefore been a series of tests developed to
ascertain whether a clinopyroxene grain was likely to be in equilibrium with opx and/or garnet based on its
geochemistry alone. We compute many of the proposed protocols and each column is given 1/0 if the analysis passes/fails.

**Analysis quality.** Poor analyses give wrong cation proportions and the thermobarometers depend on
small differences between minor elements.

- `cpx_oxide_total`: oxide total in the range 98-102 wt% (e.g., Sudholz et al., 2021). Totals outside this range
  point to either a poor analysis (e.g. surface, beam, or standardisation problems) or excessive alteration/weathering.
- `cpx_cation_total`: cations per 6 O in the range 3.98-4.02 (e.g., Nimis & Taylor, 2000; Ziberna et al., 2016).
  An ideal pyroxene should have 4. Deviations mean either a poor analysis or significant Fe<sup>3+</sup>, which the
  calculation (all Fe as Fe<sup>2+</sup>) cannot handle. Grutter (2009) allowed 3.96-4.04; Ziberna et al.
  accept wider limits in some cases but warn they may increase scatter around geotherms.
  Hoggard et al. (2026) used 3.97-4.03 (`--cpx-cation-range 3.97 4.03`), as many analyses near
  4.025 passed all other tests and yielded P-T conditions that agreed with those that did pass.

**Origin in garnet peridotite.** Single-cpx barometers generally require cpx to have equilibrated with garnet.

- `nimis98_cr_al`: requires the cpx to be inside the garnet-peridotite field of Ramsay's (1995) Cr<sub>2</sub>O<sub>3</sub> vs Al<sub>2</sub>O<sub>3</sub> diagram:
  Cr<sub>2</sub>O<sub>3</sub> > 0.5 wt%, and Al<sub>2</sub>O<sub>3</sub> <= 4.0 wt% if Cr<sub>2</sub>O<sub>3</sub> < 2.25 wt%, or <= 5.0 wt% otherwise
  (Nimis, 1998; Ziberna et al., 2016, step 3). It separates garnet-peridotite from cpxs from eclogites,
  pyroxenites, spinel peridotites, megacrysts and cognate phenocrysts, which are typically more
  aluminous or Cr-poor. It may reject some genuine garnet-peridotite cpx.
- `nimis98_al_mg`: requires cpx to be in the high-Al field of Al<sub>2</sub>O<sub>3</sub> vs MgO: Al<sub>2</sub>O<sub>3</sub> >= 0.7 wt% and
  Al<sub>2</sub>O<sub>3</sub> >= 12.175 - 0.6375 MgO (Nimis, 1998, Fig. 3; Ziberna et al., 2016, step 4). Some
  diopsides from low-P, refractory, garnet-free peridotites affected by (often carbonatitic)
  metasomatism plot in Ramsay (1995)'s garnet field but have anomalously low Al<sub>2</sub>O<sub>3</sub> for their MgO. Al is
  the main pressure indicator, so these would give misleading results. Again it may reject some
  genuine garnet peridotite cpx, especially those from the diamond stability field.
- `ptexl_na_k`: Na + K >= Al + Cr - 0.08, i.e. a low Tschermak content (Al + Cr - Na - K <=
  0.08). Grutter (2009) showed that cpx with high Al + Cr - Na generally come from spinel or spinel +
  garnet assemblages, where the Cr-in-cpx barometer underestimates P. The 0.08 cut-off is from
  PTEXL.

**Within the calibrated, well-conditioned range of the barometer.**

- `sudholz21_na`: requires Na <= 0.2 cations (Sudholz et al., 2021). In natural xenoliths, cpx with
  > 0.2 Na gave Cr-in-cpx pressures typically more than 1 GPa below opx-garnet (BKN) pressures.
  Highly sodic cpx are also typical of eclogites, and P estimates on eclogitic cpx should be
  discarded.
- `sudholz21_cr_no`: 0.10 <= Cr# <= 0.65 (Ziberna et al., 2016, step 5; adopted by Sudholz et
  al., 2021). This is the composition range covered by experiments. Ziberna et al. advise caution for
  0.50-0.65, where testing is limited.
- `sudholz21_cacrts`: a<sub>CaCrTs</sub> / Cr# >= 0.005 (Sudholz et al., 2021). When this ratio is low,
  small analytical errors in Cr, Na and Al give large P errors, and the composition lies outside
  the experiments. Ziberna et al. (2016, step 6) discarded values <= 0.011, the lower end of the
  NT00 experiments; Sudholz et al.'s new experiments extend the calibration down to 0.005.
- `ziberna16_cr2o3`: 0 < Cr<sub>2</sub>O<sub>3</sub> <= 5 wt% (Ziberna et al., 2016, step 3). Above 5 wt% the
  thermobarometer is untested, and with no Cr it cannot be applied at all.

**Equilibrium with opx.** The single-cpx enstatite-in-cpx thermometer assumes opx was present. Without it,
T is underestimated, and P with it, since the barometer depends on T.

- `ziberna16_ca_mg`: Ca / (Ca + Mg) <= 0.5 (Ziberna et al., 2016, step 7). Only about 1% of
  opx-saturated chromian diopsides lie above 0.5, so higher values suggest no coexisting opx.
- `ziberna16_T_<pair>`: T >= 600 °C (Ziberna et al., 2016, step 7). A very low calculated T is a
  strong sign the cpx was not in equilibrium with opx. Nimis & Grutter (2010) suggest 700 °C as a
  more conservative cut-off (`--T-min 700`). Since the test relies on computed T, it is evaluated for every therombarometer pairing.

`cpx_present` only records whether a cpx analysis exists. Ziberna et al. (2016) also recommend
discarding zoned grains and checking P-T plots for outliers, neither of which can be tested from a single clinopyroxene
analysis.

These tests are listed in the warnings of pairs that use cpx only (e.g. `NT00:SUD21`).

## Orthopyroxene-garnet tests

For thermobarometer pairings that use opx, cpx and garnet (e.g. `T98:NG85`), the checks of Nimis & Grutter (2010) are
applied. They compiled a reference set of well-equilibrated garnet peridotite and pyroxenite
xenoliths, and recommend the T98 thermometer with the NG85 barometer plus these checks.
The temperature used in these tests are always from the T98 thermometer, calculated at the pressure found by the chosen thermobarometer pair,
because Nimis & Grütter (2010)'s thresholds were determined using T98. This means that, for T98:NG85, this is the same temperature.
For other pairs (e.g. BKN:BKN), however, it will use a T that may differ from chosen thermometers T.

**Analysis quality**

- `ng10_opx_quality`, `ng10_cpx_quality`: orthopyroxene oxide total in the range 98.5-101.5 wt% and at least 3.99 cations
  per 6 O. Low cation sums mean a poor analysis, which matters most for the minor elements the
  pyroxene thermometers and the Al-in-opx barometer depend on.
- `ng10_gt_quality`: garnet oxide total in the range 98.5-101.5 wt%. Nimis & Grütter (2010) saw no need for a garnet cation-sum test,
  since garnet composition has relatively little influence on the results.

**Rock type.**

- `ng10_gt_type`: garnet Cr<sub>2</sub>O<sub>3</sub> > 1 wt%, or Mg/(Mg+Fe) > 0.6 (after Grutter et al., 2004). The test
  is intended to isolate garnets that come from peridotites and mantle pyroxenites and exclude eclogitic or crustal garnets.

**Plausible result.**

- `ng10_T_min`: T98 T >= 700 °C. Pyroxene thermometry is unreliable at very low T and Nimis & Grütter (2010) discarded such samples.
- `ng10_P_min`: P should be >= 1.5 GPa, which is the minimum pressure at which garnet is stable in lherzolite.
- `ng10_P_min_pyroxenite`: P shoule be >= 1.0 GPa for pyroxenites, which is the minimum pressire at which garnet is stable in olivine-free rocks.

  At lower P, the result is either likely to be wrong or the rock is not a garnet lherzolite or pyroxenite. Since lithology
  is often not recorded, both are tested for every sample. Between 1.0 and 1.5 GPa, the warning
  reads `low P (ok if pyroxenite)`, while below 1.0 GPa both tests fail and it reads `very low P`.

**Equilibrium between the minerals.** The thermobarometer pair assumes opx, cpx and garnet equilibrated
together. One way to test this is to compare independent thermometers that work on the separate minerals.

- `ng10_cpx_opx`: T98 and the Ca-in-opx thermometer (BKCa with Nimis & Grutter's 2010 correction)
  should agree within 90 °C when T98 < 900 °C, 70 °C when T98 is in the range 900-1200 °C, or 50 °C when T98 > 1200 °C. These values are about
  twice the standard error within each range. T98 depends mainly on the cpx composition and Ca-in-opx only on
  the opx, so a larger difference points to pyroxene disequilibrium or poor analyses. Note that PTEXL
  widens these limits by 10 °C (`--ng10-cpx-opx-tol 100 80 60`).
- `ng10_gt_opx`: T98 and the opx-garnet Fe-Mg thermometer (NG09) agree within 70 °C, about twice
  its standard error. T98 is based on Ca-Mg exchange between pyroxenes and NG09 on Fe-Mg
  exchange with garnet. A larger difference may mean these equilibria were decoupled,
  unusual redox conditions (all Fe is treated as Fe<sup>2+</sup>), or poor analyses.

**Accuracy limits.** These tests raise cautions rather than failures.

- `ng10_T_accuracy`: temperature from T98 should be >= 900 °C, as at cooler temperatures, the accuracy cannot be
  assured from available experiments.
- `ng10_P_accuracy`: P should be <= 6 GPa as NG85 may underestimate P above this pressure.

**Optional garnet-cpx check.** The T98 temperature can also be compared with the garnet-cpx Fe-Mg thermometer of
Krogh (1988), KR88. PTEXL flags a difference above 90 °C. However, Nimis & Grutter (2010) found all garnet-cpx
thermometers too imprecise for this comparison (differences > 200 °C in well-equilibrated xenoliths,
probably due to Fe<sup>3+</sup> in cpx) and do not recommend them, so this check is switched off by default.
`--gt-cpx-tol 90` adds it as `ptexl_gt_cpx`.

## Adding a thermometer or barometer

It is possible to add your own thermometer or barometer in xenoPT. To do so, 
write a function in `xenopt/thermobarometers.py` that takes the cations mapping and P [GPa]
(thermometer, returns T [K]) or T [K] (barometer, returns P [GPa]). Use column names, e.g.
`c["cpx"].Na`. End its docstring with a `Source:` paragraph giving the full reference and the
equation, table or page it implements (a test checks every registered calibration has one). Then
add it to `THERMOMETERS` or `BAROMETERS` with its required and minor oxides
for each mineral. Variants
are keyword arguments listed in the entry's `options`, with the published version first
(`option_oxides` adds required oxides for an option). Calibrations that need oxide wt% rather
than cations (e.g. trace-element thermometers in ppm) set `needs_wt=True` and take a `wt`
argument; `chemistry.element_ppm` converts. It can then be paired with any existing calibration
from the command line.

Run `pip install -e ".[test]" && pytest` to check the addition hasn't caused changes to other functionality. These tests also compare each
calibration with values from its paper.

## Reproducing Hoggard et al. (2026)

To analyse the 2026 AusXenos database with the settings of Hoggard et al. (2026), reading it directly from GitHub, run:

```bash
url=https://raw.githubusercontent.com/mjh217/ausxenos/v2026.11/csv_version
python -m xenopt $url/AusXenos_geochem.csv --localities $url/AusXenos_locs.csv \
    -c T98:NG85 -c NT00:SUD21 --cpx-cation-range 3.97 4.03 --blanks-as-zero -o results.csv
```

They rejected all cpx-only P-T results failing `cpx_oxide_total`, `cpx_cation_total`, `nimis98_al_mg`, `ziberna16_cr2o3`,
`ziberna16_T`, `sudholz21_na` or `sudholz21_cacrts`. They used `nimis98_cr_al`, `sudholz21_cr_no` and
`ziberna16_ca_mg` only as guides, did not use `ptexl_na_k`, and removed outliers by inspecting the P-T arrays of each site.
For opx-cpx-garnet xenoliths (`T98:NG85`), they checked only the oxide totals.

## Citation

If you use xenoPT, please cite:

> Hoggard, M. J., Sudholz, Z. J., Hawkins, R., Duvernay, T., Bagchi, S., Richards, F. D., Jaques, A. L., Yaxley, G. M., Czarnota, K., & Haynes, M. (2026). Paleothermal structure of the Australian continental lithosphere. *Geochemistry, Geophysics, Geosystems* (accepted).

Please also cite the original papers of the thermobarometers you use.

Each release of xenoPT is archived on Zenodo at [10.5281/zenodo.24386935](https://doi.org/10.5281/zenodo.24386935), which always resolves to the latest version; each release also has its own DOI (v2026.11: [10.5281/zenodo.24386936](https://doi.org/10.5281/zenodo.24386936)).

## Licence

xenoPT is released under the MIT licence. See `LICENSE`.
