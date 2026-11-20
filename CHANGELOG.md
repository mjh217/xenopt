# Changelog

Releases use dated version numbers (`vYYYY.MM`). Each GitHub release is archived on Zenodo with its own DOI.

## v2026.11

Initial public release, accompanying Hoggard et al. (2026, *Geochemistry, Geophysics, Geosystems*).

- 38 thermometers and 11 barometers, each checked against its original publication, with published variants as options.
- Iterative P-T solution for any thermometer-barometer pair.
- Required and minor oxides for each calibration: blank required oxides give no result, blank minor oxides are taken as zero with a warning.
- Clinopyroxene suitability tests (Ziberna et al., 2016; Sudholz et al., 2021) and orthopyroxene-garnet tests (Nimis & Grutter, 2010), summarised as plain-language warnings.
- Reads the AusXenos database format directly.
