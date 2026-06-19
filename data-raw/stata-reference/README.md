# Stata reference oracle

`lwdid_reference.do` regenerates the reference results used to validate `lwdidR`'s
large-N path against the canonical Stata command, `lwdid`
(<https://github.com/Soo-econ/lwdid>).

## How to regenerate

1. In Stata, install `lwdid` and fetch its example datasets:

   ```stata
   net install lwdid, from("https://raw.githubusercontent.com/Soo-econ/lwdid/main/") replace
   net get lwdid           // downloads lw_smoking.dta and lw_walmart.dta
   ```

   Copy `lw_smoking.dta` and `lw_walmart.dta` into this folder. (They are not
   committed here: `lw_walmart` already ships compressed at
   `inst/extdata/lw_walmart.rds`, and both belong to the original `lwdid` package.)

2. Run the do-file from this folder:

   ```stata
   cd ".../lwdidR/data-raw/stata-reference"
   do lwdid_reference.do
   ```

This writes (all git-ignored): `lw_walmart.csv`, `lw_smoking.csv`,
`ref_walmart_<config>.csv` (9 large-N configs), `ref_smoking_smallN.csv`, and a
logged `ref_attgt_demean_ipwra.log`.

## How the R tests use them

`tests/testthat/test-large_n.R` reads the `ref_*.csv` files if present and asserts
that `lwdidR` reproduces Stata's WATT(r)/Pre/Post point estimates to ~1e-5 and the
standard errors to within wild-bootstrap Monte-Carlo error. If the CSVs are absent
(e.g., on CI without Stata) those tests skip cleanly.

Validation status (Stata `lwdid` v2.4 vs `lwdidR` 0.2.0): all 9 large-N configs
match point estimates to <=1.5e-8 and SEs to within ~3%.
