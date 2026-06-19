*! lwdid_reference.do
*! Generate reference results from Stata lwdid v2.4.2 for validating the R port (lwdidR).
*!
*! HOW TO RUN (once, in Stata):
*!   1. Make sure lwdid v2.4.2 is installed, e.g.:
*!        net install lwdid, from("https://raw.githubusercontent.com/Soo-econ/lwdid/main/") replace
*!      (or:  ssc install lwdid, replace)
*!   2. In Stata:  cd to THIS folder, then:  do lwdid_reference.do
*!
*! OUTPUT (written to this folder):
*!   ref_walmart_<config>.csv   -- one per large-N config: WATT(r) + Pre/Post avg table
*!   ref_smoking_smallN.csv     -- small-N single (common timing) scalars + period effects
*!   ref_attgt_demean_ipwra.log -- a logged ATT(g,t) run for spot-checking cells

clear all
set more off
version 16.0

* fail early if lwdid is not installed
cap which lwdid
if _rc {
    di as error "lwdid is not installed. Install it first (see header of this do-file)."
    exit 111
}
di as result "lwdid version in use:"
which lwdid

* Export the datasets to CSV so the R tests can load the EXACT same data
* without a Stata-file reader dependency.
use "lw_walmart.dta", clear
export delimited using "lw_walmart.csv", replace
use "lw_smoking.dta", clear
export delimited using "lw_smoking.csv", replace

***********************************************************************
* PART A. LARGE-N PATH (lw_walmart): the new functionality to validate
*   y = log_wholesale_emp ; covars = x1 x2 x3 ; ivar=cid tvar=year gvar=first_year
*   save() writes a .dta with columns:
*     effect ryear watt se t_stat p_value low_ci up_ci N_cells N_units
*   Re-load each saved file and export to CSV (one CSV per config).
***********************************************************************

* helper: first token = config label; the REST = option string (gettoken keeps
* multi-token option strings intact, which `args` would not).
capture program drop _ln_run
program define _ln_run
    gettoken label 0 : 0
    use "lw_walmart.dta", clear
    lwdid log_wholesale_emp x1 x2 x3, ivar(cid) tvar(year) gvar(first_year) ///
          `0' save(_tmp_res)
    use "_tmp_res.dta", clear
    gen str40 config = "`label'"
    order config
    export delimited using "ref_walmart_`label'.csv", replace
end

* helper for the no-covariate path (outcome only)
capture program drop _ln_run_nox
program define _ln_run_nox
    gettoken label 0 : 0
    use "lw_walmart.dta", clear
    lwdid log_wholesale_emp, ivar(cid) tvar(year) gvar(first_year) ///
          `0' save(_tmp_res)
    use "_tmp_res.dta", clear
    gen str40 config = "`label'"
    order config
    export delimited using "ref_walmart_`label'.csv", replace
end

* rolling x method grid (with covariates x1 x2 x3)
_ln_run demean_ra     rolling(demean)  method(ra)
_ln_run detrend_ra    rolling(detrend) method(ra)
_ln_run demean_ipw    rolling(demean)  method(ipw)
_ln_run detrend_ipw   rolling(detrend) method(ipw)
_ln_run demean_ipwra  rolling(demean)  method(ipwra)
_ln_run detrend_ipwra rolling(detrend) method(ipwra)

* option variants
_ln_run detrend_ipwra_pre3 rolling(detrend) method(ipwra) pre(3)
_ln_run demean_ipwra_never rolling(demean)  method(ipwra) never

* no-covariate RA path
_ln_run_nox demean_ra_nox rolling(demean) method(ra)

* logged ATT(g,t) cells for one config (spot-check per-cell estimates)
use "lw_walmart.dta", clear
log using "ref_attgt_demean_ipwra.log", text replace
lwdid log_wholesale_emp x1 x2 x3, ivar(cid) tvar(year) gvar(first_year) ///
      rolling(demean) method(ipwra) attgt
log close

***********************************************************************
* PART B. SMALL-N SINGLE (lw_smoking, common timing): regression-test
*   the existing lwdidR small-N path against v2.4.2.
*   y=lcigsale ivar=state tvar=year gvar=first_year
***********************************************************************

tempname pf
postfile `pf' str24 config str24 key double value using "_ref_smoking.dta", replace

* --- demean ---
use "lw_smoking.dta", clear
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(demean)
post `pf' ("demean") ("att")    (e(att))
post `pf' ("demean") ("se_att") (e(se_att))
cap matrix A = e(ATTt)
cap matrix S = e(SEt)
if _rc == 0 {
    forvalues j = 1/`=colsof(A)' {
        post `pf' ("demean") ("ATTt_`j'") (A[1,`j'])
        post `pf' ("demean") ("SEt_`j'")  (S[1,`j'])
    }
}

* --- detrend ---
use "lw_smoking.dta", clear
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend)
post `pf' ("detrend") ("att")    (e(att))
post `pf' ("detrend") ("se_att") (e(se_att))
cap matrix A = e(ATTt)
cap matrix S = e(SEt)
if _rc == 0 {
    forvalues j = 1/`=colsof(A)' {
        post `pf' ("detrend") ("ATTt_`j'") (A[1,`j'])
        post `pf' ("detrend") ("SEt_`j'")  (S[1,`j'])
    }
}

* --- detrend, vce(hc3) ---
use "lw_smoking.dta", clear
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) vce(hc3)
post `pf' ("detrend_hc3") ("att")    (e(att))
post `pf' ("detrend_hc3") ("se_att") (e(se_att))

* --- detrend, randomization inference (fixed seed) ---
use "lw_smoking.dta", clear
lwdid lcigsale, small ivar(state) tvar(year) gvar(first_year) rolling(detrend) ///
      ri riseed(349139) rireps(999)
post `pf' ("detrend_ri") ("att")    (e(att))
post `pf' ("detrend_ri") ("se_att") (e(se_att))
post `pf' ("detrend_ri") ("p_ri")   (e(p_ri))

postclose `pf'
use "_ref_smoking.dta", clear
export delimited using "ref_smoking_smallN.csv", replace

* tidy temp files
cap erase "_tmp_res.dta"
cap erase "_ref_smoking.dta"

di as result "DONE. Reference CSVs written to this folder."
