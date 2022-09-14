set more off
set varabbrev off
set linesize 132


if `"`c(os)'"' == "MacOSX" & `"`c(username)'"' == "nix" cap cd "~/Documents/Github/RA_FALL22" 
//Collaborator 1: MAC OS
if `"`c(os)'"' == "Windows" & `"`c(username)'"' == "16083" cap cd "~/Documents/Github/RA_FALL22" 

if `"`c(os)'"' == "Windows" & `"`c(username)'"' == "natorre" cap cd "C:\Users\natorre@emory.edu\Documents\GitHub\RA_Fall22" 
* do Code/Replication_paper_SM
* do Code/Replication_paper_SM

use Data/analy_malp_paper_9_22.dta, clear

*global control "pctpoverty unemp rpcinc"
*global torts "jsl pcap necap apology"
*global atorts "pcap necap apology"

  
/*Main models*/

*poisson  malp_np   fpa $torts $control i.stfips i.year if Mstclass<3, exp(pop) cluster(stfips)
*poisson  malp_phys fpa $torts $control i.stfips i.year if Mstclass<3, exp(pop) cluster(stfips)

*foreach dvar in adv_np vunsafe_np vRX_np vSOP_np  {
	
*poisson  `dvar' fpa $atorts $control i.stfips i.year if Astclass<3, exp(pop) cluster(stfips)


*}


*foreach dvar in adv_phys vunsafe_phys vRX_phys vSOP_phys {
	
*poisson  `dvar' fpa $atorts $control i.stfips i.year if Astclass<3, exp(pop) cluster(stfips)


*}

/****** Adapting Woldrige Code ************/

/* Load the data */
use Data/analy_malp_paper_9_22.dta, clear

/* Declare a panel data */
xtset stfips year

/* Define locals controls, torts, atorts */
local control pctpoverty unemp rpcinc
local torts jsl pcap necap apology
local atorts pcap necap apology


/****** Woldrige Variables ************/

/* d=1 if recieve the treatment. */
gen d = (!missing(FPA_FULL_YEAR))
label variable d "=1 if treated"

** p=1 if post-treatment time
gen p = (pre==0)
label variable p "=1 if post-treatment period"

** w=d*p
gen w= (d*p)
label variable w "=1 if treated in post-treatment period"

** year dummies
tab year, gen(f)
rename f# f#, renumber(1998)

* varlist of year dummies
unab i_year : f1* f2*
di "`i_year'"
** cohort dummies
tab FPA_FULL_YEAR, gen(c)

** replacing missing values
foreach v of varlist c* {
    replace `v' = 0 if `v' == .
}

* varlist of cohort dummies
unab cohort : c*


** Generate the demeaned variables conditioned by cohort

foreach z in `cohort' {

	local control_dm_`z'
	local x "control_dm_`z'"

		foreach v in `control' {

			sum `v' if `z'
			gen `v'_dm_`z' = `v' - r(mean)
			label variable `v'_dm_`z' "`v' demeaned at cohort `z'"
			local `x' `control_dm_`z'' `v'_dm_`z'
		}
}

unab prueba : `control_dm_c1'-`control_dm_c14'
di "`prueba'"

/* poisson */

* Now with a covariate. The results are not the same, but the differences
* in the estimated ATTs are minor in this application.
* Also, the results change with full time dummies and full  
* of time dummy interactions with covariates.

xtpoisson malp_np w#c.c1#c.f1998-f2019 w#c.c2#c.f2000-f2019 w#c.c3#c.f2005-f2019 ///
				w#c.c4#c.f2006-f2019 w#c.c5#c.f2010-f2019 w#c.c6#c.f2011-f2019 ///
				w#c.c7#c.f2012-f2019 w#c.c8#c.f2013-f2019 w#c.c9#c.f2014-f2019 ///
				w#c.c10#c.f2015-f2019 w#c.c11#c.f2016-f2019 w#c.c12#c.f2018-f2019 w#c.c13#c.f2019 ///
				w#c.c1#c.f1998-f2019#c.`control_dm_c1' w#c.c2#c.f2000-f2019#c.`control_dm_c2' ///
				w#c.c3#c.f2005-f2019#c.`control_dm_c3' w#c.c4#c.f2006-f2019#c.`control_dm_c4' ///
				w#c.c5#c.f2010-f2019#c.`control_dm_c5' w#c.c6#c.f2011-f2019#c.`control_dm_c6' /// 
				w#c.c7#c.f2012-f2019#c.`control_dm_c7' w#c.c8#c.f2013-f2019#c.`control_dm_c8' ///
				w#c.c9#c.f2014-f2019#c.`control_dm_c9' w#c.c10#c.f2015-f2019#c.`control_dm_c10' ///
				w#c.c11#c.f2016-f2019#c.`control_dm_c11' w#c.c12#c.f2018-f2019#c.`control_dm_c12' ///
				w#c.c13#c.f2019#c.`control_dm_c13' ///
				i.year i.year#c.`control', fe
				
poisson malp_np w#c.c1#c.f1998-f2019 w#c.c2#c.f2000-f2019 w#c.c3#c.f2005-f2019 ///
				w#c.c4#c.f2006-f2019 w#c.c5#c.f2010-f2019 w#c.c6#c.f2011-f2019 ///
				w#c.c7#c.f2012-f2019 w#c.c8#c.f2013-f2019 w#c.c9#c.f2014-f2019 ///
				w#c.c10#c.f2015-f2019 w#c.c11#c.f2016-f2019 w#c.c12#c.f2018-f2019 w#c.c13#c.f2019 ///
				///
				w#c.c1#c.f1998-f2019#c.`control_dm_c1' w#c.c2#c.f2000-f2019#c.`control_dm_c2' ///
				w#c.c3#c.f2005-f2019#c.`control_dm_c3' w#c.c4#c.f2006-f2019#c.`control_dm_c4' /// 
				w#c.c5#c.f2010-f2019#c.`control_dm_c5' w#c.c6#c.f2011-f2019#c.`control_dm_c6' /// 
				w#c.c7#c.f2012-f2019#c.`control_dm_c7' w#c.c8#c.f2013-f2019#c.`control_dm_c8' ///
				w#c.c9#c.f2014-f2019#c.`control_dm_c9' w#c.c10#c.f2015-f2019#c.`control_dm_c10' ///
				w#c.c11#c.f2016-f2019#c.`control_dm_c11' w#c.c12#c.f2018-f2019#c.`control_dm_c12'/// 
				w#c.c13#c.f2019#c.`control_dm_c13' ///
				i.year i.year#c.`control'///
				`cohort' `control' c.c*#c.`control', vce(cluster stfips)			
						
				
poisson malp_np w#c.c1#c.f1998-f2019 w#c.c2#c.f2000-f2019 w#c.c3#c.f2005-f2019 ///
				w#c.c4#c.f2006-f2019 w#c.c5#c.f2010-f2019 w#c.c6#c.f2011-f2019 ///
				w#c.c7#c.f2012-f2019 w#c.c8#c.f2013-f2019 w#c.c9#c.f2014-f2019 ///
				w#c.c10#c.f2015-f2019 w#c.c11#c.f2016-f2019 w#c.c12#c.f2018-f2019 w#c.c13#c.f2019 ///
				///
				w#c.c1#c.f1998-f2019#c.`control_dm_c1' w#c.c2#c.f2000-f2019#c.`control_dm_c2' ///
				w#c.c3#c.f2005-f2019#c.`control_dm_c3' w#c.c4#c.f2006-f2019#c.`control_dm_c4' /// 
				w#c.c5#c.f2010-f2019#c.`control_dm_c5' w#c.c6#c.f2011-f2019#c.`control_dm_c6' /// 
				w#c.c7#c.f2012-f2019#c.`control_dm_c7' w#c.c8#c.f2013-f2019#c.`control_dm_c8' ///
				w#c.c9#c.f2014-f2019#c.`control_dm_c9' w#c.c10#c.f2015-f2019#c.`control_dm_c10' ///
				w#c.c11#c.f2016-f2019#c.`control_dm_c11' w#c.c12#c.f2018-f2019#c.`control_dm_c12'/// 
				w#c.c13#c.f2019#c.`control_dm_c13' ///
				i.year i.year#c.`control'///
				`cohort' `control' c.c*#c.`control' ///
				c.c*#c.t c.c*#c.t#c.`control', vce(cluster stfips)			
				
*To do : Fix Margins 
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f02 = 0 f03 = 0 f04 = 1 f05 = 0 f06 = 0) ///
	subpop(if d4 == 1) noestimcheck vce(uncond)
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f02 = 0 f03 = 0 f04 = 0 f05 = 1 f06 = 0) ///
	subpop(if d4 == 1) noestimcheck vce(uncond)  
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f02 = 0 f03 = 0 f04 = 0 f05 = 0 f06 = 1) ///
	subpop(if d4 == 1) noestimcheck vce(uncond) 
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f02 = 0 f03 = 0 f04 = 0 f05 = 1 f06 = 0) ///
	subpop(if d5 == 1) noestimcheck vce(uncond)
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f02 = 0 f03 = 0 f04 = 0 f05 = 0 f06 = 1) ///
	subpop(if d5 == 1) noestimcheck vce(uncond) 
margins, dydx(w) at(d4 = 0 d5 = 0 d6 = 1 f02 = 0 f03 = 0 f04 = 0 f05 = 0 f06 = 1) ///
	subpop(if d6 == 1) noestimcheck vce(uncond)



	
	
* Callaway and Sant'Anna (2020, Journal of Econometrics)
gen first_treat = FPA_FULL_YEAR
replace first_treat = 0 if first_treat == .

csdid malp_np pctpoverty unemp rpcinc , ivar(stfips) time(year) gvar(first_treat)














