*****************************************************************
*	Replication Exercise - staggered entry case					*
*																*
*	Implementing Two-Way Fixed Effects,							*
*	the Two-Way Mundlak Regression, 							*
*	and Difference-in-Differences Estimators. Wooldrige (2021)	*
*																*
*	Prepared for: Nixon Candiales								*
*	e-mail: natorre@emory.edu									*
*																*	
*	Prepared to: Sara Markowitz									*
*	smarko2@emory.edu											*
*	Emory University											*
*	02/27/2020													*
*****************************************************************

**# Initial Setup 

/****** Set-up ************/

set more off
set varabbrev off
set linesize 132


/****** Setting directory ************/

if `"`c(os)'"' == "MacOSX" & `"`c(username)'"' == "nix" cap cd "~/Documents/Github/RA_FALL22" 
//Collaborator 1: MAC OS
if `"`c(os)'"' == "Windows" & `"`c(username)'"' == "16083" cap cd "~/Documents/Github/RA_FALL22" 
//Collaborator 2: Windows Laptop
if `"`c(os)'"' == "Windows" & `"`c(username)'"' == "natorre" cap cd "C:\Users\natorre\Documents\GitHub\RA_Fall22" 
//Collaborator 3: Destokp Office

log using Logs/Replication_1, text replace

/****** Loading the data ************/
use Data/analy_malp_paper_9_22.dta, clear
xtset stfips year // Declaring it as a Panel

/* Define locals controls, torts, atorts */
local x pctpoverty unemp rpcinc
local torts jsl pcap necap apology
local atorts pcap necap apology

/****** Define variables as in Woldrige (2021) ************/

/* d=1 if recieve the treatment. */
gen d = (!missing(FPA_FULL_YEAR))
label variable d "=1 if treated"

/* p=1 if post-treatment time */
gen p = (pre==0)
label variable p "=1 if post-treatment period"

/* w=d*p */
gen w= (d*p)
label variable w "=1 if treated & in post-treatment period"

/* year indicators */
qui tab year, gen(f)
rename f# f#, renumber(1998) // starting year in 1998
unab i_year : f1* f2*		// varlist of year indicators

/* treatment cohort indicators */
qui tab FPA_FULL_YEAR, gen(c)
unab cohort : c*			// varlist of cohort indicators

/* replacing missing values */
foreach v of varlist c* {
    qui replace `v' = 0 if `v' == .
}


/* Generate demeaned variables
Generate new variables centered around cohort means
Note to self:   THIS MUST BE DONE OVER CONSTANT VARIABLES !!! 
				That is the reason I had a problem with my attemp. 
				I was incluing a linear combination of the controls!!!!!!!!!
				As a result the algorithm was not converging 

forvalues i = 1/14 {
		local x_dm`i'
		local y "x_dm`i'"
		foreach v in `x' {
			sum `v' if c`i'
			gen `v'_dm`i' = `v' - r(mean)
			label variable `v'_dm`i' "`v' demeaned at cohort `i'"
			local `y' `x_dm`i'' `v'_dm`i'
		}
}

*/

**# Main model - Poisson staggered entry case

/*  Estimate de model  
The first set of regresors are interactions between treatment,cohort and post-treatment indicators;
The second set of regresors are the interactions between treatment,cohort,post-treatment and controls centered around cohort means (CONSTANT CONTROLS) indicators;
The third set of regressors are indicators for year, cohort and interactions between year, cohort and controls.
*/

qui poisson malp_np i.w#c.c1#c.f1998-f2019 i.w#c.c2#c.f2000-f2019 i.w#c.c3#c.f2005-f2019 ///
				i.w#c.c4#c.f2006-f2019 i.w#c.c5#c.f2010-f2019 i.w#c.c6#c.f2011-f2019 ///
				i.w#c.c7#c.f2012-f2019 i.w#c.c8#c.f2013-f2019 i.w#c.c9#c.f2014-f2019 ///
				i.w#c.c10#c.f2015-f2019 i.w#c.c11#c.f2016-f2019 i.w#c.c12#c.f2018-f2019 i.w#c.c13#c.f2019 ///
				`i_year' c.`i_year'#c.`x' ///
				`cohort' `x' c.c*#c.`x' , noomitted vce(cluster stfips)
estimates store model1
				/*	To include if we have constant controls over time
				///
				w#c.c1#c.f1998-f2019#c.`x_dm1' w#c.c2#c.f2000-f2019#c.`x_dm2' ///
				w#c.c3#c.f2005-f2019#c.`x_dm3' w#c.c4#c.f2006-f2019#c.`x_dm4' /// 
				w#c.c5#c.f2010-f2019#c.`x_dm5' w#c.c6#c.f2011-f2019#c.`x_dm6' /// 
				w#c.c7#c.f2012-f2019#c.`x_dm7' w#c.c8#c.f2013-f2019#c.`x_dm8' ///
				w#c.c9#c.f2014-f2019#c.`x_dm9' w#c.c10#c.f2015-f2019#c.`x_dm10' ///
				w#c.c11#c.f2016-f2019#c.`x_dm11' w#c.c12#c.f2018-f2019#c.`x_dm12'/// 
				w#c.c13#c.f2019#c.`x_dm13' ///
				*/
				
**# ATT 
/* To obtain ATTs using margins.

As in CS (2021) and Sun and Abraham (2021), the treatment effects to identify are the ATTs in periods 
where the cohorts are actually subjected to the intervention. <- ATT come from pooled Poisson.->

We have to include w explicitly. Also, have to 
set the cohort (c*)and year (i_year) indicators  accordingly

The following is the example provided by Wooldrige. Using Margins after possion 
estimates the ATT at the cohort q in time t.

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
	
I tried to approach this ATT in a iterative way 
to get reproducible results. 

In a simpler way I generate some inteactions for the cohort q in time t and indclude those in the model. 
See next chunk.	
*/

* Generate  the interactions
numlist "1/14"
local y1 = "`r(numlist)'"

local n : word count `y1'
levelsof FPA_FULL_YEAR, local(y2)

 forvalues i = 1/`n' {
      local a : word `i' of `y1'
      local b : word `i' of `y2'
	  forvalues j = `b'/2019{
		qui gen c`a'f`j' = c`a'*f`j'
		qui label variable c`a'f`j' "Cohort `a'"
	  }
   }

*--------------------------------------------------------------------------------------------
*------------------------- Estimate the model DEP VAR: mlp_np *------------------------------
*--------------------------------------------------------------------------------------------
qui poisson malp_np i.c1f* i.c2f* i.c3f* i.c4f* i.c5f* i.c6f* i.c7f* i.c8f* i.c9f* i.c10f* i.c11f* i.c12f* i.c13f* `i_year' ///
					c.`i_year'#c.`x' `cohort' `x' c.c*#c.`x' , noomitted vce(cluster stfips)
					/* To include if we have constant controls over time
					i.c1f*#c.x_dm1 i.c2f*#c.x_dm3 i.c3f*#c.x_dm3 .... and so on for the demeaned variables
					*/
estimates store model2 //store the estimates

* Example How to get the ATT at cohort 13 at time 2019
	estimates restore model2
	margins, dydx(1.c13f2019) subpop(if c13f2019== 1)  noestimcheck vce(uncond)	
	margins, eydx(1.c13f2019) subpop(if c13f2019== 1)  noestimcheck vce(uncond)	
*	scalar tau13_2019 = _b[1.c13f2019]
*	di "tau13_2019 " tau13_2019

/*	* Iterative approach to get the ATT over all cohorts and post-treatment period.
*/

/* Average Marginal Effect */
numlist "1/14"
local y1 = "`r(numlist)'"
local n : word count `y1'
levelsof FPA_FULL_YEAR, local(y2)
postfile DyDx_mlp cohort time b se z pvalue ll ul using DyDx_mlp, replace 
 forvalues i = 1/`n' {
      local a : word `i' of `y1'
      local b : word `i' of `y2'
	  di "Estimating Average Marginal Effects - Cohort `a'" _cont
	  forvalues j = `b'/2019{
	  	qui estimates restore model2
		qui margins, dydx(1.c`a'f`j') subpop(if c`a'f`j'== 1) noestimcheck vce(uncond)
		local b = r(table)[1,2]
		local se = r(table)[2,2]
		local z = r(table)[3,2]
		local pvalue = r(table)[4,2]
		local ll = r(table)[5,2]
		local ul = r(table)[6,2]
		post DyDx_mlp (`a') (`j') (`b') (`se') (`z') (`pvalue') (`ll') (`ul')
		di "." _cont
	  }
	  	noi di _dup(10) "." _cont
		di "Done" _cont 	
		noi di _dup(1) "!"
   }
postclose DyDx_mlp


/* Semi Elasticity */
numlist "1/14"
local y1 = "`r(numlist)'"
local n : word count `y1'
levelsof FPA_FULL_YEAR, local(y2)
postfile EyDx_mlp cohort time b se z pvalue ll ul using EyDx_mlp, replace 
 forvalues i = 1/`n' {
      local a : word `i' of `y1'
      local b : word `i' of `y2'
	  di "Estimating Predicted Number of Events - Cohort `a'" _cont
	  forvalues j = `b'/2019{
	  	qui estimates restore model2
		qui margins, eydx(1.c`a'f`j') subpop(if c`a'f`j'== 1) noestimcheck vce(uncond)
		local b = r(table)[1,2]
		local se = r(table)[2,2]
		local z = r(table)[3,2]
		local pvalue = r(table)[4,2]
		local ll = r(table)[5,2]
		local ul = r(table)[6,2]
		post EyDx_mlp (`a') (`j') (`b') (`se') (`z') (`pvalue') (`ll') (`ul')
		di "." _cont
	  }
	  	noi di _dup(10) "." _cont
		di "Done" _cont 	
		noi di _dup(1) "!"
   }
postclose EyDx_mlp

*--------------------------------------------------------------------------------------------
*------------------------- Estimate the model DEP VAR: mlp_phy *------------------------------
*--------------------------------------------------------------------------------------------
qui poisson malp_phys i.c1f* i.c2f* i.c3f* i.c4f* i.c5f* i.c6f* i.c7f* i.c8f* i.c9f* i.c10f* i.c11f* i.c12f* i.c13f* `i_year' ///
					c.`i_year'#c.`x' `cohort' `x' c.c*#c.`x' , noomitted vce(cluster stfips)
					/* To include if we have constant controls over time
					i.c1f*#c.x_dm1 i.c2f*#c.x_dm3 i.c3f*#c.x_dm3 .... and so on for the demeaned variables
					*/
estimates store model3 //store the estimates

* Example How to get the ATT at cohort 13 at time 2019
	estimates restore model3
	margins, dydx(1.c13f2019) subpop(if c13f2019== 1)  noestimcheck vce(uncond)	
	margins, eydx(1.c13f2019) subpop(if c13f2019== 1)  noestimcheck vce(uncond)	
*	scalar tau13_2019 = _b[1.c13f2019]
*	di "tau13_2019 " tau13_2019

/*	* Iterative approach to get the ATT over all cohorts and post-treatment period.
*/

/* Average Marginal Effect */
numlist "1/14"
local y1 = "`r(numlist)'"
local n : word count `y1'
levelsof FPA_FULL_YEAR, local(y2)
postfile DyDx_mlp_phys cohort time b se z pvalue ll ul using DyDx_mlp_phys, replace 
 forvalues i = 1/`n' {
      local a : word `i' of `y1'
      local b : word `i' of `y2'
	  di "Estimating Average Marginal Effects - Cohort `a'" _cont
	  forvalues j = `b'/2019{
	  	qui estimates restore model3
		qui margins, dydx(1.c`a'f`j') subpop(if c`a'f`j'== 1) noestimcheck vce(uncond)
		local b = r(table)[1,2]
		local se = r(table)[2,2]
		local z = r(table)[3,2]
		local pvalue = r(table)[4,2]
		local ll = r(table)[5,2]
		local ul = r(table)[6,2]
		post DyDx_mlp_phys (`a') (`j') (`b') (`se') (`z') (`pvalue') (`ll') (`ul')
		di "." _cont
	  }
	  	noi di _dup(10) "." _cont
		di "Done" _cont 	
		noi di _dup(1) "!"
   }
postclose DyDx_mlp_phys


/* Semi Elasticity */
numlist "1/14"
local y1 = "`r(numlist)'"
local n : word count `y1'
levelsof FPA_FULL_YEAR, local(y2)
postfile EyDx_mlp_phys cohort time b se z pvalue ll ul using EyDx_mlp_phys, replace 
 forvalues i = 1/`n' {
      local a : word `i' of `y1'
      local b : word `i' of `y2'
	  di "Estimating Predicted Number of Events - Cohort `a'" _cont
	  forvalues j = `b'/2019{
	  	qui estimates restore model3
		qui margins, eydx(1.c`a'f`j') subpop(if c`a'f`j'== 1) noestimcheck vce(uncond)
		local b = r(table)[1,2]
		local se = r(table)[2,2]
		local z = r(table)[3,2]
		local pvalue = r(table)[4,2]
		local ll = r(table)[5,2]
		local ul = r(table)[6,2]
		post EyDx_mlp_phys (`a') (`j') (`b') (`se') (`z') (`pvalue') (`ll') (`ul')
		di "." _cont
	  }
	  	noi di _dup(10) "." _cont
		di "Done" _cont 	
		noi di _dup(1) "!"
   }
postclose EyDx_mlp_phys

   

**# Desired Output
/* Goal 
The goal of this do-file is Implementing Two-Way Fixed Effects,							
	the Two-Way Mundlak Regression, 							
	and Difference-in-Differences Estimators proposed by Wooldrige (2021)

	That is estimate the ATT per cohort and present the results in a table. 
	Also, include code to do a linear combination of ATT to get a unique effect (such an ATE)
*/

*** An example of the disired outcome is this 
*** Note I am restricting subpop(if c`i'== 1) intentioally
*** Just to get a sense of how the table should look like. 
*** The ATT estimated are not correct. Toy example. 
estimates restore model2
forvalues i = 1/13 {
	
	foreach var in c`i'f*	{
		estimates restore model2
		 margins, dydx(1.`var') subpop(if c`i'== 1) noestimcheck vce(uncond)
	}		
}		

/* Callaway and Sant'Anna (2020, Journal of Econometrics) */
*** Note I am not including covariates intentioally
*** Just to get a sense of how the table should look like. 
*** The disire outcome should append the previous tables as
gen first_treat = FPA_FULL_YEAR
replace first_treat = 0 if first_treat == .
csdid malp_np , ivar(stfips) time(year) gvar(first_treat)

   
**# To do 
/*

Reproduce the Test for non-Parallel trends in linear model. It rejects with p < 0.05,
indicating a problem with the linear model:

Try Callaway and Sant'Anna (2020, Journal of Econometrics)

gen first_treat = FPA_FULL_YEAR
replace first_treat = 0 if first_treat == .
csdid malp_np `x' , ivar(stfips) time(year) gvar(first_treat)
Check why it is not working when adding controls!!!!!!


Construct a table for ATT, and present a tabe for ATE? (or any other way of grouping them)

Make the code pretty :)
*/

/// Pruebas

* Example How to get the ATT at cohort 13 at time 2019
	estimates restore model2
	margins, dydx(1.c13f2019) subpop(if c13f2019== 1)  noestimcheck vce(uncond)
	noi di _dup(80) "-"
	local b = r(table)[1,2]
	local se = r(table)[2,2]
	local z = r(table)[3,2]
	local pvalue = r(table)[4,2]
	local ll = r(table)[5,2]
	local ul = r(table)[6,2]
	di `b' `se' `ze' `pvalue' `ll' `ul'
	scalar tau13_2019 = _b[1.c13f2019]
	di "tau13_2019 " tau13_2019


local x pctpoverty unemp rpcinc
local torts jsl pcap necap apology
local atorts pcap necap apology
jwdid malp_np , ivar(stfips) tvar(year) gvar(FPA_FULL_YEAR) method(poisson)
qui gen first_treat = FPA_FULL_YEAR
qui replace first_treat = 0 if first_treat == .
csdid malp_np , ivar(stfips) time(year) gvar(first_treat)

use https://friosavila.github.io/playingwithstata/drdid/mpdta.dta, clear
gen emp=exp(lemp)
jwdid  emp lpop, ivar(countyreal) tvar(year) gvar(first_treat) method(poisson)
estat event

	// create the results matrix
matrix A = J(4, 3, .z)

// now fill the matrix
local row = 1
foreach y in MRtb MRn {
    diagt fastib `y'
    matrix A[`row', 1] = r(sens)
    matrix A[`row', 2] = r(spec)
    local ++row
}


postfile margin b se z pvalue ll ul using prueba, replace 
	estimates restore model2
	margins, dydx(1.c13f2019) subpop(if c13f2019== 1)  noestimcheck vce(uncond)
	noi di _dup(80) "-"

post margin (`b') (`se') (`z') (`pvalue') (`ll') (`ul')
postclose margin
use prueba, clear
list



numlist "1/14"
local y1 = "`r(numlist)'"
local n : word count `y1'
levelsof FPA_FULL_YEAR, local(y2)
postfile DyDx_mlp cohort time b se z pvalue ll ul using DyDx_mlp, replace 
 forvalues i = 1/`n' {
      local a : word `i' of `y1'
      local b : word `i' of `y2'
	  di "Estimating Average Marginal Effects - Cohort `a'" _cont
	  forvalues j = `b'/2019{
	  	qui estimates restore model2
		qui margins, dydx(1.c`a'f`j') subpop(if c`a'f`j'== 1) noestimcheck vce(uncond)
		local b = r(table)[1,2]
		local se = r(table)[2,2]
		local z = r(table)[3,2]
		local pvalue = r(table)[4,2]
		local ll = r(table)[5,2]
		local ul = r(table)[6,2]
		post DyDx_mlp (`a') (`j') (`b') (`se') (`z') (`pvalue') (`ll') (`ul')
		di "." _cont
	  }
	  	noi di _dup(10) "." _cont
		di "Done" _cont 	
		noi di _dup(1) "!"
   }
postclose DyDx_mlp



numlist "1/14"
local y1 = "`r(numlist)'"
local n : word count `y1'
levelsof FPA_FULL_YEAR, local(y2)
postfile EyDx_mlp cohort time b se z pvalue ll ul using EyDx_mlp, replace 
 forvalues i = 1/`n' {
      local a : word `i' of `y1'
      local b : word `i' of `y2'
	  di "Estimating Predicted Number of Events - Cohort `a'" _cont
	  forvalues j = `b'/2019{
	  	qui estimates restore model2
		qui margins, eydx(1.c`a'f`j') subpop(if c`a'f`j'== 1) noestimcheck vce(uncond)
		local b = r(table)[1,2]
		local se = r(table)[2,2]
		local z = r(table)[3,2]
		local pvalue = r(table)[4,2]
		local ll = r(table)[5,2]
		local ul = r(table)[6,2]
		post EyDx_mlp (`a') (`j') (`b') (`se') (`z') (`pvalue') (`ll') (`ul')
		di "." _cont
	  }
	  	noi di _dup(10) "." _cont
		di "Done" _cont 	
		noi di _dup(1) "!"
   }
postclose EyDx_mlp


eydx

use EyDx_mlp, clear
list
use EyDx_mlp_phys, clear
list

jwdid 

log close
