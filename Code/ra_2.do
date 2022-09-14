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

** p=1 if post-treatment time
gen p = (pre==0)
label variable p "=1 if post-treatment period"

** w=d*p
gen w= (d*p)
label variable w "=1 if treated & in post-treatment period"

** year indicators
tab year, gen(f)
rename f# f#, renumber(1998) // starting year in 1998
unab i_year : f1* f2*		// varlist of year indicators

** treatment cohort indicators
tab FPA_FULL_YEAR, gen(c)
unab cohort : c*			// varlist of cohort indicators

** replacing missing values
foreach v of varlist c* {
    replace `v' = 0 if `v' == .
}

** Generate the demeaned variables by treatment cohort

foreach z in `cohort' {

	local x_dm_`z'
	local y "x_dm_`z'"

		foreach v in `x' {

			sum `v' if `z'
			gen `v'_dm_`z' = `v' - r(mean)
			label variable `v'_dm_`z' "`v' demeaned at cohort `z'"
			local `y' `x_dm_`z'' `v'_dm_`z'
		}
}


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

/* Run regressions in a loop */

local list
foreach var of `x_dm1'-`x_dm14' {
	
	local list "`list' `var'"
	poisson malp_np
	
}







