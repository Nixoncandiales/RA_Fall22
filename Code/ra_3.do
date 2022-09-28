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

** Initial Setup 

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

cap ssc install jwdid
cap ssc install reghdfe
cap ssc install ftools

use Data/analy_malp_paper_9_22.dta, clear

global control "pctpoverty unemp rpcinc"

global torts "jsl pcap necap apology"

global atorts "pcap necap apology"

 

 

/*Main models*/

poisson  malp_np   fpa $torts $control i.stfips i.year if Mstclass<3, exp(pop) cluster(stfips)

poisson  malp_phys fpa $torts $control i.stfips i.year if Mstclass<3, exp(pop) cluster(stfips)

 

foreach dvar in adv_np vunsafe_np vRX_np vSOP_np  {

poisson  `dvar' fpa $atorts $control i.stfips i.year if Astclass<3, exp(pop) cluster(stfips)

}

 

foreach dvar in adv_phys vunsafe_phys vRX_phys vSOP_phys{

poisson  `dvar' fpa $atorts $control i.stfips i.year if Astclass<3, exp(pop) cluster(stfips)

}



jwdid malp_np if Mstclass<3, ivar(stfips) tvar(year) gvar(first_treat) method(poisson)
estat group


jwdid malp_phys  if Mstclass<3, ivar(stfips) tvar(year) gvar(first_treat) method(poisson)
estat group

gen first_treat = FPA_FULL_YEAR
replace first_treat = 0 if first_treat == .
csdid malp_np , ivar(stfips) time(year) gvar(first_treat)