set more off
set varabbrev off
set linesize 132


if `"`c(os)'"' == "MacOSX" & `"`c(username)'"' == "nix" cap cd "~/Documents/Github/RA_FALL22" 
//Collaborator 1: MAC OS
 do Code/Replication_paper_SM

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














****** Adapting Woldrige Code ************

/* Logit */

logit y i.w#c.d#c.f04 i.w#c.d#c.f05 i.w#c.d#c.f06 ///
	i.w#c.d#c.f04#c.x i.w#c.d#c.f05#c.x i.w#c.d#c.f06#c.x ///
	f02 f03 f04 f05 f06 c.f02#c.x c.f03#c.x c.f04#c.x c.f05#c.x c.f06#c.x ///
	d x c.d#c.x, noomitted vce(cluster id)

/* poisson */

logit malp_np fpa  $torts $control i.stfips i.year if Mstclass<3, noomitted vce(cluster stfips)





























