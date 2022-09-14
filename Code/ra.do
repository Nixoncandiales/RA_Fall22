set more off
set varabbrev off
set linesize 132


if `"`c(os)'"' == "MacOSX" & `"`c(username)'"' == "nix" cap cd "~/Documents/Github/RA_FALL22" 
//Collaborator 1: MAC OS
if `"`c(os)'"' == "Windows" & `"`c(username)'"' == "16083" cap cd "~/Documents/Github/RA_FALL22" 
 do Code/Replication_paper_SM

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

** Generate demeaned variables by cohort

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


/* poisson */

poisson malp_np i.w#c.c1#c.f1998 i.w#c.c1#c.f1999 i.w#c.c1#c.f2000 ///
		i.w#c.c1#c.f2001 i.w#c.c1#c.f2002 i.w#c.c1#c.f2003 i.w#c.c1#c.f2004 ///
		i.w#c.c1#c.f2005 i.w#c.c1#c.f2006 i.w#c.c1#c.f2007 i.w#c.c1#c.f2008 ///
		i.w#c.c1#c.f2009 i.w#c.c1#c.f2010 i.w#c.c1#c.f2011 i.w#c.c1#c.f2012 ///
		i.w#c.c1#c.f2013 i.w#c.c1#c.f2014 i.w#c.c1#c.f2015 i.w#c.c1#c.f2016 ///
		i.w#c.c1#c.f2017 i.w#c.c1#c.f2018 i.w#c.c1#c.f2019 ///
		///
		i.w#c.c2#c.f2000 i.w#c.c2#c.f2001 i.w#c.c2#c.f2002 i.w#c.c2#c.f2003 ///
		i.w#c.c2#c.f2004 i.w#c.c2#c.f2005 i.w#c.c2#c.f2006 i.w#c.c2#c.f2007 ///
		i.w#c.c2#c.f2008 i.w#c.c2#c.f2009 i.w#c.c2#c.f2010 i.w#c.c2#c.f2011 ///
		i.w#c.c2#c.f2012 i.w#c.c2#c.f2013 i.w#c.c2#c.f2014 i.w#c.c2#c.f2015 ///
		i.w#c.c2#c.f2016 i.w#c.c2#c.f2017 i.w#c.c2#c.f2018 i.w#c.c2#c.f2019 ///
		///
		i.w#c.c3#c.f2005 i.w#c.c3#c.f2006 i.w#c.c3#c.f2007 i.w#c.c3#c.f2008 ///
		i.w#c.c3#c.f2009 i.w#c.c3#c.f2010 i.w#c.c3#c.f2011 i.w#c.c3#c.f2012 ///
		i.w#c.c3#c.f2013 i.w#c.c3#c.f2014 i.w#c.c3#c.f2015 i.w#c.c3#c.f2016 ///
		i.w#c.c3#c.f2017 i.w#c.c3#c.f2018 i.w#c.c3#c.f2019 ///		
		///
		i.w#c.c4#c.f2006 i.w#c.c4#c.f2007 i.w#c.c4#c.f2008 i.w#c.c4#c.f2009 ///
		i.w#c.c4#c.f2010 i.w#c.c4#c.f2011 i.w#c.c4#c.f2012 i.w#c.c4#c.f2013 ///
		i.w#c.c4#c.f2014 i.w#c.c4#c.f2015 i.w#c.c4#c.f2016 i.w#c.c4#c.f2017 ///
		i.w#c.c4#c.f2018 i.w#c.c4#c.f2019 ///		
		///
		i.w#c.c5#c.f2010 i.w#c.c5#c.f2011 i.w#c.c5#c.f2012 i.w#c.c5#c.f2013 ///
		i.w#c.c5#c.f2014 i.w#c.c5#c.f2015 i.w#c.c5#c.f2016 i.w#c.c5#c.f2017 ///
		i.w#c.c5#c.f2018 i.w#c.c5#c.f2019 ///
		///
		i.w#c.c6#c.f2011 i.w#c.c6#c.f2012 i.w#c.c6#c.f2013 i.w#c.c6#c.f2014 ///
		i.w#c.c6#c.f2015 i.w#c.c6#c.f2016 i.w#c.c6#c.f2017 i.w#c.c6#c.f2018 ///
		i.w#c.c6#c.f2019 ///
		///
		i.w#c.c7#c.f2012 i.w#c.c7#c.f2013 i.w#c.c7#c.f2014 i.w#c.c7#c.f2015 ///
		i.w#c.c7#c.f2016 i.w#c.c7#c.f2017 i.w#c.c7#c.f2018 i.w#c.c7#c.f2019 ///
		///
		i.w#c.c8#c.f2013 i.w#c.c8#c.f2014 i.w#c.c8#c.f2015 i.w#c.c8#c.f2016 ///
		i.w#c.c8#c.f2017 i.w#c.c8#c.f2018 i.w#c.c8#c.f2019 ///
		///
		i.w#c.c9#c.f2014 i.w#c.c9#c.f2015 i.w#c.c9#c.f2016 i.w#c.c9#c.f2017 ///
		i.w#c.c9#c.f2018 i.w#c.c9#c.f2019 ///
		///
		i.w#c.c10#c.f2015 i.w#c.c10#c.f2016 i.w#c.c10#c.f2017 ///
		i.w#c.c10#c.f2018 i.w#c.c10#c.f2019 ///
		///
		i.w#c.c11#c.f2016 i.w#c.c11#c.f2017 ///
		i.w#c.c11#c.f2018 i.w#c.c11#c.f2019 ///
		///
		i.w#c.c12#c.f2018 i.w#c.c12#c.f2019 ///
		///
		i.w#c.c13#c.f2019 ///
		i.w#c.c1#c.f1998#c.`control_dm_c1' i.w#c.c1#c.f1999#c.`control_dm_c1' i.w#c.c1#c.f2000#c.`control_dm_c1' ///
		i.w#c.c1#c.f2001#c.`control_dm_c1' i.w#c.c1#c.f2002#c.`control_dm_c1' i.w#c.c1#c.f2003#c.`control_dm_c1' ///
		i.w#c.c1#c.f2004#c.`control_dm_c1' i.w#c.c1#c.f2005#c.`control_dm_c1' i.w#c.c1#c.f2006#c.`control_dm_c1' ///
		i.w#c.c1#c.f2007#c.`control_dm_c1' i.w#c.c1#c.f2008#c.`control_dm_c1' i.w#c.c1#c.f2009#c.`control_dm_c1' ///
		i.w#c.c1#c.f2010#c.`control_dm_c1' i.w#c.c1#c.f2011#c.`control_dm_c1' i.w#c.c1#c.f2012#c.`control_dm_c1' ///
		i.w#c.c1#c.f2013#c.`control_dm_c1' i.w#c.c1#c.f2014#c.`control_dm_c1' i.w#c.c1#c.f2015#c.`control_dm_c1' ///
		i.w#c.c1#c.f2016#c.`control_dm_c1' i.w#c.c1#c.f2017#c.`control_dm_c1' i.w#c.c1#c.f2018#c.`control_dm_c1' ///
		i.w#c.c1#c.f2019#c.`control_dm_c1' ///
		///
		i.w#c.c2#c.f2000#c.`control_dm_c2' i.w#c.c2#c.f2001#c.`control_dm_c2' i.w#c.c2#c.f2002#c.`control_dm_c2' ///
		i.w#c.c2#c.f2003#c.`control_dm_c2' i.w#c.c2#c.f2004#c.`control_dm_c2' i.w#c.c2#c.f2005#c.`control_dm_c2' ///
		i.w#c.c2#c.f2006#c.`control_dm_c2' i.w#c.c2#c.f2007#c.`control_dm_c2' i.w#c.c2#c.f2008#c.`control_dm_c2' ///
		i.w#c.c2#c.f2009#c.`control_dm_c2' i.w#c.c2#c.f2010#c.`control_dm_c2' i.w#c.c2#c.f2011#c.`control_dm_c2' ///
		i.w#c.c2#c.f2012#c.`control_dm_c2' i.w#c.c2#c.f2013#c.`control_dm_c2' i.w#c.c2#c.f2014#c.`control_dm_c2' ///
		i.w#c.c2#c.f2015#c.`control_dm_c2' i.w#c.c2#c.f2016#c.`control_dm_c2' i.w#c.c2#c.f2017#c.`control_dm_c2' ///
		i.w#c.c2#c.f2018#c.`control_dm_c2' i.w#c.c2#c.f2019#c.`control_dm_c2' ///
		///
		i.w#c.c3#c.f2005#c.`control_dm_c3' i.w#c.c3#c.f2006#c.`control_dm_c3' i.w#c.c3#c.f2007#c.`control_dm_c3' ///
		i.w#c.c3#c.f2008#c.`control_dm_c3' i.w#c.c3#c.f2009#c.`control_dm_c3' i.w#c.c3#c.f2010#c.`control_dm_c3' ///
		i.w#c.c3#c.f2011#c.`control_dm_c3' i.w#c.c3#c.f2012#c.`control_dm_c3' i.w#c.c3#c.f2013#c.`control_dm_c3' ///
		i.w#c.c3#c.f2014#c.`control_dm_c3' i.w#c.c3#c.f2015#c.`control_dm_c3' i.w#c.c3#c.f2016#c.`control_dm_c3' 
		i.w#c.c3#c.f2017#c.`control_dm_c3' i.w#c.c3#c.f2018#c.`control_dm_c3' i.w#c.c3#c.f2019#c.`control_dm_c3' /// When Running this part my computer stops working!!!!		
		///
		i.w#c.c4#c.f2006#c.`control_dm_c4' i.w#c.c4#c.f2007#c.`control_dm_c4' i.w#c.c4#c.f2008#c.`control_dm_c4' ///
		i.w#c.c4#c.f2009#c.`control_dm_c4' i.w#c.c4#c.f2010#c.`control_dm_c4' i.w#c.c4#c.f2011#c.`control_dm_c4' ///
		i.w#c.c4#c.f2012#c.`control_dm_c4' i.w#c.c4#c.f2013#c.`control_dm_c4' i.w#c.c4#c.f2014#c.`control_dm_c4' ///
		i.w#c.c4#c.f2015#c.`control_dm_c4' i.w#c.c4#c.f2016#c.`control_dm_c4' i.w#c.c4#c.f2017#c.`control_dm_c4' ///
		i.w#c.c4#c.f2018#c.`control_dm_c4' i.w#c.c4#c.f2019#c.`control_dm_c4' ///		
		///
		i.w#c.c5#c.f2010#c.`control_dm_c5' i.w#c.c5#c.f2011#c.`control_dm_c5' i.w#c.c5#c.f2012#c.`control_dm_c5' ///
		i.w#c.c5#c.f2013#c.`control_dm_c5' i.w#c.c5#c.f2014#c.`control_dm_c5' i.w#c.c5#c.f2015#c.`control_dm_c5' ///
		i.w#c.c5#c.f2016#c.`control_dm_c5' i.w#c.c5#c.f2017#c.`control_dm_c5' i.w#c.c5#c.f2018#c.`control_dm_c5' ///
		i.w#c.c5#c.f2019#c.`control_dm_c5' ///
		///
		i.w#c.c6#c.f2011#c.`control_dm_c6' i.w#c.c6#c.f2012#c.`control_dm_c6' i.w#c.c6#c.f2013#c.`control_dm_c6' ///
		i.w#c.c6#c.f2014#c.`control_dm_c6' i.w#c.c6#c.f2015#c.`control_dm_c6' i.w#c.c6#c.f2016#c.`control_dm_c6' ///
		i.w#c.c6#c.f2017#c.`control_dm_c6' i.w#c.c6#c.f2018#c.`control_dm_c6' i.w#c.c6#c.f2019#c.`control_dm_c6' ///
		///
		i.w#c.c7#c.f2012#c.`control_dm_c7' i.w#c.c7#c.f2013#c.`control_dm_c7' i.w#c.c7#c.f2014#c.`control_dm_c7' ///
		i.w#c.c7#c.f2015#c.`control_dm_c7' i.w#c.c7#c.f2016#c.`control_dm_c7' i.w#c.c7#c.f2017#c.`control_dm_c7' ///
		i.w#c.c7#c.f2018#c.`control_dm_c7' i.w#c.c7#c.f2019#c.`control_dm_c7' ///
		///
		i.w#c.c8#c.f2013#c.`control_dm_c8' i.w#c.c8#c.f2014#c.`control_dm_c8' i.w#c.c8#c.f2015#c.`control_dm_c8' ///
		i.w#c.c8#c.f2016#c.`control_dm_c8' i.w#c.c8#c.f2017#c.`control_dm_c8' i.w#c.c8#c.f2018#c.`control_dm_c8' ///
		i.w#c.c8#c.f2019#c.`control_dm_c8' ///
		///
		i.w#c.c9#c.f2014#c.`control_dm_c9' i.w#c.c9#c.f2015#c.`control_dm_c9' i.w#c.c9#c.f2016#c.`control_dm_c9' ///
		i.w#c.c9#c.f2017#c.`control_dm_c9' i.w#c.c9#c.f2018#c.`control_dm_c9' i.w#c.c9#c.f2019#c.`control_dm_c9' ///
		///
		i.w#c.c10#c.f2015#c.`control_dm_c10' i.w#c.c10#c.f2016#c.`control_dm_c10' i.w#c.c10#c.f2017#c.`control_dm_c10' ///
		i.w#c.c10#c.f2018#c.`control_dm_c10' i.w#c.c10#c.f2019#c.`control_dm_c10' ///
		///
		i.w#c.c11#c.f2016#c.`control_dm_c11' i.w#c.c11#c.f2017#c.`control_dm_c11' ///
		i.w#c.c11#c.f2018#c.`control_dm_c11' i.w#c.c11#c.f2019#c.`control_dm_c11' ///
		///
		i.w#c.c12#c.f2018#c.`control_dm_c12' i.w#c.c12#c.f2019#c.`control_dm_c12' ///
		///
		i.w#c.c13#c.f2019#c.`control_dm_c13' , vce(cluster stfips)
		///
		c.f c.f*#c.`control' c* c.c*#c.`control' c.`control'///	
		
		
	*Wooldrige	
poisson y i.w#c.d4#c.f04 i.w#c.d4#c.f05 i.w#c.d4#c.f06 ///
	i.w#c.d5#c.f05 i.w#c.d5#c.f06 ///
	i.w#c.d6#c.f06 ///
	i.w#c.d4#c.f04#c.x i.w#c.d4#c.f05#c.x i.w#c.d4#c.f06#c.x ///
	i.w#c.d5#c.f05#c.x i.w#c.d5#c.f06#c.x ///
	i.w#c.d6#c.f06#c.x ///
	f02 f03 f04 f05 f06 ///
	c.f02#c.x c.f03#c.x c.f04#c.x c.f05#c.x c.f06#c.x ///
	d4 d5 d6 x c.d4#c.x c.d5#c.x c.d6#c.x, noomitted vce(cluster id)
	
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














