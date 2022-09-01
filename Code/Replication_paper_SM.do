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

poisson  `dvar' fpa $atorts $control i.stfips i.year if Astclass<3, exp(pop)                 cluster(stfips)

}