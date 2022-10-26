remotes::install_github("grantmcdermott/etwfe")

library(etwfe)
library(haven)
library(vroom)
library(tidyverse)

dat <- read_stata("C:/Users/16083/Documents/GitHub/RA_Fall22/Data/analy_malp_paper_9_22.dta")

dat <- dat %>% 
        group_by(stfips) %>% 
        mutate(lpop = log(pop),
               lpop.bar = mean(lpop))

dat.reg <- dat %>% 
                mutate(first_treat = ifelse(is.na(FPA_FULL_YEAR),0,FPA_FULL_YEAR)) %>% 
                filter(Mstclass<3)

mod = 
  etwfe(
    fml  = malp_phys ~ lpop,    # outcome ~ controls
    tvar = year,                # time variable
    gvar = first_treat,         # group variable
    data = dat.reg,             # dataset
    vcov = ~ stfips,            # vcov adjustment (here: clustered)
    family = poisson            # poisson  regression
    #ivar = stfips 
    )

emfx(mod, type = "group")


mod2 = 
  etwfe(
    fml  = malp_np ~ lpop,    # outcome ~ controls
    tvar = year,                # time variable
    gvar = first_treat,         # group variable
    data = dat.reg,             # dataset
    vcov = ~ stfips,            # vcov adjustment (here: clustered)
    family = poisson            # poisson  regression
    #ivar = stfips 
    )

emfx(mod2, type = "group", summary = TRUE, slope = "eydx")



##' Post-estimation treatment effects for an ETWFE regressions.
##'
##' @param object An `etwfe` model object.
##' @param type Character. The desired type of post-estimation aggregation.
##' @param summary Logical. Should the resulting marginaleffects objects be passed to [`summary`] before being returned? Defaults to TRUE, but mostly for aesthetics reasons.
##' @param ... Additional arguments passed to 
##' [`marginaleffects::marginaleffects`].
##' @return A marginaleffects object.
##' @seealso [marginaleffects::marginaleffects()]
##' @inherit etwfe return examples
##' @export
emfx = function(
    object,
    type = c("simple", "group", "calendar", "event"),
    summary = TRUE,
    ...
) {
  
  .Dtreat = NULL
  type = match.arg(type)
  gvar = attributes(object)[["etwfe"]][["gvar"]]
  tvar = attributes(object)[["etwfe"]][["tvar"]]
  
  dat = eval(object$call$data, object$call_env)
  if (".Dtreat" %in% names(dat)) dat = dat[dat[[".Dtreat"]]==1L, , drop = FALSE]
  
  if (type=="simple") by_var = ".Dtreat"
  if (type=="group") by_var = gvar
  if (type=="calendar") by_var = tvar
  if (type=="event") {
    dat[["event"]] = dat[[tvar]] - dat[[gvar]]
    by_var = "event"
  }
  
  mfx = marginaleffects::marginaleffects(
    object, 
    newdata = dat,
    variables = ".Dtreat",
    by = by_var,
    ...
  )
  
  if (type!="simple") mfx = mfx[order(mfx[[by_var]]), ]
  
  if (summary) mfx = summary(mfx)
  
  return(mfx)
}