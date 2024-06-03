### MVMR Egger
# from https://rdrr.io/cran/MendelianRandomization/man/mr_mvegger.html
# from https://rdrr.io/cran/MendelianRandomization/f/inst/doc/Vignette_MR.pdf 

my_mvmregger <- function(mvdat) {
  MRMVInputObject <- mr_mvinput(bx = mvdat$exposure_beta,
                              bxse = mvdat$exposure_se,
                              by = mvdat$outcome_beta,
                              byse = mvdat$outcome_se)
  
  # run mvmr mr-egger
  MVEgger <- mr_mvegger(MRMVInputObject, orientate = 2)
  MVEgger_df <- as.data.frame(cbind(exposure = MVEgger$Exposure, b = MVEgger$Estimate, se = MVEgger$StdError.Est, pval = MVEgger$Pvalue.Est))
  MVEgger_df$b <- as.numeric(MVEgger_df$b)
  MVEgger_df$se <- as.numeric(MVEgger_df$se)
  MVEgger_df <- generate_odds_ratios(MVEgger_df)
  return(MVEgger_df)
}
