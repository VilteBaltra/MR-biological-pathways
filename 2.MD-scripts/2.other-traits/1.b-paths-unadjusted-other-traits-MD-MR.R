
# This script obtains univariable (unadjusted) MR estimates for cholesterol traits, glycaemic traits, cortisol and CRP
# trait IDs taken from https://gwas.mrcieu.ac.uk/ 

library(tidyverse)
library(TwoSampleMR) 
library(ieugwasr)
library(cowplot)

# path
setwd("~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/DEPRESSION/other-traits/")

#---------------------------------#
#       B PATH ESTIMATES
#---------------------------------#

# define b path mr function
my.mr.b.path <- function(exposure.data, mediator){
  out<- read_outcome_data(
    snps = exposure.data$SNP,
    filename = paste0("~/Documents/Projects/MR-mediation-CM-MM/summary-stats/daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz"),
    sep = "\t", snp_col = "SNP",
    beta_col = "logOR", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "FRQ_A_116209", #  eaf FRQ_A_116209	or FRQ_U_314566 ?
    pval_col = "P")
  out$outcome <- "Depression"
  
  dat <- harmonise_data(
    exposure_dat = exposure.data, 
    outcome_dat = out)
  
  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  
  # horizontal pleiotropy
  intercept <- mr_pleiotropy_test(dat)
  
  # heterogeneity statistic
  Q <- mr_heterogeneity(dat)
  
  # scatter plot
  p1 <- mr_scatter_plot(res, dat)[[1]]
  # forest plot
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot(res_single)[[1]]
  # Leave-one-out plot
  res_loo <- mr_leaveoneout(dat)
  p3 <- mr_leaveoneout_plot(res_loo)[[1]]
  # Funnel plot
  res_single <- mr_singlesnp(dat)
  p4 <- mr_funnel_plot(res_single)[[1]]
  
  # Arrange the plots in a grid
  combined_plot <- plot_grid(p1, p2, p3, p4, ncol = 2)
  
  # Save the combined plot
  ggsave(paste0(mediator, "-MR-plot.png"), combined_plot, dpi = 600, width = 15, height = 20, units = "in")
  
  return(list(res, intercept, Q))
}

### Cholesterol traits ###
# LDL --> MD
ldl_exp_dat <- extract_instruments(outcomes = "ieu-a-300") # same as in MR from 2013 GLGC
res_ldl <- my.mr.b.path(ldl_exp_dat, mediator="ldl-ieu-a-300") 
ldl_exp_dat2 <- extract_instruments(outcomes = "ebi-a-GCST90018961") 
res_ldl2<- my.mr.b.path(ldl_exp_dat2, mediator="ldl-ebi-a-GCST90018961") 
# total cholesterol --> MD
tc_exp_dat <- extract_instruments(outcomes = "ieu-a-301") # same as in MR from 2013 GLGC
res_tc <- my.mr.b.path(tc_exp_dat, mediator="tc-ieu-a-301") 
tc_exp_dat2 <- extract_instruments(outcomes =  "ebi-a-GCST90018974") 
res_tc2 <- my.mr.b.path(tc_exp_dat2, mediator="tc-ebi-a-GCST90018974")
# triglycerides --> MD
trig_exp_dat <- extract_instruments(outcomes =  "ieu-a-302") # same as in MR from 2013 GLGC
res_trig <- my.mr.b.path(trig_exp_dat, mediator="tg-ieu-a-302") 
trig_exp_dat2 <- extract_instruments(outcomes =  "ieu-b-111") # ukbb
res_trig2 <- my.mr.b.path(trig_exp_dat2, mediator="tg-ieu-b-111") 
# HDL --> MD
hdl_exp_dat <- extract_instruments(outcomes = "ieu-a-299") # same as in MR from 2013 GLGC
res_hdl <- my.mr.b.path(hdl_exp_dat, mediator="hdl-ieu-a-299") 
hdl_exp_dat2 <- extract_instruments(outcomes = "ieu-b-109") # ukbb
res_hdl2 <- my.mr.b.path(hdl_exp_dat2, mediator = "hdl-ieu-b-109") 

### Glycaemic traits ### 
# FI --> MD
fi_exp_dat <- extract_instruments(outcomes = "ieu-b-116") 
res_fi <- my.mr.b.path(fi_exp_dat, mediator = "fi-ieu-b-116") 
fi_exp_dat2 <- extract_instruments(outcomes = "ebi-a-GCST90002238")
res_fi2 <- my.mr.b.path(fi_exp_dat2, mediator = "fi-ebi-a-GCST90002238" )
# fasting glucose --> MD
fg_exp_dat <- extract_instruments(outcomes = "ebi-a-GCST90002232") 
res_fg <- my.mr.b.path(fg_exp_dat, mediator = "fg-ebi-a-GCST90002232")
fg_exp_dat2 <- extract_instruments(outcomes = "ieu-b-113") 
res_fg2 <- my.mr.b.path(fg_exp_dat2, mediator = "fg-ieu-b-113")
# HbA1C --> MD
HbA1C_exp_dat <- extract_instruments(outcomes = "ebi-a-GCST90014006") 
res_HbA1C <- my.mr.b.path(HbA1C_exp_dat, mediator = "hba1c-ebi-a-GCST90014006")
HbA1C_exp_dat2 <- extract_instruments(outcomes = "ieu-b-4842") 
res_HbA1C2 <- my.mr.b.path(HbA1C_exp_dat2, mediator = "hba1c-ieu-b-4842") 

### Cortisol ### 
cortisol_exp_dat <- extract_instruments(outcomes = "ieu-a-1012", p1 = 5e-06, p2 = 5e-06) # only 1 snp with 5e-08
res_cortisol <- my.mr.b.path(cortisol_exp_dat, mediator = "cortisol-ieu-a-1012")

#---------------------------------#
#            SAVE OUTPUT
#---------------------------------#
# combine all output into one (b paths)
out_b_paths <- rbind(res_ldl[[1]], res_ldl2[[1]], res_tc[[1]], res_tc2[[1]], res_hdl[[1]], res_hdl2[[1]], res_trig[[1]], res_trig2[[1]], # cholesterol traits 
                     res_fi[[1]], res_fi2[[1]], res_fg[[1]], res_fg2[[1]], res_HbA1C[[1]], res_HbA1C2[[1]], # glycaemic traits 
                     res_cortisol[[1]]) # cortisol
      
# view significant effects
out_b_paths %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05) # only fasting insulin (ebi-a-GCST90002238)

out_intercepts <- rbind(res_ldl[[2]], res_ldl2[[2]], res_tc[[2]], res_tc2[[2]], res_hdl[[2]], res_hdl2[[2]], res_trig[[2]], res_trig2[[2]], # cholesterol traits 
                     res_fi[[2]], res_fi2[[2]], res_fg[[2]], res_fg2[[2]], res_HbA1C[[2]], res_HbA1C2[[2]],  # glycaemic traits 
                     res_cortisol[[2]]) # cortisol
     

out_Q <- rbind(res_ldl[[3]], res_ldl2[[3]], res_tc[[3]], res_tc2[[3]], res_hdl[[3]], res_hdl2[[3]], res_trig[[3]], res_trig2[[3]], # cholesterol traits 
                        res_fi[[3]], res_fi2[[3]], res_fg[[3]], res_fg2[[3]], res_HbA1C[[3]], res_HbA1C2[[3]], # glycaemic traits 
                        res_cortisol[[3]]) # cortisol

# save all
out.list <- list(out_b_paths, out_intercepts, out_Q)
names(out.list) <- c("unadjusted b paths", "Egger intercepts", "Q statistics") 
openxlsx::write.xlsx(out.list, file = paste0("b-paths-unadjusted-Egger-Q-stat-uni-MR-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)


