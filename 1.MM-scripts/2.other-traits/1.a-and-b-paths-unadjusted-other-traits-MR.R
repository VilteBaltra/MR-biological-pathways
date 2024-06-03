
# This script obtains univariable (unadjusted) MR estimates for cholesterol traits, glycaemic traits, and cortisol
# trait IDs taken from https://gwas.mrcieu.ac.uk/ 

library(tidyverse)
library(TwoSampleMR) 
library(ieugwasr)
library(cowplot)

# a path
setwd("~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/MULTIMORBIDITY/other-traits/unadjusted-a-and-b-paths")

#---------------------------------#
#       READ IN CM GWAS
#---------------------------------#
# make sure gwas sumstats are consistently named as my_mvmr function is expecting! 
# also note that sometimes mvmr function does not run due to traffic (?), but ends up running after a few tries 
# Because of this I edited it to use ld_clump(), instead of remote clumping function. Requires library(ieugwasr)

# first makes sure CM GWAS is read in
CM_gwas <- read.delim("../../../../../summary-stats/Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

# format exposure data
CM_gwas$Phenotype <- "Maltreatment"
CM_gwas$N <- 185414

# rename columns
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N") 
names(CM_gwas) # correct

# reduce the size of the dataset to only suggestive SNPs (makes it faster to clump later)
CM_exp_dat <- CM_gwas %>% filter(P < 0.000005) 
# save list of suggestive SNPs for later
length(CM_exp_dat$SNP) # 1698

CM_exp_dat <- format_data(CM_exp_dat, type = "exposure",
                          snp_col = "SNP", beta_col = "BETA",
                          se_col = "SE", effect_allele_col = "A1",
                          other_allele_col = "A2", pval_col = "P",
                          samplesize_col = "N", min_pval = 1e-200,
                          #z_col = "Z", #info_col = "INFO",
                          chr_col = "CHR", pos_col = "BP") 


CM_exp_dat_clumped <- clump_data(CM_exp_dat) # clump

#---------------------------------#
#       A PATH ESTIMATES
#---------------------------------#
# define a path mr function
my.mr.a.path <- function(outcome.id, mediator){

  # extract snps from outcome
  out_dat <- extract_outcome_data(snps = CM_exp_dat_clumped$SNP, outcomes = outcome.id)
  
  dat <- harmonise_data(
    exposure_dat = CM_exp_dat_clumped, 
    outcome_dat = out_dat)
  
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
  ggsave(paste0(mediator, "-MR-plot-a-path.png"), combined_plot, dpi = 600, width = 15, height = 20, units = "in")
  
  return(list(res, intercept, Q))
}


### Cholesterol traits ###
# CM --> LDL
res_ldl_a_path <- my.mr.a.path(outcome.id = "ieu-a-300", mediator = "LDL-ieu-a-300") # 	GLGC
res_ldl_a_path2 <- my.mr.a.path(outcome.id = 'ebi-a-GCST90018961', mediator = 'LDL-ebi-a-GCST90018961')
# res_ldl_a_path3 <- my.mr.a.path(outcome.id = 'ieu-b-110')
# CM --> total cholesterol
res_tc_a_path <- my.mr.a.path(outcome.id = "ieu-a-301", mediator = "TC-ieu-a-301") # 	GLGC
res_tc_a_path2 <- my.mr.a.path(outcome.id = "ebi-a-GCST90018974", mediator = "TC-ebi-a-GCST90018974") 
#res_tc_a_path3 <- my.mr.a.path(outcome.id = "ebi-a-GCST90025953") 
# CM --> triglycerides
res_tg_a_path <- my.mr.a.path(outcome.id = "ieu-a-302", mediator = "TG-ieu-a-302") # 	GLGC
res_tg_a_path2 <- my.mr.a.path(outcome.id = "ieu-b-111", mediator = "TG-ieu-b-111") #  UKBB
# CM --> HDL
res_hdl_a_path <- my.mr.a.path(outcome.id = "ieu-a-299", mediator = "HDL-ieu-a-299") # 	GLGC
res_hdl_a_path2 <- my.mr.a.path(outcome.id = "ieu-b-109", mediator = "HDL-ieu-b-109")  # 	UKBB

### Glycaemic traits ###
# CM --> FI
res_fi_a_path <- my.mr.a.path(outcome.id = "ieu-b-116", mediator = "FI-ieu-b-116") # only 64,421 snps in this gwas
res_fi_a_path2 <- my.mr.a.path(outcome.id = "ebi-a-GCST90002238", mediator = "FI-ebi-a-GCST90002238")
# CM --> fasting glucose
res_fg_a_path <- my.mr.a.path(outcome.id = "ebi-a-GCST90002232", mediator = "FG-ebi-a-GCST90002232") 
res_fg_a_path2 <- my.mr.a.path(outcome.id = "ieu-b-113", mediator = "FG-ieu-b-113") #	Cortisol
# CM --> HbA1C
res_HbA1C_a_path <- my.mr.a.path(outcome.id = "ebi-a-GCST90014006", mediator = "HbA1C-ebi-a-GCST90014006") 
res_HbA1C_a_path2 <- my.mr.a.path(outcome.id = "ieu-b-4842", mediator = "HbA1C-ieu-b-4842") 

### Cortisol ###
# CM --> cortisol 
res_cortisol_a_path <- my.mr.a.path(outcome.id = "ieu-a-1012", mediator = "Cortisol-ieu-a-1012") # Plasma cortisol


#---------------------------------#
#       B PATH ESTIMATES
#---------------------------------#
# define b path mr function
my.mr.b.path <- function(exposure.data, mediator){
  out<- read_outcome_data(
    snps = exposure.data$SNP,
    filename = paste0("~/Documents/Projects/MR-mediation-CM-MM/summary-stats/PCM_multimorbidity_summary_stats_noUKBB_non-het.txt.gz"),
    sep = "\t", snp_col = "SNP",
    beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "MAF",
    pval_col = "P")
  out$outcome <- "MM"
  
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
  ggsave(paste0(mediator, "-MR-plot-b-path.png"), combined_plot, dpi = 600, width = 15, height = 20, units = "in")
  
  return(list(res, intercept, Q))
}

### Cholesterol traits ###
# LDL --> MM
ldl_exp_dat <- extract_instruments(outcomes = "ieu-a-300") # same as in MR from 2013 GLGC
res_ldl <- my.mr.b.path(ldl_exp_dat, "LDL-ieu-a-300") # 73 SNPs
ldl_exp_dat2 <- extract_instruments(outcomes = "ebi-a-GCST90018961") 
res_ldl2<- my.mr.b.path(ldl_exp_dat2, "LDL-ebi-a-GCST90018961")
# ldl_exp_dat3 <- extract_instruments(outcomes = 'ieu-b-110') 
# res_ldl3<- my.mr.b.path(ldl_exp_dat3)
# total cholesterol --> MM
tc_exp_dat <- extract_instruments(outcomes = "ieu-a-301") # same as in MR from 2013 GLGC
res_tc <- my.mr.b.path(tc_exp_dat, "TC-ieu-a-301") 
tc_exp_dat2 <- extract_instruments(outcomes =  "ebi-a-GCST90018974") 
res_tc2 <- my.mr.b.path(tc_exp_dat2, "TC-ebi-a-GCST90018974") 
# triglycerides --> MM
trig_exp_dat <- extract_instruments(outcomes =  "ieu-a-302") # same as in MR from 2013 GLGC
res_trig <- my.mr.b.path(trig_exp_dat, "TG-ieu-a-302") 
trig_exp_dat2 <- extract_instruments(outcomes =  "ieu-b-111") # ukbb
res_trig2 <- my.mr.b.path(trig_exp_dat2, "TG-ieu-b-111") 
# HDL --> MM
hdl_exp_dat <- extract_instruments(outcomes =  "ieu-a-299") # same as in MR from 2013 GLGC
res_hdl <- my.mr.b.path(hdl_exp_dat, "HDL-ieu-a-299") 
hdl_exp_dat2 <- extract_instruments(outcomes =  "ieu-b-109") # ukbb
res_hdl2 <- my.mr.b.path(hdl_exp_dat2, "HDL-ieu-b-109") 

### Glycaemic traits ###
# FI --> MM
fi_exp_dat <- extract_instruments(outcomes = "ieu-b-116") 
res_fi <- my.mr.b.path(fi_exp_dat, "FI-ieu-b-116") # 13 snps
fi_exp_dat2 <- extract_instruments(outcomes = "ebi-a-GCST90002238") 
res_fi2 <- my.mr.b.path(fi_exp_dat2, "FI-ebi-a-GCST90002238") 
# fasting glucose --> MM
fg_exp_dat <- extract_instruments(outcomes = "ebi-a-GCST90002232") 
res_fg <- my.mr.b.path(fg_exp_dat,  "FG-ebi-a-GCST90002232") 
fg_exp_dat2 <- extract_instruments(outcomes = "ieu-b-113") 
res_fg2 <- my.mr.b.path(fg_exp_dat2, "FG-ieu-b-113") 
# HbA1C --> MM
HbA1C_exp_dat <- extract_instruments(outcomes = "ebi-a-GCST90014006") 
res_HbA1C <- my.mr.b.path(HbA1C_exp_dat,"HbA1C-ebi-a-GCST90014006")
HbA1C_exp_dat2 <- extract_instruments(outcomes = "ieu-b-4842") 
res_HbA1C2 <- my.mr.b.path(HbA1C_exp_dat2, "HbA1C-ieu-b-4842") 

### Cortisol ###
cortisol_exp_dat <- extract_instruments(outcomes = "ieu-a-1012", p1 = 5e-06, p2 = 5e-06) # Plasma cortisol (only 1 snp, so setting 5e-06)
res_cortisol <- my.mr.b.path(cortisol_exp_dat, "Cortisol-ieu-a-1012") 


#---------------------------------#
#            SAVE OUTPUT
#---------------------------------#
# combine all output into one and save
## a paths
out_a_paths <- rbind(res_ldl_a_path[[1]], res_ldl_a_path2[[1]], res_tc_a_path[[1]], res_tc_a_path2[[1]], res_hdl_a_path[[1]], res_hdl_a_path2[[1]], res_tg_a_path[[1]], res_tg_a_path2[[1]], # cholesterol traits
                     res_fi_a_path[[1]], res_fi_a_path2[[1]], res_fg_a_path[[1]], res_fg_a_path2[[1]], res_HbA1C_a_path[[1]], res_HbA1C_a_path2[[1]], # glycaemic traits
                     res_cortisol_a_path[[1]]) # cortisol

out_intercepts_a <- rbind(res_ldl_a_path[[2]], res_ldl_a_path2[[2]], res_tc_a_path[[2]], res_tc_a_path2[[2]], res_hdl_a_path[[2]], res_hdl_a_path2[[2]], res_tg_a_path[[2]], res_tg_a_path2[[2]], # cholesterol traits
                     res_fi_a_path[[2]], res_fi_a_path2[[2]], res_fg_a_path[[2]], res_fg_a_path2[[2]], res_HbA1C_a_path[[2]], res_HbA1C_a_path2[[2]], # glycaemic traits
                     res_cortisol_a_path[[2]]) # cortisol

out_Q_a <- rbind(res_ldl_a_path[[3]], res_ldl_a_path2[[3]], res_tc_a_path[[3]], res_tc_a_path2[[3]], res_hdl_a_path[[3]], res_hdl_a_path2[[3]], res_tg_a_path[[3]], res_tg_a_path2[[3]], # cholesterol traits
                     res_fi_a_path[[3]], res_fi_a_path2[[3]], res_fg_a_path[[3]], res_fg_a_path2[[3]], res_HbA1C_a_path[[3]], res_HbA1C_a_path2[[3]], # glycaemic traits
                     res_cortisol_a_path[[3]]) # cortisol

# view significant effects (a path)
out_a_paths %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05) 
# HDL cholesterol || id:ieu-a-299, HDL cholesterol || id:ieu-b-109, Triglycerides || id:ieu-a-302, triglycerides || id:ieu-b-111,
# ebi-a-GCST90014006 Glycated haemoglobin HbA1c levels (UKB data field 30750) || id:ebi-a-GCST90014006, HbA1c || id:ieu-b-4842

## b paths
out_b_paths <- rbind(res_ldl[[1]], res_ldl2[[1]], res_tc[[1]], res_tc2[[1]], res_hdl[[1]], res_hdl2[[1]], res_trig[[1]], res_trig2[[1]], # cholesterol traits 
                     res_fi[[1]], res_fi2[[1]], res_fg[[1]], res_fg2[[1]], res_HbA1C[[1]], res_HbA1C2[[1]], # glycaemic traits 
                     res_cortisol[[1]]) # cortisol

out_intercepts_b <- rbind(res_ldl[[2]], res_ldl2[[2]], res_tc[[2]], res_tc2[[2]], res_hdl[[2]], res_hdl2[[2]], res_trig[[2]], res_trig2[[2]], # cholesterol traits 
                        res_fi[[2]], res_fi2[[2]], res_fg[[2]], res_fg2[[2]], res_HbA1C[[2]], res_HbA1C2[[2]],  # glycaemic traits 
                        res_cortisol[[2]]) # cortisol


out_Q_b <- rbind(res_ldl[[3]], res_ldl2[[3]], res_tc[[3]], res_tc2[[3]], res_hdl[[3]], res_hdl2[[3]], res_trig[[3]], res_trig2[[3]], # cholesterol traits 
               res_fi[[3]], res_fi2[[3]], res_fg[[3]], res_fg2[[3]], res_HbA1C[[3]], res_HbA1C2[[3]], # glycaemic traits 
               res_cortisol[[3]]) # cortisol

# view significant effects (b path)
out_b_paths %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05) # only fasting insulin (ebi-a-GCST90002238)
# LDL, TC, HDL, TG, FI, FG, and HbA1c (both data sources for all these)

# save a and b paths
out.list <- list(out_a_paths, out_b_paths, out_intercepts_a, out_intercepts_b, out_Q_a, out_Q_b)
names(out.list) <- c("a paths","b paths", "Egger intercepts a", "Egger intercepts b", "Q stats a", "Q stats b") 
openxlsx::write.xlsx(out.list, file = paste0("a-and-b-paths-unadjusted-uni-MR-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)


