
# This script obtains univariable (unadjusted) MR estimates for cholesterol traits, glycaemic traits, cortisol and CRP
# trait IDs taken from https://gwas.mrcieu.ac.uk/ 

library(tidyverse)
library(TwoSampleMR) 
library(ieugwasr)
library(cowplot)

# path
setwd("")

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
ldl_exp_dat2 <- extract_instruments(outcomes = "ebi-a-GCST90018961") 
res_ldl2<- my.mr.b.path(ldl_exp_dat2, mediator="ldl-ebi-a-GCST90018961") 
# total cholesterol --> MD
tc_exp_dat2 <- extract_instruments(outcomes =  "ebi-a-GCST90018974") 
res_tc2 <- my.mr.b.path(tc_exp_dat2, mediator="tc-ebi-a-GCST90018974")
# triglycerides --> MD
trig_exp_dat2 <- extract_instruments(outcomes =  "ieu-b-111") # ukbb
res_trig2 <- my.mr.b.path(trig_exp_dat2, mediator="tg-ieu-b-111") 
# HDL --> MD
hdl_exp_dat2 <- extract_instruments(outcomes = "ieu-b-109") # ukbb
res_hdl2 <- my.mr.b.path(hdl_exp_dat2, mediator = "hdl-ieu-b-109") 

### Glycaemic traits ### 
# FI --> MD
fi_exp_dat2 <- extract_instruments(outcomes = "ebi-a-GCST90002238")
res_fi2 <- my.mr.b.path(fi_exp_dat2, mediator = "fi-ebi-a-GCST90002238" )
# fasting glucose --> MD
fg_exp_dat <- extract_instruments(outcomes = "ebi-a-GCST90002232") 
res_fg <- my.mr.b.path(fg_exp_dat, mediator = "fg-ebi-a-GCST90002232")

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
out_b_paths <- rbind(res_ldl2[[1]], res_tc2[[1]], res_hdl2[[1]], res_trig2[[1]], # cholesterol traits 
                     res_fi2[[1]], res_fg[[1]], res_HbA1C[[1]], res_HbA1C2[[1]], # glycaemic traits 
                     res_cortisol[[1]]) # cortisol
      
# view significant effects
out_b_paths %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05) # only fasting insulin (ebi-a-GCST90002238)

out_intercepts <- rbind(res_ldl2[[2]], res_tc2[[2]], res_hdl2[[2]], res_trig2[[2]], # cholesterol traits 
                     res_fi2[[2]], res_fg[[2]], res_HbA1C[[2]], res_HbA1C2[[2]],  # glycaemic traits 
                     res_cortisol[[2]]) # cortisol
     

out_Q <- rbind(res_ldl2[[3]], res_tc2[[3]], res_hdl2[[3]], res_trig2[[3]], # cholesterol traits 
                        res_fi2[[3]], res_fg[[3]], res_HbA1C[[3]], res_HbA1C2[[3]], # glycaemic traits 
                        res_cortisol[[3]]) # cortisol

# save all
out.list <- list(out_b_paths, out_intercepts, out_Q)
names(out.list) <- c("unadjusted b paths", "Egger intercepts", "Q statistics") 
openxlsx::write.xlsx(out.list, file = paste0("b-paths-unadjusted-Egger-Q-stat-uni-MR-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)



#------------------------------------#
###  UNIVARIABLE MR LOOP (b path)  ### 
#------------------------------------#
# ### in response to co-author suggestion adding new, larger GWAS for lipids from Nightingale GWAS ###

### READ IN ALL LIPIDS IN A LOOP ###

# create three empty dataframes for b path results
df_b <- data.frame()
egger_intercept_df_b <- data.frame()
Q_df_b <- data.frame()

library(data.table)

file_names <- c("HDL", "LDL", "TC", "nonHDL", "logTG")

# version with UKBB for MD
# using GC corrected p-value (as for path a)
for(trait in file_names){

  cat("Reading in", trait, "GWAS\n")
  gwas <- fread(paste0(trait, "_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"))

  # format exposure data
  gwas$Phenotype <- paste0(trait, "withUKBB")

  # rename columns
  # columns should be named "SNP", "CHR","POS", "BETA", "SE", "EA", "NEA", "P", "N", "Z", "INFO"
  # [1] "rsID"                "CHROM"               "POS_b37"
  # [4] "REF"                 "ALT"                 "N"
  # [7] "N_studies"           "POOLED_ALT_AF"       "EFFECT_SIZE"
  # [10] "SE"                  "pvalue_neg_log10"    "pvalue"
  # [13] "pvalue_neg_log10_GC" "pvalue_GC"           "Phenotype"

  colnames(gwas) <- c("SNP", "CHR", "BP", "NEA", "EA", "N", "N_studies", "EAF", "BETA", "SE", "pvalue_neg_log10",  "pvalue", "pvalue_neg_log10_GC", "P", "Phenotype")

  # reduce the size of the dataset to only suggestive SNPs (makes it faster to clump later)
  gwas <- gwas %>% filter(P < 0.000005)
  # save list of suggestive SNPs for later
  length(gwas$SNP)
  gwas <- data.frame(gwas)
  exp_dat <- format_data(gwas, type = "exposure",
                         snp_col = "SNP", beta_col = "BETA",
                         se_col = "SE", effect_allele_col = "EA",
                         other_allele_col = "NEA", pval_col = "P",
                         samplesize_col = "N", min_pval = 1e-200,
                         eaf_col = 'EAF',
                         chr_col = "CHR", pos_col = "BP")

  # if problems connecting to server, can use the chunk below (for local clumping version; requires 1kg.v3/EUR ref files to be downloaded)
  #exp_dat_clumped <- clump_data(exp_dat) # clump_kb = 10000, clump_r2 = 0.001

  clumped_snps <- ld_clump(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
                           plink_bin = genetics.binaRies::get_plink_binary(),
                           bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p = 0.00000005,
                           pop = "EUR" )
  # keep only clumped snps
  exp_dat_clumped <- subset(exp_dat, SNP %in% clumped_snps$rsid)

  out_mm <- read_outcome_data(
    snps = exp_dat_clumped$SNP,
    filename = paste0("../metabolites/daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz"), # no eaf in sumstats
    sep = "\t", snp_col = "SNP",
    beta_col = "logOR", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "FRQ_A_116209", #  eaf FRQ_A_116209	or FRQ_U_314566 ?
    pval_col = "P")
  out_mm$outcome <- "Depression"

  dat <- harmonise_data(
    exposure_dat = exp_dat_clumped,
    outcome_dat = out_mm)

  # Horizontal pleiotropy
  # The intercept term in MR Egger regression can be a useful indication of whether directional horizontal pleiotropy is driving the results of an MR analysis.
  egger_intercept_b <- mr_pleiotropy_test(dat)

  # Heterogeneity statistic
  Q_b <- mr_heterogeneity(dat)

  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

  df_b <- rbind(df_b, res)
  print(df_b)

  # combine
  egger_intercept_df_b <- rbind(egger_intercept_df_b, egger_intercept_b)
  Q_df_b <- rbind(Q_df_b, Q_b)

  rm(res,gwas,dat)
  gc()
}

# save all output
openxlsx::write.xlsx(list(res_b_path = df_b, egger_intercept_df_b = egger_intercept_df_b, Q_df_b = Q_df_b), file = paste0("lipids-large-withUKBB-DEPRESSION-uniMR-b-paths-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)

# check which IVW estimates are significant 
df_b %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05)


# version without UKBB
# using GC corrected p-value (as for path a)
file_names <- c("HDL", "LDL", "TC", "nonHDL", "logTG")
for(trait in file_names){

  cat("Reading in", trait, "GWAS\n")
  gwas <- fread(paste0("without_UKB_", trait, "_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"))

  # format exposure data
  gwas$Phenotype <- trait

  # rename columns
  # columns should be named "SNP", "CHR","POS", "BETA", "SE", "EA", "NEA", "P", "N", "Z", "INFO"
  # [1] "rsID"                "CHROM"               "POS_b37"
  # [4] "REF"                 "ALT"                 "N"
  # [7] "N_studies"           "POOLED_ALT_AF"       "EFFECT_SIZE"
  # [10] "SE"                  "pvalue_neg_log10"    "pvalue"
  # [13] "pvalue_neg_log10_GC" "pvalue_GC"           "Phenotype"

  colnames(gwas) <- c("SNP", "CHR", "BP", "NEA", "EA", "N", "N_studies", "EAF", "BETA", "SE", "pvalue_neg_log10",  "pvalue", "pvalue_neg_log10_GC", "P", "Phenotype")

  # reduce the size of the dataset to only suggestive SNPs (makes it faster to clump later)
  gwas <- gwas %>% filter(P < 0.000005)
  # save list of suggestive SNPs for later
  length(gwas$SNP)
  gwas <- data.frame(gwas)
  exp_dat <- format_data(gwas, type = "exposure",
                         snp_col = "SNP", beta_col = "BETA",
                         se_col = "SE", effect_allele_col = "EA",
                         other_allele_col = "NEA", pval_col = "P",
                         samplesize_col = "N", min_pval = 1e-200,
                         eaf_col = 'EAF',
                         chr_col = "CHR", pos_col = "BP")

  # if problems connecting to server, can use the chunk below (for local clumping version; requires 1kg.v3/EUR ref files to be downloaded)
  #exp_dat_clumped <- clump_data(exp_dat) # clump_kb = 10000, clump_r2 = 0.001

  clumped_snps <- ld_clump(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
                           plink_bin = genetics.binaRies::get_plink_binary(),
                           bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p = 0.00000005,
                           pop = "EUR" )
  # keep only clumped snps
  exp_dat_clumped <- subset(exp_dat, SNP %in% clumped_snps$rsid)
  
  out_mm <- read_outcome_data(
    snps = exp_dat_clumped$SNP,
    filename = paste0("../metabolites/daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz"), # no eaf in sumstats
    sep = "\t", snp_col = "SNP",
    beta_col = "logOR", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "FRQ_A_116209", #  eaf FRQ_A_116209	or FRQ_U_314566 ?
    pval_col = "P")
  out_mm$outcome <- "Depression"

  dat <- harmonise_data(
    exposure_dat = exp_dat_clumped,
    outcome_dat = out_mm)

  # Horizontal pleiotropy
  # The intercept term in MR Egger regression can be a useful indication of whether directional horizontal pleiotropy is driving the results of an MR analysis.
  egger_intercept_b <- mr_pleiotropy_test(dat)

  # Heterogeneity statistic
  Q_b <- mr_heterogeneity(dat)

  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

  df_b <- rbind(df_b, res)
  print(df_b)

  # combine
  egger_intercept_df_b <- rbind(egger_intercept_df_b, egger_intercept_b)
  Q_df_b <- rbind(Q_df_b, Q_b)

  rm(res,gwas,dat)
  gc()
}


# save all output
openxlsx::write.xlsx(list(res_b_path = df_b, egger_intercept_df_b = egger_intercept_df_b, Q_df_b = Q_df_b), file = paste0("lipids-large-with-and-withoutUKBB-DEPRESSION-uniMR-b-paths-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)


