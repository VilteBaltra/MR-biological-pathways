# MR SCRIPT - METABOLITE MARKERS (using depression as outcome)
# This script obtains unadjusted b path estimates (metabolites-->MD) that are later used to calculate the indirect effect (using the product of coefficients method)
# no need to obtain a path (CM-->metabolites) as already obtained as part of multimorbidity scripts
#---------------------------------#
#           SET UP
#---------------------------------#
# clean session
rm(list=ls())
gc()

# load libraries 
library(TwoSampleMR)
library(tidyverse)
library(ieugwasr)

# define path to working directory 
setwd(" ")

# define path to sumstats
SUMSTATS=" "

file_names_long <- dir("sumstats")
file_names <- sub("\\.h\\.tsv\\.gz$", "", file_names_long)

# first makes sure CM GWAS is read in
CM_gwas <- read.delim("Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

# format exposure data
CM_gwas$Phenotype <- "Maltreatment"
CM_gwas$N <- 185414

# rename columns
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N") 
names(CM_gwas) # correct

#------------------------------------#
###  UNIVARIABLE MR LOOP (b path)  ###
#------------------------------------#

### READ IN ALL METABOLITE MARKERS IN A LOOP ###

# create three empty dataframes for b path results
df <- data.frame()
egger_intercept_df_b <- data.frame()
Q_df_b <- data.frame()

set.seed(123)
for(trait in file_names){
  
  cat("Reading in", trait, "GWAS\n")
  gwas <- read.delim(paste0(SUMSTATS, trait, ".h.tsv.gz"))
  
  # format exposure data
  gwas$Phenotype <- trait
  
  # reduce the size of the dataset to only genome-wide significant SNPs (makes it faster to clump later)
  gwas <- gwas %>% filter(p_value < 0.00000005) 
  # save list of suggestive SNPs for later
  length(gwas$rsid) 
  
  exp_dat <- format_data(gwas, type = "exposure",
                         snp_col = "rsid", beta_col = "beta",
                         se_col = "standard_error", effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele", pval_col = "p_value",
                         #samplesize_col = "N", 
                         min_pval = 1e-200,
                         eaf_col = 'effect_allele_frequency',
                         chr_col = "chromosome", pos_col = "base_pair_location") 
  
  # if problems connecting to server, can use the chunk below (for local clumping version; requires 1kg.v3/EUR ref files to be downloaded)
  #exp_dat_clumped <- clump_data(exp_dat) # clump_kb = 10000, clump_r2 = 0.001
  
  clumped_snps2 <- ld_clump(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
                            plink_bin = genetics.binaRies::get_plink_binary(),
                            bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
                            clump_kb = 10000,
                            clump_r2 = 0.001,
                            clump_p = 0.00000005,
                            pop = "EUR" )
  # keep only clumped snps
  exp_dat_clumped <- subset(exp_dat, SNP %in% clumped_snps2$rsid)
  exp_dat_clumped$exposure <- trait
 
  out_dep <- read_outcome_data(
    snps = exp_dat_clumped$SNP,
    filename = paste0("daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz"), 
    sep = "\t", snp_col = "SNP",
    beta_col = "logOR", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", eaf_col = 'FRQ_A_116209') #  eaf FRQ_A_116209	or FRQ_U_314566 ?
  out_dep$outcome <- "Depression-noUKBB"
  
  dat <- harmonise_data(
    exposure_dat = exp_dat_clumped, 
    outcome_dat = out_dep)
  
  # Horizontal pleiotropy
  # The intercept term in MR Egger regression can be a useful indication of whether directional horizontal pleiotropy is driving the results of an MR analysis. 
  egger_intercept_b <- mr_pleiotropy_test(dat)
  
  # Heterogeneity statistic
  Q_b <- mr_heterogeneity(dat)
  
  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  
  df <- rbind(df, res)
  print(df)
  write.csv(res, file = paste0(trait, "-MD-univariable-MR-unadjusted-", Sys.Date(), ".csv"), quote = F)
  # clean up
  
  # combine
  egger_intercept_df_b <- rbind(egger_intercept_df_b, egger_intercept_b)
  Q_df_b <- rbind(Q_df_b, Q_b)
  
  rm(res,gwas,dat) 
  gc()
}

# save b paths
out.list <- list(df, egger_intercept_df_b, Q_df_b)
names(out.list) <- c("b paths", "Egger intercepts b", "Q stats b") 
openxlsx::write.xlsx(out.list, file = paste0("metabolites-b-paths-unadjusted-uni-MR-MD-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)

# check which IVW estimates are significant
df %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05)
