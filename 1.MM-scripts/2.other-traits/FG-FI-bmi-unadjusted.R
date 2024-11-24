.libPaths(c("/campaign/VB-FM5HPC-001/Vilte/pkgs", "/apps/hc/software/2021a/R/4.1.0-foss-2021b/lib64/R/library"))

# load libraries 
library(TwoSampleMR)
library(tidyverse)
library(ieugwasr)

setwd("/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites")

#SUMSTATS="/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/sumstats/"

# first makes sure CM GWAS is read in
CM_gwas <- read.delim("Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

# format exposure data
CM_gwas$Phenotype <- "Maltreatment"
CM_gwas$N <- 185414

# rename columns
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N") 
names(CM_gwas) # correct

#---------------------------------------------#
#      UNIVARIABLE MR LOOP (a path)
#---------------------------------------------#
# this loop will obtain a path estimate for each olink protein 
# univariable MR

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

clumped_snps <- ld_clump(
  dplyr::tibble(rsid=CM_exp_dat$SNP, pval=CM_exp_dat$pval.exposure, id=CM_exp_dat$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.000005,
  pop = "EUR" )

# keep only clumped snps
CM_exp_dat2 <- subset(CM_exp_dat, SNP %in% clumped_snps$rsid)

# create three empty dataframes where results will be stored
df_a_paths <- data.frame()
egger_intercept_df <- data.frame()
Q_df <- data.frame()

# outcome names
file_names_long <- c("FG_combined_1000G_density_formatted_21-03-29.txt.gz", "FI_combined_1000G_density_formatted_21-03-29.txt.gz")
file_names <- sub("\\_combined_1000G_density_formatted_21-03-29.txt\\.gz$", "", file_names_long)

# run loop
set.seed(123)
for(trait in file_names){
  
  cat(paste("Running univariable MR between CM and", trait, "\n"))
  
  out <- read_outcome_data(
    snps = CM_exp_dat2$SNP,
    filename = paste0(trait, "_combined_1000G_density_formatted_21-03-29.txt.gz"),
    sep = "\t", snp_col = "rsid",
    beta_col = "beta", se_col = "se",
    effect_allele_col = "a1", other_allele_col = "a2",
    pval_col = "p-value", 
    log_pval = FALSE,  # log_pval = TRUE when the pval is -log10(P). Our pval is log10(P), not negative
    #my_log_pval = TRUE, # added my_log_pval command to account for non-negative log(P) values
    #id_col = trait,  
    #chr_col = "chromosome", pos_col = "base_pair_location",
    #samplesize_col = ,
    eaf_col = "maf") 
  out$outcome <- trait
  
  dat <- harmonise_data(
    exposure_dat = CM_exp_dat2, 
    outcome_dat = out)
  
  # Horizontal pleiotropy
  # The intercept term in MR Egger regression can be a useful indication of whether directional horizontal pleiotropy is driving the results of an MR analysis. 
  egger_intercept <- mr_pleiotropy_test(dat)
  
  # Heterogeneity statistic
  Q <- mr_heterogeneity(dat)
  
  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  
  df_a_paths <- rbind(res, df_a_paths)
  print(df_a_paths)
  #write.csv(res, file = paste0("CM-", trait, "-univariable-MR-", Sys.Date(), ".csv"), quote = F)
  
  # combine
  egger_intercept_df <- rbind(egger_intercept_df, egger_intercept)
  Q_df <- rbind(Q_df, Q)
  
  # clear
  rm(out, dat, res)
  gc() # releases memory
}


# save a paths
out.list <- list(df_a_paths, egger_intercept_df, Q_df)
names(out.list) <- c("a paths", "Egger intercepts a", "Q stats a") 
openxlsx::write.xlsx(out.list, file = paste0("FI-FG-a-paths-bmi-unadjusted-uni-MR-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)


#------------------------------------#
###  UNIVARIABLE MR LOOP (b path)  ###
#------------------------------------#

# create three empty dataframes for b path results
df <- data.frame()
egger_intercept_df_b <- data.frame()
Q_df_b <- data.frame()

set.seed(123)
for(trait in file_names){
  
  cat("Reading in", trait, "GWAS\n")
  gwas <- read.delim(paste0(trait, "_combined_1000G_density_formatted_21-03-29.txt.gz"))
  
  # format exposure data
  gwas$Phenotype <- trait
  
  # reduce the size of the dataset to only genome-wide significant SNPs (makes it faster to clump later)
  gwas <- gwas %>% filter(p.value < 0.00000005) 
  # save list of suggestive SNPs for later
  length(gwas$rsid) 
  
  exp_dat <- format_data(gwas, type = "exposure",
                         snp_col = "rsid", beta_col = "beta",
                         se_col = "se", effect_allele_col = "a1",
                         other_allele_col = "a2", pval_col = "p.value",
                         #samplesize_col = "N", 
                         min_pval = 1e-200,
                         eaf_col = 'maf') 
  
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
  
  out_mm <- read_outcome_data(
    snps = exp_dat_clumped$SNP,
    filename = paste0("PCM_multimorbidity_summary_stats_noUKBB_non-het.txt.gz"), # no eaf in sumstats
    sep = "\t", snp_col = "SNP",
    beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", eaf_col = 'MAF')
  out_mm$outcome <- "multimorbidity"
  
  dat <- harmonise_data(
    exposure_dat = exp_dat_clumped, 
    outcome_dat = out_mm)
  
  # Horizontal pleiotropy
  # The intercept term in MR Egger regression can be a useful indication of whether directional horizontal pleiotropy is driving the results of an MR analysis. 
  egger_intercept_b <- mr_pleiotropy_test(dat)
  
  # Heterogeneity statistic
  Q_b <- mr_heterogeneity(dat)
  
  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  
  df <- rbind(df, res)
  print(df)
  #write.csv(res, file = paste0(trait, "-MM-univariable-MR-unadjusted-", Sys.Date(), ".csv"), quote = F)
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
openxlsx::write.xlsx(out.list, file = paste0("FI-FG-b-paths-bmi-unadjusted-uni-MR-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)

