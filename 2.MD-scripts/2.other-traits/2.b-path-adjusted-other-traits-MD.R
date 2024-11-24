# MR SCRIPT - OTHER TRAITS (outcome = MD)
# This script obtains adjusted b paths that are used to calculate the indirect effect (using the product of coefficients method)

#---------------------------------#
#           SET UP
#---------------------------------#
# clean session
rm(list=ls())
gc()

# load libraries 
library(TwoSampleMR) # for two-sample MR
library(tidyverse)
library(ieugwasr)

# for vcf files I'll use a vcf package 
# from https://mrcieu.github.io/gwasvcf/articles/guide.html
#remotes::install_github("mrcieu/gwasvcf")
library(gwasvcf)
#set_bcftools('/path/to/bcftools')
#remotes::install_github('mrcieu/genetics.binaRies')
#set_plink()
#set_bcftools()

suppressWarnings(suppressPackageStartupMessages({
  library(gwasvcf)
  library(VariantAnnotation)
  library(dplyr)
  library(magrittr)
}))

# from https://rdrr.io/cran/MendelianRandomization/man/mr_mvegger.html
# https://rdrr.io/cran/MendelianRandomization/f/inst/doc/Vignette_MR.pdf 
library(MendelianRandomization) # for MVMR Egger
library(MVMR) # for conditional F stat and modified form of Cochran's Q statistic

# define path to working directory 
PATH = " "
setwd(PATH)
# define path to functions
FUNCTIONS = " " 

# source functions
source(paste0(FUNCTIONS, "my_mvmr_pval_olink_md.R")) # performs MVMR analysis. Replaced clump_data() with local ld_clump() for faster performance (+ no server connection issues)
source(paste0(FUNCTIONS, "my_mvmregger.R")) #  MVMR Egger function

# define function to print order of traits (for MVMR Egger output)
obtain_exposure_names <- function(dat) {
  # print order of traits (relevant for MR Egger output)
  exposure1_name <- dat$expname %>% filter(id.exposure %in% colnames(dat$exposure_beta)[1])
  exposure2_name <- dat$expname %>% filter(id.exposure %in% colnames(dat$exposure_beta)[2])
  cat(paste0(" exposure 1 is ", exposure1_name[[2]], " (", colnames(dat$exposure_beta)[1], ")", "\n exposure 2 is ", exposure2_name[[2]], " (", colnames(dat$exposure_beta)[2], ")"))
  return(c(exposure1_name[[2]], exposure2_name[[2]]))
}

#---------------------------------#
#       READ IN CM GWAS
#---------------------------------#
# make sure gwas sumstats are consistently named as my_mvmr function is expecting! 
# also note that sometimes mvmr function does not run due to traffic (?), but ends up running after a few tries 
# Because of this I edited it to use ld_clump(), instead of remote clumping function. Requires library(ieugwasr)

# first makes sure CM GWAS is read in
CM_gwas <- read.delim("../../../../summary-stats/Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

# format exposure data
CM_gwas$Phenotype <- "Maltreatment"
CM_gwas$N <- 185414

# rename columns
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N") 
names(CM_gwas) # correct


#---------------------------------#
#       READ IN CRP GWAS
#---------------------------------#
# read in CRP gwas 
crp_gwas <- data.table::fread(paste0("../../../../summary-stats/CRP-Said-2022/GCST90029070_buildGRCh37.tsv.gz"), header = TRUE, sep = "\t")
# format exposure data
crp_gwas$Phenotype <- "CRP"
# make sure pvalue is numeric
crp_gwas$p_value <- as.numeric(crp_gwas$p_value)
# add sample size
crp_gwas$N <- 575531 # sample size taken from https://www.ebi.ac.uk/gwas/studies/GCST90029070 

# rename columns
names(crp_gwas) <- c('SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA', 'SE', 'P', 'Phenotype', 'N')

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it uising mvmr_pval = 5e-6 
mvdat_crp <- my_mvmr_pval_olink_md(exposure1.gwas = CM_gwas, exposure2.gwas = crp_gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                          exposure.names = c('cm', 'CRP'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                          return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
# set return_mvdat=TRUE to obtain input data for mvmr Egger analysis
# set rerun=TRUE if input files for MVMR my_mv_extract_exposures_local() function inside my_mvmr_pval() have already been obtained


#---------------------------------#
#     RUN MVMR (with UKBB)
#---------------------------------#
library(data.table)
CM_gwas <- fread("Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

# format exposure data
CM_gwas$Phenotype <- "Maltreatment"
CM_gwas$N <- 185414
CM_gwas <- as.data.frame(CM_gwas)
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N") 

file_names <- "logTG"

# create three empty dataframes for b path results
df_res <- data.frame()
sensitivity_df <- data.frame()

for(trait in file_names){
  
  cat("Reading in", paste0(trait, "withUKBB"), "GWAS\n")
  gwas <- fread(paste0("../lipids/",trait, "_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"))
  gwas <- as.data.frame(gwas)
  gwas$Phenotype <- paste0(trait, "withUKBB")
  # columns should contain c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N")
  colnames(gwas) <- c("SNP", "CHR", "BP", "A2", "A1", "N", "N_studies", "EAF", "BETA", "SE", "pvalue_neg_log10",  "pvalue", "pvalue_neg_log10_GC", "P", "Phenotype") 
  
  # run it using mvmr_pval = 5e-6 
  mvdat <- my_mvmr_pval_olink_md(exposure1.gwas = CM_gwas, exposure2.gwas = gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                        exposure.names = c('cm', paste0(trait, "withUKBB")), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                        return_mvdat=TRUE, rerun=FALSE)
  # pval_1 = pvalue threshold for exposure1.gwas hits
  # pval_2 = pvalue threshold for exposure2.gwas hits
  # mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
  # set return_mvdat=TRUE to obtain input data for mvmr Egger analysis
  # set rerun=TRUE if input files for MVMR my_mv_extract_exposures_local() function inside my_mvmr_pval() have already been obtained
  
  # run MVMR
  # res_mvmr <- mv_multiple(mvdat[[2]], pval_threshold=5e-6) # same as below
  res_mvmr <- mvdat[[1]]
  df_res <- rbind(df_res, as.data.frame(res_mvmr))
  
  # run MVMR Egger 
  mvmr_egger <- my_mvmregger(mvdat[[2]])
  mvmr_egger$exposure <- obtain_exposure_names(mvdat[[2]])
  
  # calculate conditional F statistics for instrument strength
  my.F  <- format_mvmr(BXGs = mvdat[[2]]$exposure_beta,
                       BYG = mvdat[[2]]$outcome_beta,
                       seBXGs = mvdat[[2]]$exposure_se,
                       seBYG = mvdat[[2]]$outcome_se,
                       RSID = rownames(mvdat[[2]]$exposure_beta))
  sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
  sres <- as.data.frame(t(sres))
  sres$exposure <- obtain_exposure_names(mvdat[[2]])
  
  # calculates modified form of Cochran's Q statistic 
  pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)
  
  # save MVMR Egger with Fstat (per trait) and Q stat (overall)
  out <- merge(mvmr_egger, sres, by = "exposure")
  out <- cbind(out, pres)
  sensitivity_df <- rbind(sensitivity_df, out)
  print(df_res)
  print(sensitivity_df)
  #write.csv(out, file = paste0('cm-', paste0(trait, "withUKBB"), '-md/results_mvmr_egger_CM.', paste0(trait, "withUKBB"), '.MD.csv'))
  rm(gwas, mvdat, out)
  gc()
}


#---------------------------------#
#     RUN MVMR (without UKBB) 
#---------------------------------#

for(trait in file_names){
  
  cat("Reading in", trait, "GWAS\n")
  gwas <- fread(paste0("../lipids/without_UKB_",trait, "_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"))
  gwas <- as.data.frame(gwas)
  gwas$Phenotype <- trait
  # columns should contain c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N")
  colnames(gwas) <- c("SNP", "CHR", "BP", "A2", "A1", "N", "N_studies", "EAF", "BETA", "SE", "pvalue_neg_log10",  "pvalue", "pvalue_neg_log10_GC", "P", "Phenotype") 
  
  # run it using mvmr_pval = 5e-6 
  mvdat <- my_mvmr_pval_olink_md(exposure1.gwas = CM_gwas, exposure2.gwas = gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                                 exposure.names = c('cm', trait), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                                 return_mvdat=TRUE, rerun=FALSE)
  # pval_1 = pvalue threshold for exposure1.gwas hits
  # pval_2 = pvalue threshold for exposure2.gwas hits
  # mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
  # set return_mvdat=TRUE to obtain input data for mvmr Egger analysis
  # set rerun=TRUE if input files for MVMR my_mv_extract_exposures_local() function inside my_mvmr_pval() have already been obtained
  
  # run MVMR
  # res_mvmr <- mv_multiple(mvdat[[2]], pval_threshold=5e-6) # same as below
  res_mvmr <- mvdat[[1]]
  df_res <- rbind(df_res, as.data.frame(res_mvmr))
  
  # run MVMR Egger 
  mvmr_egger <- my_mvmregger(mvdat[[2]])
  mvmr_egger$exposure <- obtain_exposure_names(mvdat[[2]])
  
  # calculate conditional F statistics for instrument strength
  my.F  <- format_mvmr(BXGs = mvdat[[2]]$exposure_beta,
                       BYG = mvdat[[2]]$outcome_beta,
                       seBXGs = mvdat[[2]]$exposure_se,
                       seBYG = mvdat[[2]]$outcome_se,
                       RSID = rownames(mvdat[[2]]$exposure_beta))
  sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
  sres <- as.data.frame(t(sres))
  sres$exposure <- obtain_exposure_names(mvdat[[2]])
  
  # calculates modified form of Cochran's Q statistic 
  pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)
  
  # save MVMR Egger with Fstat (per trait) and Q stat (overall)
  out <- merge(mvmr_egger, sres, by = "exposure")
  out <- cbind(out, pres)
  sensitivity_df <- rbind(sensitivity_df, out)
  print(df_res)
  print(sensitivity_df)
  #write.csv(out, file = paste0('cm-', trait), '-md/results_mvmr_egger_CM.', paste0(trait, "withUKBB"), '.MD.csv'))
  rm(gwas, mvdat, out)
  gc()

}

# save output
openxlsx::write.xlsx(list(df_res=df_res, sensitivity_df=sensitivity_df), file = paste0("lipids-b-paths-adjusted-with-and-withoutUKBB-DEPRESSION-", Sys.Date(), ".xlsx"), 
                     rowNames = T, overwrite = T)


