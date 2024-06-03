# MR SCRIPT - OLINK MARKERS
# This script obtains beta1 (a path) and beta2 (b path) paths that are used to calculate the indirect effect (using the product of coefficients method)

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
remotes::install_github('mrcieu/genetics.binaRies')
set_plink()
set_bcftools()

suppressWarnings(suppressPackageStartupMessages({
  library(gwasvcf)
  library(VariantAnnotation)
  library(dplyr)
  library(magrittr)
}))
set_bcftools()

# from https://rdrr.io/cran/MendelianRandomization/man/mr_mvegger.html
# https://rdrr.io/cran/MendelianRandomization/f/inst/doc/Vignette_MR.pdf 
library(MendelianRandomization) # for MVMR Egger
library(MVMR) # for conditional F stat and modified form of Cochran's Q statistic

# define path to working directory 
PATH = "/Users/vb506/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/"
setwd(PATH)
# define path to functions
FUNCTIONS = "/Users/vb506/Documents/Projects/MR-mediation-CM-MM/scripts/TIDY-WORKFLOW/source/" 

# source functions
source(paste0(FUNCTIONS, "my_mvmr_pval.R")) # performs MVMR analysis. Replaced clump_data() with local ld_clump() for faster performance (+ no server connection issues)
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
CM_gwas <- read.delim("../../summary-stats/Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

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
crp_gwas <- data.table::fread(paste0("../../summary-stats/CRP-Said-2022/GCST90029070_buildGRCh37.tsv.gz"), header = TRUE, sep = "\t")
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
mvdat_crp <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = crp_gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
        exposure.names = c('cm', 'CRP'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
        return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
# set return_mvdat=TRUE to obtain input data for mvmr Egger analysis
# set rerun=TRUE if input files for MVMR my_mv_extract_exposures_local() function inside my_mvmr_pval() have already been obtained

# run MVMR Egger 
mvmr_egger <- my_mvmregger(mvdat_crp)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_crp)

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_crp$exposure_beta,
                     BYG = mvdat_crp$outcome_beta,
                     seBXGs = mvdat_crp$exposure_se,
                     seBYG = mvdat_crp$outcome_se,
                     RSID = rownames(mvdat_crp$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_crp)

# calculates modified form of Cochran's Q statistic 
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-CRP-mm/results_mvmr_egger_CM.CRP.MM.csv')

#---------------------------------#
#    READ IN TG GWAS (ieu-b-111)
#---------------------------------#
# read in TG (ieu-b-111.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf <- readVcf("../../summary-stats/sumstats-2024/ieu-b-111.vcf.gz")

# view key info
vcf
# dim: 12321875 1 
#  List of length 9: ES, SE, LP, AF, SS, EZ, SI, NC, ID

# explore it
header(vcf)
samples(header(vcf)) # "ieu-b-111"
rowRanges(vcf)

# Create a data frame:
TG_vcf <- vcf_to_granges(vcf) %>% dplyr::as_tibble()
# ID is rsid 
# ES = effect size
# SE = standard error 
# AF = alternative allele freq
# SS = sample size 
# LP A      Float  -log10 p-value for effect estimate   

# add sample size
TG_vcf$N <- 441016 # sample size taken from https://gwas.mrcieu.ac.uk/datasets/ieu-b-111/

# select key columns
TG_vcf_sub <- TG_vcf[, c('ID', 'REF', 'ALT', 'ES', 'SE', 'AF', 'LP', 'id', 'N')]
head(TG_vcf_sub)

# rename columns
names(TG_vcf_sub) <- c('SNP', 'A2', 'A1', 'BETA', 'SE', 'EAF', 'LP', 'Phenotype', 'N')
head(TG_vcf_sub)

# convert to normal p-value
TG_vcf_sub$P <- 10^(-TG_vcf_sub$LP)

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it uising mvmr_pval = 5e-6 
mvdat_tg <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = TG_vcf_sub, # note, I am using EU ancestry as ref for clumping here (can change inside function)
             exposure.names = c('cm', 'TG-ieu-b-111'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
             return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
# set return_mvdat=TRUE to obtain input data for mvmr Egger analysis
# set rerun=TRUE if input files for MVMR my_mv_extract_exposures_local() function inside my_mvmr_pval() have already been obtained

# run MVMR Egger 
mvmr_egger <- my_mvmregger(mvdat_tg)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_tg)

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_tg$exposure_beta,
                     BYG = mvdat_tg$outcome_beta,
                     seBXGs = mvdat_tg$exposure_se,
                     seBYG = mvdat_tg$outcome_se,
                     RSID = rownames(mvdat_tg$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_tg)

# calculates modified form of Cochran's Q statistic 
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-TG-ieu-b-111-mm/results_mvmr_egger_CM.TG-ieu-b-111.MM.csv')


#---------------------------------#
#    READ IN TG GWAS (ieu-a-302)
#---------------------------------#
# read in TG (ieu-a-302.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf2 <- readVcf("../../summary-stats/sumstats-2024/ieu-a-302.vcf.gz")
# dim: 2415963 1 
#  List of length 9:  ES, SE, LP, AF, SS, EZ, SI, NC, ID

# ensure correct ID
samples(header(vcf2)) # "ieu-a-302"

# Create a data frame:
TG_vcf2 <- vcf_to_granges(vcf2) %>% dplyr::as_tibble()
# ID is rsid 
# ES = effect size
# SE = standard error 
# AF = alternative allele freq
# SS = sample size 
# LP A      Float  -log10 p-value for effect estimate   

# add sample size
TG_vcf2$N <- 177861 # sample size taken from https://gwas.mrcieu.ac.uk/datasets/ieu-a-302/

# select key columns
TG_vcf_sub2 <- TG_vcf2[, c('ID', 'REF', 'ALT', 'ES', 'SE', 'AF', 'LP', 'id', 'N')]
head(TG_vcf_sub2)

# rename columns
names(TG_vcf_sub2) <- c('SNP', 'A2', 'A1', 'BETA', 'SE', 'EAF', 'LP', 'Phenotype', 'N')
head(TG_vcf_sub2)

# convert to normal p-value
TG_vcf_sub2$P <- 10^(-TG_vcf_sub2$LP)

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it uising mvmr_pval = 5e-6 
mvdat_tg2 <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = TG_vcf_sub2, # note, I am using EU ancestry as ref for clumping here (can change inside function)
             exposure.names = c('cm', 'TG-ieu-a-302'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
             return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)

# run MVMR Egger 
mvmr_egger <- my_mvmregger(mvdat_tg2)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_tg2)

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_tg2$exposure_beta,
                     BYG = mvdat_tg2$outcome_beta,
                     seBXGs = mvdat_tg2$exposure_se,
                     seBYG = mvdat_tg2$outcome_se,
                     RSID = rownames(mvdat_tg2$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_tg2)

# calculates modified form of Cochran's Q statistic 
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-TG-ieu-a-302-mm/results_mvmr_egger_CM.TG-ieu-a-302.MM.csv')


#---------------------------------#
#   READ IN HbA1C GWAS (ieu-b-4842)
#---------------------------------#
# read in TG (ieu-b-4842.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf_hba1c <- readVcf("../../summary-stats/sumstats-2024/ieu-b-4842.vcf.gz")
vcf_hba1c
# dim:  9696819 1 
#  List of length 9:  ES, SE, LP, AF, SS, EZ, SI, NC, ID

# ensure correct ID
samples(header(vcf_hba1c)) # "ieu-b-4842"

# Create a data frame:
hba1c_gwas <- vcf_to_granges(vcf_hba1c) %>% dplyr::as_tibble()
# ID is rsid 
# ES = effect size
# SE = standard error 
# AF = alternative allele freq
# SS = sample size 
# LP A      Float  -log10 p-value for effect estimate   

# didn't add sample size as already available in SS
# hba1c_gwas$N <- 45734 # sample size taken from https://gwas.mrcieu.ac.uk/datasets/ieu-b-4842/

# select key columns
hba1c_gwas_sub <- hba1c_gwas[, c('ID', 'REF', 'ALT', 'ES', 'SE', 'AF', 'LP', 'id', 'SS')]
head(hba1c_gwas_sub)

# rename columns
names(hba1c_gwas_sub) <- c('SNP', 'A2', 'A1', 'BETA', 'SE', 'EAF', 'LP', 'Phenotype', 'N')
head(hba1c_gwas_sub)

# convert to normal p-value
hba1c_gwas_sub$P <- 10^(-hba1c_gwas_sub$LP)

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it uising mvmr_pval = 5e-6 
mvdat_hba1c <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = hba1c_gwas_sub, # note, I am using EU ancestry as ref for clumping here (can change inside function)
             exposure.names = c('cm', 'HbA1c-ieu-b-4842'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
             return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)

# run MVMR Egger 
mvmr_egger <- my_mvmregger(mvdat_hba1c)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_hba1c)

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_hba1c$exposure_beta,
                     BYG = mvdat_hba1c$outcome_beta,
                     seBXGs = mvdat_hba1c$exposure_se,
                     seBYG = mvdat_hba1c$outcome_se,
                     RSID = rownames(mvdat_hba1c$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_hba1c)

# calculates modified form of Cochran's Q statistic 
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-HbA1c-ieu-b-4842-mm/results_mvmr_egger_CM.HbA1c-ieu-b-4842.MM.csv')

#---------------------------------#
#   READ IN HbA1c GWAS (ebi-a-GCST90014006)
#---------------------------------#
# read in HbA1c (GCST90014006_buildGRCh38.tsv.rsid.txt.gz) gwas 
# added rsid to sumstats as was not present
# first makes sure CM GWAS is read in
HbA1c2_gwas <- read.delim("../../summary-stats/sumstats-2024/GCST90014006_buildGRCh38.tsv.rsid.txt.gz")

# format exposure data
HbA1c2_gwas$Phenotype <- "HbA1c-ebi-a-GCST90014006"
HbA1c2_gwas$N <- 389889 # sample sizes taken from https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90014006/

# rename columns
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N") 
names(HbA1c2_gwas) <- c("CHR", "BP", "P", "A1", "A2", "BETA", "SE", "SNP", "Phenotype", "N") 
head(HbA1c2_gwas)

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it using mvmr_pval = 5e-6 
mvdat_HbA1c2 <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = HbA1c2_gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                             exposure.names = c('cm', 'HbA1c-ebi-a-GCST90014006'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                             return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)

# run MVMR Egger 
mvmr_egger <- my_mvmregger(mvdat_HbA1c2)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_HbA1c2) 

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_HbA1c2$exposure_beta,
                     BYG = mvdat_HbA1c2$outcome_beta,
                     seBXGs = mvdat_HbA1c2$exposure_se,
                     seBYG = mvdat_HbA1c2$outcome_se,
                     RSID = rownames(mvdat_HbA1c2$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_HbA1c2) 

# calculates modified form of Cochran's Q statistic 
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-HbA1c-ebi-a-GCST90014006-mm/results_mvmr_egger_CM.HbA1c-ebi-a-GCST90014006.MM.csv')

#---------------------------------#
#   READ IN HDL GWAS (id:ieu-b-109)
#---------------------------------#
# read in HDL (ieu-b-109.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf_hdl <- readVcf("../../summary-stats/sumstats-2024/ieu-b-109.vcf.gz")
vcf_hdl
# dim:  12321875 1 
#  List of length 9:  ES, SE, LP, AF, SS, EZ, SI, NC, ID

# ensure correct ID
samples(header(vcf_hdl)) # "ieu-b-109"

# Create a data frame:
hdl_gwas <- vcf_to_granges(vcf_hdl) %>% dplyr::as_tibble()
# ID is rsid 
# ES = effect size
# SE = standard error 
# AF = alternative allele freq
# SS = sample size 
# LP A      Float  -log10 p-value for effect estimate   

# add sample size
hdl_gwas$N <- 403943 # sample size taken from https://gwas.mrcieu.ac.uk/datasets/ieu-b-109/

# select key columns
hdl_gwas_sub <- hdl_gwas[, c('ID', 'REF', 'ALT', 'ES', 'SE', 'AF', 'LP', 'id', 'N')]
head(hdl_gwas_sub)

# rename columns
names(hdl_gwas_sub) <- c('SNP', 'A2', 'A1', 'BETA', 'SE', 'EAF', 'LP', 'Phenotype', 'N')
head(hdl_gwas_sub)

# convert to normal p-value
hdl_gwas_sub$P <- 10^(-hdl_gwas_sub$LP)

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it using mvmr_pval = 5e-6 
mvdat_hdl <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = hdl_gwas_sub, # note, I am using EU ancestry as ref for clumping here (can change inside function)
             exposure.names = c('cm', 'HDL-ieu-b-109'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
             return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)

# run MVMR Egger 
mvmr_egger <- my_mvmregger(mvdat_hdl)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_hdl)

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_hdl$exposure_beta,
                     BYG = mvdat_hdl$outcome_beta,
                     seBXGs = mvdat_hdl$exposure_se,
                     seBYG = mvdat_hdl$outcome_se,
                     RSID = rownames(mvdat_hdl$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_hdl)

# calculates modified form of Cochran's Q statistic 
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-HDL-ieu-b-109-mm/results_mvmr_egger_CM.HDL-ieu-b-109.MM.csv')

#---------------------------------#
#   READ IN HDL GWAS (id:ieu-a-299)
#---------------------------------#
# read in HDL (ieu-a-299.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf_hdl2 <- readVcf("../../summary-stats/sumstats-2024/ieu-a-299.vcf.gz")
vcf_hdl2

# ensure correct ID
samples(header(vcf_hdl2)) # "ieu-a-299"

# Create a data frame:
hdl2_gwas <- vcf_to_granges(vcf_hdl2) %>% dplyr::as_tibble()
# ID is rsid 
# ES = effect size
# SE = standard error 
# AF = alternative allele freq
# SS = sample size 
# LP A      Float  -log10 p-value for effect estimate   

# add sample size
hdl2_gwas$N <-  hdl2_gwas$SS

# select key columns
hdl2_gwas_sub <- hdl2_gwas[, c('ID', 'REF', 'ALT', 'ES', 'SE', 'AF', 'LP', 'id', 'N')]
head(hdl2_gwas_sub)

# rename columns
names(hdl2_gwas_sub) <- c('SNP', 'A2', 'A1', 'BETA', 'SE', 'EAF', 'LP', 'Phenotype', 'N')
head(hdl2_gwas_sub)

# convert to normal p-value
hdl2_gwas_sub$P <- 10^(-hdl2_gwas_sub$LP)

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it using mvmr_pval = 5e-6
mvdat_hdl2 <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = hdl2_gwas_sub, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                          exposure.names = c('cm', 'HDL-ieu-a-299'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                          return_mvdat=TRUE, rerun=FALSE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)

# run MVMR Egger
mvmr_egger <- my_mvmregger(mvdat_hdl2)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_hdl2)

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_hdl2$exposure_beta,
                     BYG = mvdat_hdl2$outcome_beta,
                     seBXGs = mvdat_hdl2$exposure_se,
                     seBYG = mvdat_hdl2$outcome_se,
                     RSID = rownames(mvdat_hdl2$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_hdl2)

# calculates modified form of Cochran's Q statistic
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-HDL-ieu-a-299-mm/results_mvmr_egger_CM.HDL-ieu-a-299.MM.csv')

#---------------------------------#
#   READ IN LDL GWAS (ebi-a-GCST90018961)
#---------------------------------#
# read in LDL (GCST90018961_buildGRCh37.tsv.rsid.txt.gz) gwas 
# added rsid to sumstats as was not present
# first makes sure CM GWAS is read in
ldl2_gwas <- read.delim("../../summary-stats/sumstats-2024/GCST90018961_buildGRCh37.tsv.rsid.txt.gz")

# format exposure data
ldl2_gwas$Phenotype <- "LDL-ebi-a-GCST90018961"
ldl2_gwas$N <- 343621 # sample sizes taken from https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90018961/ 

# rename columns
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N")
names(ldl2_gwas) <- c("CHR", "BP", "A1", "A2", "EAF", "BETA", "SE", "P", "variant_id", "SNP", "Phenotype", "N")
head(ldl2_gwas)

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on LDL (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for LDL

# run it using mvmr_pval = 5e-6 
mvdat_ldl2 <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = ldl2_gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                             exposure.names = c('cm', 'LDL-ebi-a-GCST90018961'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                             return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)

# run MVMR Egger 
mvmr_egger <- my_mvmregger(mvdat_ldl2)
mvmr_egger$exposure <- obtain_exposure_names(mvdat_ldl2)

# calculate conditional F statistics for instrument strength
my.F  <- format_mvmr(BXGs = mvdat_ldl2$exposure_beta,
                     BYG = mvdat_ldl2$outcome_beta,
                     seBXGs = mvdat_ldl2$exposure_se,
                     seBYG = mvdat_ldl2$outcome_se,
                     RSID = rownames(mvdat_ldl2$exposure_beta))
sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
sres <- as.data.frame(t(sres))
sres$exposure <- obtain_exposure_names(mvdat_ldl2)

# calculates modified form of Cochran's Q statistic 
pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)

# save MVMR Egger with Fstat (per trait) and Q stat (overall)
out <- merge(mvmr_egger, sres, by = "exposure")
out <- cbind(out, pres)
write.csv(out, file = 'cm-LDL-ebi-a-GCST90018961-mm/results_mvmr_egger_CM.LDL-ebi-a-GCST90018961.MM.csv')


#---------------------------------#
#        COMBINE MVMR
#---------------------------------#
# combine MVMR output into a single document
df <- data.frame() # initiate empty df 

for(trait in c('TG-ieu-a-302',  'TG-ieu-b-111', 'HDL-ieu-b-109', 'HDL-ieu-a-299', 'LDL-ebi-a-GCST90018961', 'HbA1c-ieu-b-4842', 'HbA1c-ebi-a-GCST90014006',  'CRP')){
  
  cat("Obtaining estimates for", trait, "\n")
  
  BETA2 = paste0("cm-", trait, "-mm/cm-", trait, "-mm", "-MVMR-results-") # not full name as dates differ for some files
  
  # Wrap the file reading and product.coef calculation in tryCatch
  result <- tryCatch({
    b_df <- read.csv(dir(dirname(BETA2), full.names=T, pattern=paste("^", basename(BETA2), sep='')), row.names=1)
    b_df$mediator <- trait # add name of mediator
    
    df <- rbind(df, b_df)
    rm(b_df)
    
    # Print a success message
    cat("Estimates obtained successfully for", trait, "\n")
  }, error = function(err) {
    cat("Error occurred for", trait, ":", conditionMessage(err), "\n")
    # Handle the error as needed
  })
}


#---------------------------------#
#         SAVE ALL OUTPUT
#---------------------------------#

# save adjusted b paths output
openxlsx::write.xlsx(df, file = paste0("b-paths-adjusted-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)

# extract a path estimates and combine with adjusted b path estimates obtained earlier
a_df <- openxlsx::read.xlsx("a-and-b-paths-unadjusted-uni-MR-2024-01-31.xlsx", sheet = 1) # unadjusted a paths
# a_df "ebi-a-GCST90014006" (HBA1C), "ebi-a-GCST90018961" (LDL) renamed to "HbA1c-ebi-a-CST90014006" (HBA1C), "LDL-ebi-a-GCST90018961" to match adjusted b data (in a path of excel file)
# and also added CRP to a and b paths of excel file
b_df <- df # adjusted b paths
#b_df_unadjusted <- openxlsx::read.xlsx("a-and-b-paths-unadjusted-uni-MR-2024-01-31.xlsx", sheet = 2) # unadjusted b paths

# select matching a path ids 
a_path_ids <- c('CRP', 'ieu-b-4842', 'ieu-b-111', 'ieu-a-302', 'ieu-b-109', 'ieu-a-299', 'HbA1c-ebi-a-CST90014006', 'LDL-ebi-a-GCST90018961')

# filter to matching traits
a <- (a_df %>% filter(method == "Inverse variance weighted") %>% filter(id.outcome %in% a_path_ids))
b <- (b_df %>% filter(exposure != "Maltreatment")) %>% filter(exposure %in% a_path_ids)

# ensure order of traits for a and b path are identical
a_b_merged <- merge(a, b, by.x = 'id.outcome', by.y = 'exposure')

# combine output into same dataframe
a_and_b_adjusted <- cbind(a = a_b_merged$b.x, b = a_b_merged$b.y, 
                          a.se=a_b_merged$se.x, b.se=a_b_merged$se.y, 
                          a.pval=a_b_merged$pval.x, b.pval=a_b_merged$pval.y)
rownames(a_and_b_adjusted) <- a_b_merged$mediator # add mediator name as row name


# save output
openxlsx::write.xlsx(as.data.frame(a_and_b_adjusted), file = paste0("selected-a-and-b-paths-adjusted-", Sys.Date(), ".xlsx"), 
                     rowNames = T, overwrite = T)



