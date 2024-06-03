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
PATH = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/DEPRESSION/other-traits/"
setwd(PATH)
# define path to functions
FUNCTIONS = "/Users/vb506/Documents/Projects/MR-mediation-CM-MM/scripts/TIDY-WORKFLOW/source/" 

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
#    READ IN TG GWAS (ieu-b-111)
#---------------------------------#
# read in TG (ieu-b-111.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf <- readVcf("../../../../summary-stats/sumstats-2024/ieu-b-111.vcf.gz")

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

# clean up 
rm(vcf, TG_vcf, crp_gwas)
gc()
#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it uising mvmr_pval = 5e-6 
mvdat_tg <- my_mvmr_pval_olink_md(exposure1.gwas = CM_gwas, exposure2.gwas = TG_vcf_sub, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                         exposure.names = c('cm', 'TG-ieu-b-111'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                         return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
# set return_mvdat=TRUE to obtain input data for mvmr Egger analysis
# set rerun=TRUE if input files for MVMR my_mv_extract_exposures_local() function inside my_mvmr_pval() have already been obtained

#---------------------------------#
#    READ IN TG GWAS (ieu-a-302)
#---------------------------------#
# read in TG (ieu-a-302.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf2 <- readVcf("../../../../summary-stats/sumstats-2024/ieu-a-302.vcf.gz")
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

# clean up
rm(vcf2, TG_vcf2)
gc()

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it uising mvmr_pval = 5e-6 
mvdat_tg2 <- my_mvmr_pval_olink_md(exposure1.gwas = CM_gwas, exposure2.gwas = TG_vcf_sub2, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                          exposure.names = c('cm', 'TG-ieu-a-302'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                          return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)


#---------------------------------#
#   READ IN HDL GWAS (id:ieu-b-109)
#---------------------------------#
# read in HDL (ieu-b-109.vcf.gz) gwas 

# read in vcf gwas sumstats
vcf_hdl <- readVcf("../../../../summary-stats/sumstats-2024/ieu-b-109.vcf.gz")
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

# clean up
rm(vcf_hdl, hdl_gwas)
gc()

#---------------------------------#
#            RUN MVMR
#---------------------------------#
# obtains direct effect of CM on CRP (to be used as adjusted b path in prod.coef mediation)
# also obtains direct effect of CM on MM, controlling for CRP

# run it using mvmr_pval = 5e-6 
mvdat_hdl <- my_mvmr_pval_olink_md(exposure1.gwas = CM_gwas, exposure2.gwas = hdl_gwas_sub, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                          exposure.names = c('cm', 'HDL-ieu-b-109'), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                          return_mvdat=TRUE, rerun=TRUE)
# pval_1 = pvalue threshold for exposure1.gwas hits
# pval_2 = pvalue threshold for exposure2.gwas hits
# mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)


#-----------------------------------------------#
#     COMBINE MVMR with Sensitivity analyses
#-----------------------------------------------#
# create empty df for MVMR results 
df_res <- data.frame()
# create empty df for MVMR Egger results 
df_egger <- data.frame()

# specify function to add sensitivity analyses to mvmr results
# saves output in df_res and df_egger
all.output <- function(output.list, trait){
  
  # select res_mvmr output
  res_mvmr <- output.list[[1]]
  res_mvmr$mediator <- trait # add name of mediator
  
  # select mvdat output
  mvdat <- output.list[[2]]
  
  # run MVMR Egger 
  mvmr_egger <- my_mvmregger(mvdat)
  mvmr_egger$exposure <- obtain_exposure_names(mvdat) 
  
  # calculate conditional F statistics for instrument strength
  my.F  <- format_mvmr(BXGs = mvdat$exposure_beta,
                       BYG = mvdat$outcome_beta,
                       seBXGs = mvdat$exposure_se,
                       seBYG = mvdat$outcome_se,
                       RSID = rownames(mvdat$exposure_beta))
  sres <- strength_mvmr(r_input = my.F, gencov = 0) # Fixing covariance at 0
  sres <- as.data.frame(t(sres))
  sres$exposure <- obtain_exposure_names(mvdat)
  
  # calculates modified form of Cochran's Q statistic 
  pres <- pleiotropy_mvmr(r_input = my.F, gencov = 0)
  
  # save MVMR Egger with Fstat (per trait) and Q stat (overall)
  out <- merge(mvmr_egger, sres, by = "exposure")
  out <- cbind(out, pres)
  out$mediator <- trait # add name of mediator
  print(out)
  #return
  return(list(res_mvmr, out))

}

# our output: mvdat_crp, mvdat_tg, mvdat_tg2, mvdat_hba1c, mvdat_hdl, mvdat_hdl2, mvdat_fi2
results_all1 <- all.output(output.list = mvdat_crp, trait = "CRP")
results_all2 <- all.output(output.list = mvdat_tg, trait = "TG-ieu-b-111")
results_all3 <- all.output(output.list = mvdat_tg2, trait = "TG-ieu-a-302")
results_all4 <- all.output(output.list = mvdat_hdl, trait = "HDL-ieu-b-109")

# combine all output
res_combined <- rbind(results_all1[[1]], results_all2[[1]], results_all3[[1]], results_all4[[1]])
sensitivity_combined <- rbind(results_all1[[2]], results_all2[[2]], results_all3[[2]], results_all4[[2]])

# save all
openxlsx::write.xlsx(res_combined, file = paste0('other-traits-MD-MVMR-adjusted-b-paths-', Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)
openxlsx::write.xlsx(sensitivity_combined, file =  paste0('other-traits-MD-MVMR-Egger-b-paths-', Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)

# view
res_combined %>% filter(exposure != "Maltreatment") %>% filter(pval < 0.05) # none
