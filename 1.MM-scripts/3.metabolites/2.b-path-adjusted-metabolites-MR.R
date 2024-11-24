# MR SCRIPT - METABOLITES
# This script obtains adjusted b path that is later used to calculate the indirect effect (using the product of coefficients method)

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
# from https://rdrr.io/cran/MendelianRandomization/man/mr_mvegger.html
# https://rdrr.io/cran/MendelianRandomization/f/inst/doc/Vignette_MR.pdf 
library(MendelianRandomization) # for MVMR Egger
library(MVMR) # for conditional F stat and modified form of Cochran's Q statistic

# define path to working directory 
PATH = " "
setwd(PATH)
# define path to functions
FUNCTIONS = " " 
SUMSTATS=" "

# # source function
# source(paste0(FUNCTIONS, "my_mvmr_pval_olink.R")) # VB: performs MVMR analysis. Replaced clump_data() with local ld_clump() for faster performance (+ no server connection issues)
# # uses pval_threshold=5e-6 in mv_multiple(mvdat, pval_threshold=5e-6)
# source(paste0(FUNCTIONS, "my_mvmregger.R")) #  MVMR Egger function
# # define function to print order of traits (for MVMR Egger output)

# source functions
source(paste0(FUNCTIONS, "my_mvmr_pval.R")) # performs MVMR analysis. Replaced clump_data() with local ld_clump() for faster performance (+ no server connection issues)
source(paste0(FUNCTIONS, "my_mvmregger.R")) #  MVMR Egger function

obtain_exposure_names <- function(dat) {
  # print order of traits (relevant for MR Egger output)
  exposure1_name <- dat$expname %>% filter(id.exposure %in% colnames(dat$exposure_beta)[1])
  exposure2_name <- dat$expname %>% filter(id.exposure %in% colnames(dat$exposure_beta)[2])
  cat(paste0(" exposure 1 is ", exposure1_name[[2]], " (", colnames(dat$exposure_beta)[1], ")", "\n exposure 2 is ", exposure2_name[[2]], " (", colnames(dat$exposure_beta)[2], ")"))
  return(c(exposure1_name[[2]], exposure2_name[[2]]))
}

# MVMR steps taken from TwoSampleMR
# 1. Get instruments for each exposure
# 2. Combine these into a set of all instruments
# 3. Clump these to avoid the possibility that e.g. a variant for exposure 1 is in LD with a variant for exposure 2
# 4. Re-extract all the final clumped SNPs from (3) from all of the exposures
# 5. Harmonise them all to be on the same effect allele
# 6. Use the multivariable MR method against these harmonised data

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


#---------------------------------#
###   LOOP THROUGH OLINK GWASs  ###
#---------------------------------#
# this loop will obtain b path estimate / direct effect of CM-MM

# specify traits that passed early mediator selection steps
file_names <- dir(SUMSTATS)

### READ IN ALL INFLAMMATORY MARKERS IN A LOOP ###
setwd(PATH)

# create empty df for MVMR results 
df_res <- data.frame()

for(trait in file_names){
  
  cat("Reading in", trait, "GWAS\n")
  gwas <- read.delim(paste0(SUMSTATS, trait))
  
  # format exposure data
  gwas$Phenotype <- trait

  # rename columns
  # columns should be named "SNP", "CHR","POS", "BETA", "SE", "EA", "NEA", "P", "N", "Z", "INFO"
  names(gwas) <- c("CHR", "BP", "A1", "A2", "BETA", "SE", "eaf", "P", "variant_id", "SNP", "direction",         
                   "hetisq", "hetchisq", "hetdf", "hetpval", "hm_coordinate_conversion", "hm_code", "Phenotype")
  
  gwas_subset <- gwas[, c("CHR", "BP", "SNP", "P", "A1" , "A2", "BETA",  "SE", "Phenotype", "eaf")]
  head(gwas_subset)
  
  if(min(gwas_subset$P) > 0.000005){
    message(paste(trait, "minimum p-value is above 5e-6, so MR analysis will not run"))
  }
  
  # mvmr
  setwd(PATH)
  # run it uising mvmr_pval = 5e-6 
  mvdat <- my_mvmr_pval(exposure1.gwas = CM_gwas, exposure2.gwas = gwas_subset, # note, I am using EU ancestry as ref for clumping here (can change inside function)
                            exposure.names = c('cm', trait), exposure.number = 2, pval_1 = 5e-6, pval_2 = 5e-8, mvmr_pval = 5e-6,
                            return_mvdat=TRUE, rerun=FALSE)
  # pval_1 = pvalue threshold for exposure1.gwas hits
  # pval_2 = pvalue threshold for exposure2.gwas hits
  # mvmr_pval should be equivalent to the most lenient p-value in both exposures (used in mvmr analysis)
  # set return_mvdat=TRUE to obtain input data for mvmr Egger analysis
  # set rerun=TRUE if input files for MVMR my_mv_extract_exposures_local() function inside my_mvmr_pval() have already been obtained
  
  # run MVMR
  res_mvmr <- mv_multiple(mvdat, pval_threshold=5e-06) # fits all exposures together (recommended as default in TwoSampleMR)
  
  # generate OR
  res_mvmr <- generate_odds_ratios(res_mvmr$result)
  print(res_mvmr)
  
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
  colnames(mvmr_egger) <- paste0(colnames(mvmr_egger), "_egger")
  tmp <- merge(mvmr_egger, sres, by.x = "exposure_egger", by.y = "exposure")
  out <- cbind(tmp, pres)
  out <- merge(res_mvmr, out, by.x = "exposure", by.y = "exposure_egger")
  
  # combine 
  df_res <- rbind(df_res, out)
  print(df_res)
  
  # clear
  rm(gwas)
  rm(gwas_subset)
  gc() # releases memory
}

# save all
openxlsx::write.xlsx(df_res, file = paste0('metabolites-MM-MVMR-adjusted-b-paths-', Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)


