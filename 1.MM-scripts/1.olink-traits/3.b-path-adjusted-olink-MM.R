# MR SCRIPT - OLINK MARKERS
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

# source function
source(paste0(FUNCTIONS, "my_mvmr_pval_olink.R")) # VB: performs MVMR analysis. Replaced clump_data() with local ld_clump() for faster performance (+ no server connection issues)
# uses pval_threshold=5e-6 in mv_multiple(mvdat, pval_threshold=5e-6)
source(paste0(FUNCTIONS, "my_mvmregger.R")) #  MVMR Egger function
# define function to print order of traits (for MVMR Egger output)

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
# Get a list of all files in the working directory
# file_names <- c("VEGF.A-1", "TWEAK-1", "TNFSF14-1", "TGF.alpha-1", "TNFB-1", "TNFRSF9-1", "ST1A1-1", "PD.L1-1", "SIRT2-1", "SCF-1", "MCP.2-1", "MCP.3-1", "MCP.4-1", "MIP.1.alpha-1", "OSM-1",
#                 "CASP.8-1", "CCL11-1", "CCL19-1", "CCL20-1", "CCL23-1", "CCL25-1", "CCL4-1", "CD244-1", "CD40-1", "CD5-1", "CD6-1", "CSF.1-1", "CST5-1", "CXCL1-1", "CXCL10-1", "CXCL11-1", "CXCL9-1", "DNER-1", "EN.RAGE-1",
#                 "FGF.19-1", "FGF.23-1", "FGF.5-1", "GDNF-1", "IL.1.alpha-1", "IL.10RB-1", "IL.12B-1", "IL.15RA-1", "IL.17C-1", "IL.18-1", "IL.18R1-1", "IL.6-1", "IL.7-1", "IL.8-1",
#                 "LIF.R-1", "MCP.1-1", "MMP.1-1", "MMP.10-1", "OPG-1", "SLAMF1-1", "TRAIL-1", "uPA-1", "4E.BP1-1", "ADA-1", "CDCP1-1", "CX3CL1-1", "CXCL5-1", "CXCL6-1", "Flt3L-1", "HGF-1", "IL.10-1", "TRANCE-1", "FGF.21-1",
#                 "TSLP-1", "TNF-1", "STAMPB-1", "NRTN-1", "NT.3-1", "AXIN1-1", "Beta.NGF-1", "CCL28-1", "IFN.gamma-1", "IL.10RA-1", "IL.13-1", "IL.17A-1", 
#                 "IL.2-1", "IL.20-1", "IL.20RA-1", "IL.22.RA1-1", "IL.24-1", "IL.2RB-1", "IL.33-1", "IL.4-1", "IL.5-1", "LIF-1", "LAP.TGF.beta.1-1", "ARTN-1") 

# specify traits that passed early mediator selection steps
file_names <- c("MIP.1.alpha-1", "HGF-1", "FGF.21-1", "CDCP1-1", "SIRT2-1", "LIF.R-1", "CSF.1-1")


### READ IN ALL INFLAMMATORY MARKERS IN A LOOP ###
setwd(PATH)

# create empty df for MVMR results 
df_res <- data.frame()
# create empty df for MVMR Egger results 
df_egger <- data.frame()

for(trait in file_names){
  
  cat("Reading in", trait, "GWAS\n")
  gwas <- read.delim(paste0(SUMSTATS, trait, ".tbl.rsid.txt.gz"))
  
  # format exposure data
  gwas$Phenotype <- trait

  # rename columns
  # columns should be named "SNP", "CHR","POS", "BETA", "SE", "EA", "NEA", "P", "N", "Z", "INFO"
  names(gwas) <- c("CHR", "BP", "MarkerName", "A1", "A2", "eaf", "FreqSE", "MinFreq", "MaxFreq",  "BETA", "SE", "log.P.", "Direction", "HetISq", "HetChiSq", "HetDf", "logHetP", "N", "SNP", "Phenotype")
  
  # tranform log10 p value to standard p value
  gwas$P <- 10^(gwas$log.P.) 
  
  gwas_subset <- gwas[, c("CHR", "BP", "SNP", "P", "A1" , "A2", "BETA",  "SE", "Phenotype", "N", "eaf")]
  head(gwas_subset)
  
  if(min(gwas_subset$P) > 0.000005){
    message(paste(trait, "minimum p-value is above 5e-6, so MR analysis will not run"))
  }
  
  # mvmr
  setwd(PATH)
  output.list <- my_mvmr_pval_olink(exposure1.gwas = CM_gwas, exposure2.gwas = gwas_subset,
                              exposure.names = c('cm', trait), exposure.number = 2, 
                              pval_1 = 5e-6, pval_2 = 5e-6, mvmr_pval = 5e-6,
                              return_mvdat=TRUE, rerun=FALSE)
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
  #write.csv(out, file = paste0('results_mvmr_egger_', trait, '_MM.csv'))
  # combine 
  df_egger <- rbind(df_egger, out)
  print(df_egger)
  # combine 
  df_res <- rbind(df_res, res_mvmr)
  print(df_res)
  # clear
  rm(gwas)
  rm(gwas_subset)
  gc() # releases memory
}

# save all
openxlsx::write.xlsx(list('mvmr' = df_res, 'mvmr_egger' = df_egger), file =  paste0('olink-MM-MVMR-Egger-b-paths-', Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)


