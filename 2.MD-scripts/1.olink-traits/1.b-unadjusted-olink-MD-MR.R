# MR SCRIPT - OLINK MARKERS (using depression as outcome)
# This script obtains unadjusted b path estimates (olink-->MD) that are later used to calculate the indirect effect (using the product of coefficients method)
# no need to obtain a path (CM-->olink) as already obtained as part of multimorbidity scripts
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
PATH = ""
setwd(PATH)
# define path to functions
FUNCTIONS = "" 
SUMSTATS=""

#---------------------------------#
###   SPECIFY OLINK GWASs NAMES ###
#---------------------------------#
# will loop through each of these GWAS for a and b paths
file_names <- c("VEGF.A-1", "TWEAK-1", "TNFSF14-1", "TGF.alpha-1", "TNFB-1", "TNFRSF9-1", "ST1A1-1", "PD.L1-1", "SIRT2-1", "SCF-1", "MCP.2-1", "MCP.3-1", "MCP.4-1", "MIP.1.alpha-1", "OSM-1",
                "CASP.8-1", "CCL11-1", "CCL19-1", "CCL20-1", "CCL23-1", "CCL25-1", "CCL4-1", "CD244-1", "CD40-1", "CD5-1", "CD6-1", "CSF.1-1", "CST5-1", "CXCL1-1", "CXCL10-1", "CXCL11-1", "CXCL9-1", "DNER-1", "EN.RAGE-1",
                "FGF.19-1", "FGF.23-1", "FGF.5-1", "GDNF-1", "IL.1.alpha-1", "IL.10RB-1", "IL.12B-1", "IL.15RA-1", "IL.17C-1", "IL.18-1", "IL.18R1-1", "IL.6-1", "IL.7-1", "IL.8-1",
                "LIF.R-1", "MCP.1-1", "MMP.1-1", "MMP.10-1", "OPG-1", "SLAMF1-1", "TRAIL-1", "uPA-1", "4E.BP1-1", "ADA-1", "CDCP1-1", "CX3CL1-1", "CXCL5-1", "CXCL6-1", "Flt3L-1", "HGF-1", "IL.10-1", "TRANCE-1", "FGF.21-1",
                "TSLP-1", "TNF-1", "STAMPB-1", "NRTN-1", "NT.3-1", "AXIN1-1", "Beta.NGF-1", "CCL28-1", "IFN.gamma-1", "IL.10RA-1", "IL.13-1", "IL.17A-1",
                "IL.2-1", "IL.20-1", "IL.20RA-1", "IL.22.RA1-1", "IL.24-1", "IL.2RB-1", "IL.33-1", "IL.4-1", "IL.5-1", "LIF-1", "LAP.TGF.beta.1-1", "ARTN-1")

#------------------------------------#
###  UNIVARIABLE MR LOOP (b path)  ###
#------------------------------------#

### READ IN ALL INFLAMMATORY MARKERS IN A LOOP ###

# create three empty dataframes for b path results
df <- data.frame()
egger_intercept_df_b <- data.frame()
Q_df_b <- data.frame()

for(trait in file_names){
  
  cat("Reading in", trait, "GWAS\n")
  gwas <- read.delim(paste0(SUMSTATS, trait, ".tbl.rsid.txt.gz"))
  
  # format exposure data
  gwas$Phenotype <- trait
  
  # rename columns
  # columns should be named "SNP", "CHR","POS", "BETA", "SE", "EA", "NEA", "P", "N", "Z", "INFO"
  colnames(gwas) <- c("CHR", "BP", "MarkerName", "A1", "A2", "eaf", "FreqSE", "MinFreq", "MaxFreq",  "BETA", "SE", "log.P.", "Direction", "HetISq", "HetChiSq", "HetDf", "logHetP", "N", "SNP", "Phenotype")
  
  # tranform log10 p value to standard p value
  gwas$P <- 10^(gwas$log.P.) 
  
  # reduce the size of the dataset to only suggestive SNPs (makes it faster to clump later)
  gwas <- gwas %>% filter(P < 0.000005) 
  # save list of suggestive SNPs for later
  length(gwas$SNP) 
  
  exp_dat <- format_data(gwas, type = "exposure",
                         snp_col = "SNP", beta_col = "BETA",
                         se_col = "SE", effect_allele_col = "A1",
                         other_allele_col = "A2", pval_col = "P",
                         samplesize_col = "N", min_pval = 1e-200,
                         eaf_col = 'eaf',
                         chr_col = "CHR", pos_col = "BP") 
  
  # if problems connecting to server, can use the chunk below (for local clumping version; requires 1kg.v3/EUR ref files to be downloaded)
  #exp_dat_clumped <- clump_data(exp_dat) # clump_kb = 10000, clump_r2 = 0.001
  
  clumped_snps <- ld_clump(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
                           plink_bin = genetics.binaRies::get_plink_binary(),
                           bfile = "/Users/vb506/Documents/Projects/MR-mediation-CM-MM/1kg.v3/EUR",
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p = 0.000005,
                           pop = "EUR" )
  # keep only clumped snps
  exp_dat_clumped <- subset(exp_dat, SNP %in% clumped_snps$rsid)
  
  out_dep <- read_outcome_data(
    snps = exp_dat_clumped$SNP,
    filename = paste0("../../../summary-stats/daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz"), 
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
  #write.csv(res, file = paste0(trait, "-MM-univariable-MR-unadjusted-", Sys.Date(), ".csv"), quote = F)
  # clean up
  
  # combine
  egger_intercept_df_b <- rbind(egger_intercept_df_b, egger_intercept_b)
  Q_df_b <- rbind(Q_df_b, Q_b)
  
  rm(res,gwas,dat) 
  gc()
}

# save all output
openxlsx::write.xlsx(df, file = paste0("olink-MD-uniMR-unadjusted-b-paths-", Sys.Date(), ".xlsx"),rowNames = F, overwrite = T)
openxlsx::write.xlsx(egger_intercept_df_b, file = paste0("olink-MD-egger-intercept-b-paths-", Sys.Date(), ".xlsx"),rowNames = F, overwrite = T)
openxlsx::write.xlsx(Q_df_b, file = paste0("olink-MD-Q-stat-b-paths-", Sys.Date(), ".xlsx"),rowNames = F, overwrite = T)

# check which IVW estimates are significant
df %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05)
# for MM was: TWEAK-1, SIRT2-1, CCL11-1, CSF.1-1, LIF.R-1, IL.33-1 
# for MD:  CASP.8-1, CXCL10-1,  FGF.5-1, IL.6-1,  TNF-1, IFN.gamma-1, IL.2-1, LIF-1


#---------------------------------------------#
#     UNIVARIABLE MR for CRP (b path)
#---------------------------------------------#
### CRP ### --> add to inflammation script
# CRP --> MM
# read in CRP gwas
crp_gwas <- data.table::fread(paste0("../../../summary-stats/CRP-Said-2022/GCST90029070_buildGRCh37.tsv.gz"), header = TRUE, sep = "\t")
# format exposure data
crp_gwas$Phenotype <- "CRP"
# make sure pvalue is numeric
crp_gwas$p_value <- as.numeric(crp_gwas$p_value)
# add sample size
crp_gwas$N <- 575531 # sample size taken from https://www.ebi.ac.uk/gwas/studies/GCST90029070

# rename columns
names(crp_gwas) <- c('SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA', 'SE', 'P', 'Phenotype', 'N')

# reduce the size of the dataset to only suggestive SNPs (makes it faster to clump later)
crp_exp_dat <- crp_gwas %>% filter(P < 0.00000005)
# save list of suggestive SNPs for later
length(crp_exp_dat$SNP) # 48919

crp_exp_dat <- format_data(crp_exp_dat, type = "exposure",
                           snp_col = "SNP", beta_col = "BETA",
                           se_col = "SE", effect_allele_col = "A1",
                           other_allele_col = "A2", pval_col = "P",
                           samplesize_col = "N", min_pval = 1e-200,
                           #z_col = "Z", #info_col = "INFO",
                           chr_col = "CHR", pos_col = "BP")

# replaced clump_data() with ld_clump() due to "Server code: 502; Server is possibly experiencing traffic, trying again..."
clumped_snps <- ld_clump(
  dplyr::tibble(rsid=crp_exp_dat$SNP, pval=crp_exp_dat$pval.exposure, id=crp_exp_dat$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/Users/vb506/Documents/Projects/MR-mediation-CM-MM/1kg.v3/EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  pop = "EUR" )

# keep only clumped snps
crp_exp_dat_clumped <- subset(crp_exp_dat, SNP %in% clumped_snps$rsid)

out_dep <- read_outcome_data(
  snps = crp_exp_dat_clumped$SNP,
  filename = paste0("../../../summary-stats/daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz"), 
  sep = "\t", snp_col = "SNP",
  beta_col = "logOR", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  pval_col = "P", eaf_col = 'FRQ_A_116209') #  eaf FRQ_A_116209	or FRQ_U_314566 ?
out_dep$outcome <- "Depression-noUKBB"

crp_dat <- harmonise_data(
  exposure_dat = crp_exp_dat_clumped,
  outcome_dat = out_dep)

# Horizontal pleiotropy
# The intercept term in MR Egger regression can be a useful indication of whether directional horizontal pleiotropy is driving the results of an MR analysis. 
egger_intercept_crp_b <- mr_pleiotropy_test(crp_dat)

# Heterogeneity statistic
Q_crp_b <- mr_heterogeneity(crp_dat)

# MR
res_crp_b <- mr(crp_dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

# combine output into one
res_crp_b1 <- merge(res_crp_b, Q_crp_b[, c('method', 'Q', 'Q_df', 'Q_pval')], by = 'method', all.x = TRUE)
# and two empty rows to cbind with output
egger_intercept_crp_b[2, ] <- NA
egger_intercept_crp_b[3, ] <- NA
res_crp_b2 <- cbind(res_crp_b1, egger_intercept_crp_b[, c('egger_intercept', 'se', 'pval')])

# save CRP a and b paths
openxlsx::write.xlsx(res_crp_b2, file = paste0("olink/CRP-MD-univariable-MR-unadjusted-b-", Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)
