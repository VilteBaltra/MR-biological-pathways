

#---------------------------------#
#       READ IN CM GWAS
#---------------------------------#
# make sure gwas sumstats are consistently named as my_mvmr function is expecting! 
# also note that sometimes mvmr function does not run due to traffic (?), but ends up running after a few tries 
# Because of this I edited it to use ld_clump(), instead of remote clumping function. Requires library(ieugwasr)

# first makes sure CM GWAS is read in
CM_gwas <- read.delim("../summary-stats/Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

# format exposure data
CM_gwas$Phenotype <- "Maltreatment"
CM_gwas$N <- 185414

# rename columns
# columns should be named c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Phenotype", "N") 
names(CM_gwas) # correct

#---------------------------------------------#
#  univariable MR for MD (total effect)
#---------------------------------------------#

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

CM_exp_dat <- clump_data(CM_exp_dat) # clump

out_dep <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz"), #  eaf FRQ_A_116209	or FRQ_U_314566 ?
  sep = "\t", snp_col = "SNP",
  beta_col = "logOR", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  pval_col = "P", 
  chr_col = "CHR", pos_col = "BP") 

dat <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = out_dep)

res <- mr(dat) 

# save output
write.csv(res, file = paste0("CM-MD-univariable-MR-", Sys.Date(), ".csv"), quote = F)

#---------------------------------------------#
#  univariable MR for CAD (total effect)
#---------------------------------------------#

out_cad <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/cad.add.160614.website.txt"), 
  sep = "\t", snp_col = "markername",
  beta_col = "beta", se_col = "se_dgc",
  effect_allele_col = "effect_allele", other_allele_col = "noneffect_allele",
  eaf_col = "effect_allele_freq",
  pval_col = "p_dgc", 
  chr_col = "chr", pos_col = "BP") 

dat <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = out_cad)

res <- mr(dat) 

# save output
write.csv(res, file = paste0("CM-CAD-univariable-MR-", Sys.Date(), ".csv"), quote = F)

#---------------------------------------------#
#  univariable MR for T2D (total effect)
#---------------------------------------------#

out_t2d <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/METAANALYSIS_DIAGRAM_SE1_mapped.txt"),
  sep = "\t", snp_col = "rsid",
  beta_col = "Effect", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  eaf_col = "effect_allele_freq",
  pval_col = "P")

dat <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = out_t2d)

res <- mr(dat) 

# save output
write.csv(res, file = paste0("CM-T2D-univariable-MR-", Sys.Date(), ".csv"), quote = F)

#---------------------------------------------#
#  univariable MR for MM (total effect)
#---------------------------------------------#

out_mm <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/PCM_multimorbidity_summary_stats_noUKBB_non-het.txt.gz"),
  sep = "\t", snp_col = "SNP",
  beta_col = "BETA", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  eaf_col = "MAF",
  pval_col = "P")

dat <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = out_mm)

res <- mr(dat) 

# save output
write.csv(res, file = paste0("CM-MM-univariable-MR-", Sys.Date(), ".csv"), quote = F)


