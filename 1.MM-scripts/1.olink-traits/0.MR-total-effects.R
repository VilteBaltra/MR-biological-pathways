
setwd("~/Documents/Projects/MR-mediation-CM-MM/output/")

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
write.csv(res, file = paste0("CM-MD-univariable-MR-5e-8-", Sys.Date(), ".csv"), quote = F)

#------------------------------------------------#
#  univariable MR for MD (total effect) with UKBB
#------------------------------------------------#
out_dep <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/PGC_UKB_23andMe_depression_genome-wide.txt.gz"), #  eaf FRQ_A_116209	or FRQ_U_314566 ?
  sep = " ", snp_col = "SNP",
  beta_col = "LOGOR", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  pval_col = "P",
  eaf_col = "FREQ1") 

dat <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = out_dep)

res <- mr(dat) 

# save output
write.csv(res, file = paste0("CM-MDwithUKBB-univariable-MR-", Sys.Date(), ".csv"), quote = F)

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
write.csv(res, file = paste0("CM-CAD-univariable-MR-5e-8-", Sys.Date(), ".csv"), quote = F)

#-------------------------------------------------#
#  univariable MR for CAD (total effect) withUKBB
#-------------------------------------------------#

out_cad <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz"), 
  sep = "\t", snp_col = "snptestid",
  beta_col = "logOR", se_col = "se_gc",
  effect_allele_col = "effect_allele", other_allele_col = "noneffect_allele",
  eaf_col = "effect_allele_freq",
  pval_col = "p-value_gc",
  chr_col = "chr", pos_col = "bp_hg19") 
out_cad$outcome <- "CADwithUKBB"

dat <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = out_cad)

res <- mr(dat) 

# save output
write.csv(res, file = paste0("CM-CADwithUKBB-univariable-MR-", Sys.Date(), ".csv"), quote = F)

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
write.csv(res, file = paste0("CM-T2D-univariable-MR-5e-8-", Sys.Date(), ".csv"), quote = F)

#------------------------------------------------#
#  univariable MR for T2D (total effect) withUKBB
#------------------------------------------------#

out_t2d <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt.gz"),
  sep = "\t", snp_col = "RSID",
  beta_col = "Beta", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  eaf_col = "EAF",
  pval_col = "P")

dat <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = out_t2d)

res <- mr(dat) 

# save output
write.csv(res, file = paste0("CM-T2DwithUKBB-univariable-MR-", Sys.Date(), ".csv"), quote = F)


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
write.csv(res, file = paste0("CM-MM-univariable-MR-5e-8-", Sys.Date(), ".csv"), quote = F)


#-----------------------------------------------#
#  univariable MR for MM (total effect) withUKBB
#-----------------------------------------------#
# this uses MM outcome with UKBB
out_mm <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = paste0("../summary-stats/PCM_multimorbidity_summary_stats_withUKBB_non-het.txt.gz"),
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
write.csv(res, file = paste0("CM-MMwithUKBB-univariable-MR-", Sys.Date(), ".csv"), quote = F)
write.csv(res, file = paste0("CM-MMwithUKBB-univariable-MR-5e-8-", Sys.Date(), ".csv"), quote = F)



## Section below obtains the effect of CM on MM controlling for individual outcomes

#---------------------------------#
#       MVMR for MD (beta2)
#---------------------------------#
# make sure gwas sumstats are consistently named as my_mvmr function is expecting! 
setwd("../output/large-gwas-output-beta2")

### MD ###
md_gwas <- read.delim("../summary-stats/daner_pgc_mdd_meta_w2_rmUKBB_full_logOR.gz")
head(md_gwas)

# format exposure data
md_gwas$Phenotype <- "MD"

# rename columns
names(md_gwas) <- c("CHR", "SNP", "BP", "A1", "A2", 'FRQ_A_116209', 'FRQ_U_314566', 'INFO', 'OR', 'SE', "P", 'ngt', 'Direction', 'HetISqt', 'HetDf', 'HetPVa', 'Nca', 'Nco', 'Neff_half', "BETA", "Phenotype")       

# mvmr CM --> MM, controlling for MD
my_mvmr(exposure1.gwas = CM_gwas, exposure2.gwas = md_gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
        exposure.names = c('cm', 'MD'), exposure.number = 2, pval = 0.000005)


#---------------------------------#
#       MVMR for CAD (beta2)
#---------------------------------#
### CAD ###
cad_gwas <- read.delim("../summary-stats/cad.add.160614.website.txt")
head(cad_gwas)

# format exposure data
cad_gwas$Phenotype <- "CAD"

# rename columns
# current: markername	chr	bp_hg19	effect_allele	noneffect_allele	effect_allele_freq	median_info	model	beta	se_dgc	p_dgc	het_pvalue	n_studies
names(cad_gwas) <- c('SNP','CHR', 'bp_hg19',	'A1',	'A2',	'EAF',	'median_info',	'model',	'BETA',	'SE',	'P',	'het_pvalue',	'n_studies', 'Phenotype')

# mvmr CM --> MM, controlling for CAD
my_mvmr(exposure1.gwas = CM_gwas, exposure2.gwas = cad_gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
        exposure.names = c('cm', 'CAD'), exposure.number = 2, pval = 0.000005)


#---------------------------------#
#       MVMR for T2D (beta2)
#---------------------------------#
### T2D ###
t2d_gwas <- read.delim("../summary-stats/METAANALYSIS_DIAGRAM_SE1_mapped.txt")
head(t2d_gwas)

# format exposure data
t2d_gwas$Phenotype <- "T2D"

# rename columns
names(t2d_gwas) <- c('SNP', 'Chr.Position', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Phenotype')

# mvmr CM --> MM, controlling for T2D
my_mvmr(exposure1.gwas = CM_gwas, exposure2.gwas = t2d_gwas, # note, I am using EU ancestry as ref for clumping here (can change inside function)
        exposure.names = c('cm', 'T2D'), exposure.number = 2, pval = 0.000005)



