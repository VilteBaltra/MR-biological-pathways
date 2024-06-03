
### Expects GWAS files to have these headers 
### columns should be named "SNP", "CHR","POS", "BETA", "SE", "EA", "NEA", "P", "N", 1e-200, "Z", "INFO"


# ADAPTED mv_extract_exposures_local() function
# clump_data() inside mv_extract_exposures_local() DOES NOT ALWAYS CONNECT TO SERVER AND THAT LEADS TO ISSUES WITH EXPOSURE EXTRACTION 
# adapted version below replaces clump_data() with a local ld_clump() function from "ieugwasr" library 
# more about ld_clump() available here https://mrcieu.github.io/ieugwasr/articles/local_ld.html
# for ld_clump() need library(ieugwasr)
my_mv_extract_exposures_local <- function(filenames_exposure, sep = " ", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se", 
                                          eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval", 
                                          units_col = "units", ncase_col = "ncase", ncontrol_col = "ncontrol", 
                                          samplesize_col = "samplesize", gene_col = "gene", id_col = "id", min_pval = 1e-200, log_pval = FALSE, 
                                          pval_threshold=5e-8, clump_r2=0.001, clump_kb=10000, harmonise_strictness=2)
{
  message("WARNING: Experimental function")
  l_full <- list()
  l_inst <- list()
  for(i in 1:length(filenames_exposure))
  {
    print(filenames_exposure[i])
    l_full[[i]] <- read_outcome_data(filenames_exposure[i], 
                                     sep = sep,
                                     phenotype_col = phenotype_col,
                                     snp_col = snp_col,
                                     beta_col = beta_col,
                                     se_col = se_col,
                                     eaf_col = eaf_col,
                                     effect_allele_col = effect_allele_col,
                                     other_allele_col = other_allele_col,
                                     pval_col = pval_col,
                                     units_col = units_col,
                                     ncase_col = ncase_col,
                                     ncontrol_col = ncontrol_col,
                                     samplesize_col = samplesize_col,
                                     gene_col = gene_col,
                                     id_col = id_col,
                                     min_pval = min_pval,
                                     log_pval = log_pval
    )
    l_inst[[i]] <- subset(l_full[[i]], pval.outcome < pval_threshold)
    l_inst[[i]] <- convert_outcome_to_exposure(l_inst[[i]])
    l_inst[[i]] <- subset(l_inst[[i]], pval.exposure < pval_threshold)
    # l_inst[[i]] <- clump_data(l_inst[[i]], clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb) #VB: replacing with local plink version
    clumped_snps <- ld_clump(
      dplyr::tibble(rsid=l_inst[[i]]$SNP, pval=l_inst[[i]]$pval.exposure, id=l_inst[[i]]$id.exposure),
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = "/Users/vb506/Documents/Projects/MR-mediation-CM-MM/1kg.v3/EUR",
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = pval_threshold,
      pop = "EUR" )
    # VB: keep only clumped snps
    l_inst[[i]] <- subset(l_inst[[i]], SNP %in% clumped_snps$rsid)
  }
  
  exposure_dat <- dplyr::bind_rows(l_inst)
  id_exposure <- unique(exposure_dat$id.exposure)
  temp <- exposure_dat
  temp$id.exposure <- 1
  temp <- temp[order(temp$pval.exposure, decreasing=FALSE), ]
  temp <- subset(temp, !duplicated(SNP))
  #temp <- clump_data(temp, clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb)
  clumped_snps2 <- ld_clump(
    dplyr::tibble(rsid=temp$SNP, pval=temp$pval.exposure, id=temp$id.exposure),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "/Users/vb506/Documents/Projects/MR-mediation-CM-MM/1kg.v3/EUR",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = pval_threshold,
    pop = "EUR" )
  # VB: keep only clumped snps
  temp <- subset(temp, SNP %in% clumped_snps2$rsid)
  exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)
  
  d1 <- lapply(l_full, function(x) {
    subset(x, SNP %in% exposure_dat$SNP)
  }) %>% dplyr::bind_rows()
  
  stopifnot(length(unique(d1$id)) == length(unique(id_exposure)))
  d1 <- subset(d1, mr_keep.outcome)
  d2 <- subset(d1, id.outcome != id_exposure[1])
  d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == id_exposure[1]))
  
  # Harmonise against the first id
  d <- harmonise_data(d1, d2, action=harmonise_strictness)
  
  # Only keep SNPs that are present in all
  tab <- table(d$SNP)
  keepsnps <- names(tab)[tab == length(id_exposure)-1]
  d <- subset(d, SNP %in% keepsnps)
  
  # Reshape exposures
  dh1 <- subset(d, id.outcome == id.outcome[1], select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure))
  dh2 <- subset(d, select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome))
  names(dh2) <- gsub("outcome", "exposure", names(dh2))
  dh <- rbind(dh1, dh2)
  return(dh)
}

# MY MVMR FUNCTION
my_mvmr <- function(exposure1.gwas, exposure2.gwas, exposure3.gwas, exposure4.gwas, 
                    exposure5.gwas, exposure6.gwas,
                    exposure.names, exposure.number, pval = 0.000005){
  
  # Specify the folder name you want to create
  folder_name <- paste0(paste0(exposure.names, collapse = "-"), "-mm")
  
  # Create the folder if it doesn't exist
  if (!file.exists(folder_name)) {
    dir.create(folder_name)
  }
  # Set the newly created folder as the working directory
  setwd(folder_name)
  
  # Initialize a list to store results
  exposure_hits_list <- list()
  
  for(i in 1:exposure.number){
    # select snps suggestive of significance 
    exposure_hits <- get(paste0("exposure", i, ".gwas")) %>% filter(P < pval)
    
    # format exposure
    exposure_hits <- format_data(exposure_hits, snp_col = "SNP",
                                 type = "exposure", beta_col = "BETA",
                                 se_col = "SE", effect_allele_col = "A1",
                                 other_allele_col = "A2",pval_col = "P",
                                 samplesize_col = "N", min_pval = 1e-200,
                                 #z_col = "Z", info_col = "INFO",
                                 chr_col = "CHR", pos_col = "BP")
    
    # Store the results in the list with dynamically generated names
    exposure_hits_list[[paste0("exposure", i, ".gwas.hits")]] <- exposure_hits
  }
  
  # Initialize an empty list to store snp_set for each exposure
  snp_set <- character()
  
  for(i in 1:exposure.number){
    # Create snp_set for the current exposure
    snp_set2 <- c(exposure_hits_list[[paste0("exposure", i, ".gwas.hits")]]$SNP)
    cat(length(snp_set2), "SNPs with p <", pval, "in", exposure.names[[i]], "\n")
    
    # Store the snp_set in the list with a descriptive name
    snp_set <- unique(c(snp_set, snp_set2))
  }
  cat(length(snp_set), "unique SNPs in snp_set\n")
  
  exposure.list <- list()
  for(i in 1:exposure.number){
    # extract snp_set from exposures
    exposure.list[[i]] <- get(paste0("exposure", i, ".gwas")) %>% filter(SNP %in% snp_set)
    cat("Total of", nrow(exposure.list[[i]]), "SNPs extracted from", exposure.names[[i]], "\n")
    write.csv(exposure.list[[i]], file = paste0(exposure.names[[i]], '_mvmr.csv'), row.names = F)
  }
  if(length(exposure.names) == 2){
    exposure_dat  <- my_mv_extract_exposures_local(filenames_exposure = c( paste0(exposure.names[[1]], '_mvmr.csv'),
                                                                        paste0(exposure.names[[2]], '_mvmr.csv') ),
                                                phenotype_col = 'Phenotype',
                                                snp_col = "SNP",
                                                beta_col = "BETA",
                                                se_col = "SE",
                                                effect_allele_col = "A1",
                                                other_allele_col = "A2",
                                                pval_col = "P",
                                                pval_threshold = pval,
                                                harmonise_strictness=2,sep = ',') 
    
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = "~/Documents/Projects/MR-mediation-CM-MM/summary-stats/PCM_multimorbidity_summary_stats_noUKBB_non-het.txt.gz",
      sep = "\t",
      #sep = " ",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "MAF",
      pval_col = "P"
    )
    outcome_dat$outcome <- 'Multimorbidity' 
    
    # Once the data has been obtained, harmonise so that all are on the same reference allele.
    mvdat <- mv_harmonise_data(exposure_dat, outcome_dat) 
    res_mvmr <- mv_multiple(mvdat, pval_threshold=5e-6) # fits all exposures together (recommended as default in TwoSampleMR)
    
    # generate OR
    res_mvmr <- generate_odds_ratios(res_mvmr$result)
    write.csv(res_mvmr, file = paste0(exposure.names[[1]], "-", exposure.names[[2]], "-", "mm-MVMR-results-5e-6", Sys.Date(), ".csv"))
    
  }
  setwd("../")
}

