library(tidyverse)

#---------------------------------------------#
#      CALCULATE INDIRECT EFFECT
#---------------------------------------------#
# Calculating indirect effect of all traits in olink-output-mediation folder 

# define path to functions
FUNCTIONS = " " 
# path MM
PATH_mm = " "

# path to data 
PATH_A = paste0(PATH_mm, "../MULTIMORBIDITY/metabolites/unadjusted-a-and-b-path-metabolites/metabolites-a-paths-unadjusted-uni-MR-2024-08-09-withnames-nothirdsource.xlsx") # same as for MM
PATH_B_adjusted = paste0(PATH_mm, "metabolites/metabolites-MM-MVMR-adjusted-b-paths-DEPRESSION-2024-09-15.xlsx")

# read in a and b path data 
a_df <- openxlsx::read.xlsx(PATH_A, sheet = 1, rowNames=F)
b_df_adjusted <- openxlsx::read.xlsx(PATH_B_adjusted, sheet = 1, rowNames=F)
b_df_adjusted$exposure <- sub("\\.tsv\\.gz$", "", b_df_adjusted$exposure)

# selected mediators based on previous ansalyses (only IL.6-1 was significant in MVMR)
selected_mediators <- c('GCST90301955.h')

# merge
a_df_selected <- a_df %>% filter(method == "Inverse variance weighted") %>% filter(outcome.id %in% selected_mediators)
b_df_adjusted_selected <- b_df_adjusted %>% filter(exposure %in% selected_mediators)

# ensure order of traits for a and b path are identical (and inspect visually to ensure this)
# (note: does nothing here as only one trait)
a_ordered <- a_df_selected[order(a_df_selected$outcome.id), ]
b_ordered <- b_df_adjusted_selected[order(b_df_adjusted_selected$exposure), ]

# combine output into same dataframe
a_and_b_adjusted <- as.data.frame(cbind(a = a_ordered$b, b = b_ordered$b, 
                          a.se=a_ordered$se, b.se=b_ordered$se, 
                          a.pval=a_ordered$pval, b.pval=b_ordered$pval))
# add rownames
rownames(a_and_b_adjusted) = a_ordered$outcome

# source product of coef function
source(paste0(FUNCTIONS, "product.coef.R")) 

mediation <- product.coef(a = a_and_b_adjusted$a, b = a_and_b_adjusted$b, 
                          a.se=a_and_b_adjusted$a.se, b.se=a_and_b_adjusted$b.se, 
                          a.pval=a_and_b_adjusted$a.pval, b.pval=a_and_b_adjusted$b.pval,
                          total.effect = 0.45344289) # total effect for CM-MD added for proportion mediated calculation
rownames(mediation) <- rownames(a_and_b_adjusted)
print(mediation) # print all output
mediation[mediation$ab.pval < 0.05,] # no significant mediators

# save output
openxlsx::write.xlsx(mediation, file = paste0(PATH_mm, "mediation-results-CM-metabolites-MD-adjusted-", Sys.Date(), ".xlsx"), rowNames = T, overwrite = T)
