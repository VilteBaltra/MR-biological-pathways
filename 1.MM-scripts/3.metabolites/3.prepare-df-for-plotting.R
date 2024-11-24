rm(list=ls())
library(tidyverse)

#---------------------------------------------#
#      OBTAIN a_b_formatted
#---------------------------------------------#
# Calculating indirect effect of all traits in olink-output-mediation folder 

# define path to functions
FUNCTIONS = " " 
# path MM
PATH_mm = " /"

# path to data 
PATH_A = paste0(PATH_mm, "unadjusted-a-and-b-path-metabolites/metabolites-a-paths-unadjusted-uni-MR-2024-08-09-withnames-nothirdsource.xlsx")
PATH_B = paste0(PATH_mm, "unadjusted-a-and-b-path-metabolites/metabolites-b-paths-unadjusted-uni-MR-2024-08-09-nothirdsource.xlsx")

# read in a and b path data 
a_df <- openxlsx::read.xlsx(PATH_A, sheet = 1, rowNames=F)
b_df_adjusted <- openxlsx::read.xlsx(PATH_B, sheet = 1, rowNames=F)
# # remove '.tsv.gz' extension from b_df_adjusted$exposure
# b_df_adjusted$exposure <- sub("\\.tsv\\.gz$", "", b_df_adjusted$exposure)

# merge
a_df_selected <- a_df %>% filter(method == "Inverse variance weighted") 
b_df_adjusted_selected <- b_df_adjusted %>% filter(method == "Inverse variance weighted") 

# ensure order of traits for a and b path are identical (and inspect visually to ensure this)
a_ordered <- a_df_selected[order(a_df_selected$outcome), ]
b_ordered <- b_df_adjusted_selected[order(b_df_adjusted_selected$exposure), ]
all(a_ordered$outcome == b_ordered$exposure) # must be TRUE

# combine output into same dataframe
a_and_b_adjusted <- as.data.frame(cbind(a = a_ordered$b, b = b_ordered$b, 
                                        a.se=a_ordered$se, b.se=b_ordered$se, 
                                        a.pval=a_ordered$pval, b.pval=b_ordered$pval))

# add rownames
rownames(a_and_b_adjusted) = a_ordered$outcome

# source product of coef function --> ONLY USING THIS NOW TO GET THE DF FORMAT (NOT RELEVANT FOR MEDIATION AS USING A AND B PATHS UNADJUSTED)
source(paste0(FUNCTIONS, "product.coef.R")) 

mediation <- product.coef(a = a_and_b_adjusted$a, b = a_and_b_adjusted$b, 
                          a.se=a_and_b_adjusted$a.se, b.se=a_and_b_adjusted$b.se, 
                          a.pval=a_and_b_adjusted$a.pval, b.pval=a_and_b_adjusted$b.pval,
                          total.effect = 0.15975) # total effect for CM-MM added for proportion mediated calculation
rownames(mediation) <- rownames(a_and_b_adjusted)
print(mediation) # print all output
mediation[mediation$ab.pval < 0.05,] # no significant metabolite mediators

# save output
saveRDS(mediation, file = paste0(PATH_mm, "data_a_b_paths_formatted_metabolites.rds")) 
