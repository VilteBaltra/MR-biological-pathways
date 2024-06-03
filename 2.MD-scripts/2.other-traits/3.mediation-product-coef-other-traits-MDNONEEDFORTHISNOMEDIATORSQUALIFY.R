library(tidyverse)

#---------------------------------------------#
#      CALCULATE INDIRECT EFFECT
#---------------------------------------------#
# Calculating indirect effect of all traits in olink-output-mediation folder 

# define path to functions
FUNCTIONS = "~/Documents/Projects/MR-mediation-CM-MM/scripts/TIDY-WORKFLOW/source/" 
# path MM
PATH_mm = "~/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/DEPRESSION/"

# path to data 
PATH_A = paste0(PATH_mm, "../MULTIMORBIDITY/other-traits/unadjusted-a-and-b-paths/a-and-b-paths-unadjusted-uni-MR-2023-10-12.xlsx") # same as for MM
PATH_B_adjusted = paste0(PATH_mm, "other-traits/adjusted-b-path-md/other-traits-MD-MVMR-adjusted-b-paths-2024-01-30.xlsx")

# read in a and b path data 
a_df <- openxlsx::read.xlsx(PATH_A, sheet = 1, rowNames=F)
b_df_adjusted <- openxlsx::read.xlsx(PATH_B_adjusted, sheet = 1, rowNames=F)

# selected mediators based on previous analyses (only FI-ebi-a-GCST90002238 was significant in MVMR)
selected_mediators <- c("Fasting insulin || id:ebi-a-GCST90002238")

# merge
a_df_selected <- a_df %>% filter(method == "Inverse variance weighted") %>% filter(outcome %in% selected_mediators)
b_df_adjusted_selected <- b_df_adjusted %>% filter(exposure %in% "FI-ebi-a-GCST90002238")
# name same way
b_df_adjusted_selected$exposure <- "Fasting insulin || id:ebi-a-GCST90002238"

# ensure order of traits for a and b path are identical (and inspect visually to ensure this)
# (note: here only one trait so does nothing)
a_ordered <- a_df_selected[order(a_df_selected$outcome), ]
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
                          total.effect = 0.45344289) # total effect for CM-MM added for proportion mediated calculation
rownames(mediation) <- rownames(a_and_b_adjusted)
print(mediation) # print all output
mediation[mediation$ab.pval < 0.05,] # no significant mediators

# save output
openxlsx::write.xlsx(mediation, file = paste0(PATH_mm, "other-traits/mediation-other-traits-md/mediation-results-CM-other-traits-MD-adjusted-", Sys.Date(), ".xlsx"), rowNames = T, overwrite = T)
