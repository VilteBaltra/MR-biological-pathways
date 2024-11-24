library(tidyverse)

# define path to functions
FUNCTIONS = " " 
# path MM
PATH_mm = " "
setwd(PATH_mm)

# path to data 
PATH_A = paste0(PATH_mm, "unadjusted-a-and-b-path-metabolites/metabolites-a-paths-unadjusted-uni-MR-2024-08-09-withnames-nothirdsource.xlsx")
PATH_B_adjusted = paste0(PATH_mm, "adjusted-b-path-metabolites/metabolites-MM-MVMR-adjusted-b-paths-2024-08-30.xlsx")

# read in a and b path data 
a_df <- openxlsx::read.xlsx(PATH_A, sheet = 1, rowNames=F)
a_df$outcome.id <- sub("\\.h$", "", a_df$outcome.id)
b_df_adjusted <- openxlsx::read.xlsx(PATH_B_adjusted, sheet = 1, rowNames=F)
# remove '.tsv.gz' extension from b_df_adjusted$exposure
b_df_adjusted$exposure <- sub("\\.h\\.tsv\\.gz$", "", b_df_adjusted$exposure)

# selected mediators based on previous ansalyses 
selected_mediators <- c("GCST90302074", "GCST90301963", "GCST90301957", "GCST90301946") 
# remove '.tsv.gz' extension from selected_mediators
#selected_mediators <- sub("\\.tsv\\.gz$", "", selected_mediators)

# merge
a_df_selected <- a_df %>% filter(method == "Inverse variance weighted") %>% filter(outcome.id %in% selected_mediators)
b_df_adjusted_selected <- b_df_adjusted %>% filter(exposure %in% selected_mediators)

# ensure order of traits for a and b path are identical (and inspect visually to ensure this)
a_ordered <- a_df_selected[order(a_df_selected$outcome.id), ]
b_ordered <- b_df_adjusted_selected[order(b_df_adjusted_selected$exposure), ]

# combine output into same dataframe
a_and_b_adjusted <- as.data.frame(cbind(a = a_ordered$b, b = b_ordered$b, 
                                        a.se=a_ordered$se, b.se=b_ordered$se, 
                                        a.pval=a_ordered$pval, b.pval=b_ordered$pval))
# add rownames
a_ordered$outcome <- c("Apolipoprotein B", "Esterified cholesterol", "Free cholesterol", "Remnant cholesterol")
rownames(a_and_b_adjusted) = a_ordered$outcome

# source product of coef function
source(paste0(FUNCTIONS, "product.coef.R")) 

mediation <- product.coef(a = a_and_b_adjusted$a, b = a_and_b_adjusted$b, 
                          a.se=a_and_b_adjusted$a.se, b.se=a_and_b_adjusted$b.se, 
                          a.pval=a_and_b_adjusted$a.pval, b.pval=a_and_b_adjusted$b.pval,
                          total.effect = 0.15975) # total effect for CM-MM added for proportion mediated calculation
rownames(mediation) <- rownames(a_and_b_adjusted)
print(mediation) # print all output
mediation[mediation$ab.pval < 0.05,] # no significant metabolite mediators

# save output
openxlsx::write.xlsx(mediation, file = paste0(PATH_mm, "mediation-metabolites/mediation-results-CM-metabolites-MM-", Sys.Date(), ".xlsx"), rowNames = T, overwrite = T)

