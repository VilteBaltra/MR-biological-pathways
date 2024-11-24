
setwd("/Users/vb506/Documents/Projects/MR-mediation-CM-MM/output/OUTPUT-2024/DEPRESSION/other-traits/unadjusted-b-path-md/")

# b path combine
matabolites.b <- openxlsx::read.xlsx("lipids-large-with-and-withoutUKBB-DEPRESSION-uniMR-b-paths-2024-09-15.xlsx", sheet = 1)
matabolites.b.intercept <- openxlsx::read.xlsx("lipids-large-with-and-withoutUKBB-DEPRESSION-uniMR-b-paths-2024-09-15.xlsx", sheet = 2)
matabolites.b.Q <- openxlsx::read.xlsx("lipids-large-with-and-withoutUKBB-DEPRESSION-uniMR-b-paths-2024-09-15.xlsx", sheet = 3)


head(matabolites.b)
head(matabolites.b.intercept)
head(matabolites.b.Q)

colnames(matabolites.b.intercept) <- paste0(colnames(matabolites.b.intercept), ".intercept")

# combine output for supplementary tables
merged1b <- merge(matabolites.b, matabolites.b.intercept[, c('exposure.intercept', 'egger_intercept.intercept', 'se.intercept', 'pval.intercept')], by.x = 'exposure', by.y = 'exposure.intercept')
merged2b <- merge(merged1b, matabolites.b.Q[, c( 'exposure', 'method', 'Q', 'Q_df', 'Q_pval')], by.x = c('exposure', 'method'), by.y = c('exposure', 'method'), all = TRUE)

ind.b <- which(merged2b$method != "Inverse variance weighted")

merged2b[ind.b, c('egger_intercept.intercept', 'se.intercept', 'pval.intercept')] <- NA
head(merged2b)

openxlsx::write.xlsx(merged2b, file =  paste0('combined-uniMR-b-paths-new-lipids-DEPRESSION-', Sys.Date(), ".xlsx"), rowNames = F, overwrite = T)


