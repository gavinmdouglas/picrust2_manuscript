### Commands to write the DA results to an excel spreadsheet, which is easier to explore.

library(openxlsx)

full_DA <- list()

for(genes in names(go_enrichment_results)) {
  full_GO[[genes]] <- go_enrichment_results[[genes]]$all_results
}

write.xlsx(full_GO, file = "tables/gene_set_full_GO_enrichment.xlsx")

