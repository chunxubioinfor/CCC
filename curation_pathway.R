library(rjson)
library(jsonlite)
setwd('~/Desktop/CCI_DK/ccc/')
# Curation of KEGG pathways
gsea_kegg <- fromJSON("./kegg.json",simplifyDataFrame = TRUE)
kegg_symbols <- names(gsea_kegg)
kegg_ids <- c()
for (i in 1:length(kegg_symbols)){
  kegg_ids <- c(kegg_ids,gsea_kegg[[kegg_symbols[i]]]$exactSource)
}
kegg_symbol_id <- data.frame(kegg_symbols,kegg_ids)
write.csv(kegg_symbol_id,'./tmp.csv')

# Curation of WIKI pathway
gsea_wiki <- fromJSON("./wiki.json",simplifyDataFrame = TRUE)
wiki_symbols <- names(gsea_wiki)
wiki_ids <- c()
for (i in 1:length(wiki_symbols)){
  wiki_ids <- c(wiki_ids,gsea_wiki[[wiki_symbols[i]]]$exactSource)
}
wiki_symbol_id <- data.frame(wiki_symbols,wiki_ids)
write.csv(wiki_symbol_id,'./tmp.csv')
