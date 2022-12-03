library(rjson)
library(jsonlite)
library(gdata)
library(dplyr)
library(tidyr)
setwd('~/Desktop/CCI_DK/ccc/')
options(stringsAsFactors = FALSE)
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


# Then manually curate the KEGG pathway and WIKi pathway
# This step took me aroud two weeks
# The results were saved in curation_pathway.csv
# Read into the file and reshape the dataset
curation_pathway_kegg <- read.xls('./curation_pathway.xlsx',sheet = 1,na.strings=c("NA","#DIV/0!",''))
curation_pathway_kegg <- curation_pathway_kegg %>% filter(receptor != 'NA') %>%
  select(gsea_symbol,kegg_id,receptor)  ## Filter out the invalid rows 
curation_pathway_wiki <- read.xls('./curation_pathway.xlsx',sheet = 2,na.strings=c("NA","#DIV/0!",''))
curation_pathway_wiki <- curation_pathway_wiki %>% filter(receptor != 'NA') %>%
  select(gsea_symbol,wiki_id,receptor)
dim(curation_pathway_kegg)
dim(curation_pathway_wiki)
print(paste('The number of valid pathways is',nrow(curation_pathway_kegg)+nrow(curation_pathway_wiki)))

# Reshape the dateset
# 
curation_pathway_kegg <- curation_pathway_kegg %>% as_tibble() %>% separate_rows(receptor,sep = ',')
curation_pathway_wiki <- curation_pathway_wiki %>% as_tibble() %>% separate_rows(receptor,sep = ',')
