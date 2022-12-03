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
  select(gsea_symbol,id,receptor)  ## Filter out the invalid rows 
curation_pathway_wiki <- read.xls('./curation_pathway.xlsx',sheet = 2,na.strings=c("NA","#DIV/0!",''))
curation_pathway_wiki <- curation_pathway_wiki %>% filter(receptor != 'NA') %>%
  select(gsea_symbol,id,receptor)
dim(curation_pathway_kegg)
dim(curation_pathway_wiki)
print(paste('The number of valid pathways is',nrow(curation_pathway_kegg)+nrow(curation_pathway_wiki)))

# Reshape the dateset
curation_pathway_kegg <- curation_pathway_kegg %>% as_tibble() %>% separate_rows(receptor,sep = ',')
curation_pathway_wiki <- curation_pathway_wiki %>% as_tibble() %>% separate_rows(receptor,sep = ',')
curation_pathway <- rbind(curation_pathway_kegg,curation_pathway_wiki)

# Filter out the non-receptor (maybe it is caused by the curation mistakes)
# 
cor_matrix_spearman <- readRDS("~/Desktop/CCI_DK/ccc/cor_matrix_spearman.rds")
receptor_ref <- colnames(cor_matrix_spearman)
curation_pathway_filtered <- filter(curation_pathway,receptor %in% receptor_ref)
write.csv(curation_pathway_filtered,'./curation_pathway_filtered.csv')

# Calculation of average ranking
# Define a function to calculate the average ranking given a correlation data.frame and a validation set

ranking_cal <- function(ranking_df,curation_pathway,file_name){
  average_ranking <- c()
  for (i in 1:ncol(ranking_df)){
    receptor <- colnames(ranking_df)[i]
    receptor_curation_df <- filter(curation_pathway,receptor = receptor)
    pathway_2_receptor <- unique(receptor_curation_df$gsea_symbol)
    ranking <- c()
    for (j in 1:length(pathway_2_receptor)){
      pathway <- pathway_2_receptor[j]
      ranking <- c(ranking,nrow(ranking_df) +1 - rank(ranking_df[,i])[pathway])
    }
    ranking <- mean(ranking)
    average_ranking <- c(average_ranking,ranking)
  }
  print(mean(average_ranking,na.rm = TRUE))
  print(matrix(average_ranking,nrow = 18))
}


ranking_cal <- function(ranking_df,kegg_val,file_name){
  average_ranking <- c()
  for (i in 1:length(colnames(ranking_df))){
    receptor <- colnames(ranking_df)[i]
    receptor_cellcall_df <- kegg_val[kegg_val$Receptor_Symbol == receptor,]
    pathway_cellcall_tmp <- receptor_cellcall_df$pathway_ID
    pathway_cellcall <- c()
    for (k in 1:length(pathway_cellcall_tmp)){
      pathway_cellcall <- c(pathway_cellcall,unlist(strsplit(pathway_cellcall_tmp[k],',')))
    }
    pathway_cellcall <- intersect(unique(pathway_cellcall),kegg_df$id)
    ranking <- c()
    for (j in 1:length(pathway_cellcall)){
      pathway <- pathway_cellcall[j]
      ranking <- c(ranking,nrow(ranking_df) +1 - rank(ranking_df[,i])[pathway])
    }
    ranking <- mean(ranking)
    average_ranking <- c(average_ranking,ranking)
  }
  print(mean(average_ranking,na.rm = TRUE))
  print(matrix(average_ranking,nrow = 18))
  average_ranking_df <- data.frame(receptor = colnames(ranking_df),average_ranking)
  p <- ggplot(average_ranking_df, aes(x = average_ranking)) +
    geom_density(color = 'black', fill = 'gray') +
    geom_vline(xintercept = mean(average_ranking,na.rm = TRUE),color = 'red')
  p
  h <- ggplot(average_ranking_df,aes(x = average_ranking)) + 
    geom_histogram(bins = 40) + 
    geom_vline(xintercept = mean(average_ranking,na.rm = TRUE),color = 'red')
  write_fig(h,paste('./',file_name,'.png',sep = ''))
  return(average_ranking)
}

