# 

library(stringr)

setwd('~/Desktop/CCI_DK/ccc/')
cor_matrix_pearson <- readRDS("~/Desktop/CCI_DK/ccc/cor_matrix_pearson.rds")
cor_matrix_spearman <- readRDS("~/Desktop/CCI_DK/ccc/cor_matrix_spearman.rds")
kegg_index <- c()
for (i in 1:nrow(cor_matrix_pearson)) {
  if (str_detect(rownames(cor_matrix_pearson)[i],'KEGG')){
    kegg_index <- append(kegg_index,i)
}
}

# Split the matrix data into a train and test set
train_data_index <- sample(receptor,trunc(length(receptor)/2))
test_data_index <- setdiff(receptor,train_data_index)
train_data <- cor_matrix_pearson[kegg_index,train_data_index]
test_data <- cor_matrix_pearson[kegg_index,test_data_index]
print(dim(train_data))
print(dim(test_data))

# Import the curated validation data
kegg_symbol_id <- read.table('./KEGG_SYMBOL_ID.txt',header = FALSE,sep = '\t',col.names = c('ID','symbol','category1','category2'))
kegg_val <- read.table('./new_ligand_receptor_TFs.txt',header = TRUE)

# Transform the KEGG symbols into KEGG IDs
kegg_symbols <- row.names(train_data)
kegg_ids <- vector(length = nrow(train_data))
for (i in 1:length(kegg_symbols)){
  kegg_symbols[i] <- str_replace_all(str_replace_all(kegg_symbols[i],'KEGG_',''),'_',' ')
  kegg_symbols[i] <- str_to_sentence(str_to_lower(kegg_symbols[i]))
  tryCatch({
    kegg_ids[i] <- kegg_symbol_id[kegg_symbol_id$symbol == kegg_symbols[i],]$ID
  }, error = function(e){
    print(e)
    kegg_ids[i] <- kegg_symbols[i]
  })
}
val_df <- data.frame(symbol = kegg_symbols,id = kegg_ids)
# There are some pathways are not transformed automatically
# Manually correct the data.frame
val_df <- read.csv('./val_df.csv',header = TRUE)
val_df$symbol <- row.names(train_data)

# Get the intersection of unique KEGG pathways in the curated validation data and KEGG pathways of GSEA 
# Get the intersection of receptors
kegg_pathway <- c()
for (i in 1:length(kegg_pathway_tmp)){
  kegg_pathway <- c(kegg_pathway,unlist(strsplit(kegg_pathway_tmp[i],',')))
}
kegg_pathway <- unique(kegg_pathway)
length(kegg_pathway)
val_kegg_pathway <- intersect(kegg_pathway,val_df$id)
print(paste('There are',length(val_kegg_pathway),'available pathways to validate!'))

# Get the 
