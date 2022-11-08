# Apply machine learning for determination of correlation mathod and cut-off

library(stringr)
library(ggplot2)
library(Ipaper)

setwd('~/Desktop/CCI_DK/ccc/')
cor_matrix_pearson <- readRDS("~/Desktop/CCI_DK/ccc/cor_matrix_pearson.rds")
cor_matrix_spearman <- readRDS("~/Desktop/CCI_DK/ccc/cor_matrix_spearman.rds")
kegg_symbol_id <- read.table('./KEGG_SYMBOL_ID.txt',header = FALSE,sep = '\t',col.names = c('ID','symbol','category1','category2'))
cellcall <- read.table('./new_ligand_receptor_TFs.txt',header = TRUE)

# Get the available KEGG pathways for validation
# The available pathways occur both in cor matrix and cellcall set
# Get the KEGG pathways in cor matrix
kegg_index <- c()
for (i in 1:nrow(cor_matrix_pearson)) {
  if (str_detect(rownames(cor_matrix_pearson)[i],'KEGG')){
    kegg_index <- append(kegg_index,i)
  }
}
# Transform the KEGG symbols into KEGG IDs
kegg_symbols <- row.names(cor_matrix_pearson)[kegg_index]
kegg_ids <- vector(length = length(kegg_symbols))
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
kegg_df <- data.frame(symbol = kegg_symbols,id = kegg_ids)
# There are some pathways are not transformed automatically
# Manually correct the data.frame
kegg_df <- read.csv('./val_df.csv',header = TRUE)
kegg_df$symbol <- kegg_symbols
# Get the KEGG pathways from cellcall set
kegg_pathway_tmp <- cellcall$pathway_ID
kegg_pathway_cellcall <- c()
for (i in 1:length(kegg_pathway_tmp)){
  kegg_pathway_cellcall <- c(kegg_pathway_cellcall,unlist(strsplit(kegg_pathway_tmp[i],',')))
}
kegg_pathway_cellcall <- unique(kegg_pathway_cellcall)
length(kegg_pathway_cellcall)

# Interset to get the available patwhays for validation 
kegg_pathway_val <- intersect(kegg_pathway_cellcall,kegg_df$id)
print(paste('There are',length(kegg_pathway_val),'available pathways to validate!'))


# Get the available receptors for validation
# The available receptors occur both in cor matrix and cellcall set
# Get the receptors in the cellcall set
receptor_cellcall <- cellcall$Receptor_Symbol
# Interset to get the available receptors for validation
receptor_val <- intersect(receptor_cellcall,receptor)
print(paste('There are',length(receptor_val),'available receptors to validate!'))

# Firstly, for cor_matrix_pearson
# Split the matrix data into a train and test set
train_data_index <- sample(receptor_val,trunc(length(receptor_val)/2))
test_data_index <- setdiff(receptor_val,train_data_index)
kegg_pathway_val_idx <- which(kegg_df$id %in% kegg_pathway_val)

# Define a function to calculate the average ranking given a correlation data.frame and a validation set
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

# There are three ranking ranges:
# 1. The intersected pathways (32)
train_data_pearson <- cor_matrix_pearson[kegg_index,][kegg_pathway_val_idx,train_data_index]
test_data_pearson <- cor_matrix_pearson[kegg_index,][kegg_pathway_val_idx,test_data_index]
print(dim(train_data_pearson))
print(dim(test_data_pearson))
row.names(train_data_pearson) <- kegg_df$id[kegg_pathway_val_idx]
row.names(test_data_pearson) <- kegg_df$id[kegg_pathway_val_idx]
average_ranking <- ranking_cal(train_data_pearson,kegg_val)
average_ranking <- ranking_cal(test_data_pearson,kegg_val)

# 2. All KEGG pathways in data.frame
cor_matrix_pearson_kegg <- cor_matrix_pearson[kegg_index,train_data_index]
row.names(cor_matrix_pearson_kegg) <- kegg_df$id
average_ranking <- ranking_cal(cor_matrix_pearson_kegg,kegg_val)

# 3. All pathways in data.frame
cor_matrix_pearson_all <- cor_matrix_pearson[,train_data_index]
row.names(cor_matrix_pearson_all)[kegg_index] <- kegg_df$id
average_ranking <- ranking_cal(cor_matrix_pearson_all,kegg_val)

# Then for cor_matrix_spearman,
# 1. The intersected pathways (32)
train_data_spearman <- cor_matrix_spearman[kegg_index,][kegg_pathway_val_idx,train_data_index]
test_data_spearman <- cor_matrix_spearman[kegg_index,][kegg_pathway_val_idx,test_data_index]
row.names(train_data_spearman) <- kegg_df$id[kegg_pathway_val_idx]
row.names(test_data_spearman) <- kegg_df$id[kegg_pathway_val_idx]
average_ranking <- ranking_cal(train_data_spearman,kegg_val)
average_ranking <- ranking_cal(test_data_spearman,kegg_val)

# 2. All KEGG pathways in data.frame
cor_matrix_spearman_kegg <- cor_matrix_spearman[kegg_index,train_data_index]
row.names(cor_matrix_spearman_kegg) <- kegg_df$id
average_ranking <- ranking_cal(cor_matrix_spearman_kegg,kegg_val)

# 3. All pathways in data.frame
cor_matrix_spearman_all <- cor_matrix_spearman[,train_data_index]
row.names(cor_matrix_spearman_all)[kegg_index] <- kegg_df$id
average_ranking <- ranking_cal(cor_matrix_spearman_all,kegg_val)

# No train or test data because the validation data is not enough
data_spearman <- cor_matrix_spearman[kegg_index,][kegg_pathway_val_idx,receptor_val]
row.names(data_spearman) <- kegg_df$id[kegg_pathway_val_idx]
average_ranking <- ranking_cal(data_spearman,kegg_val,1)

cor_matrix_spearman_kegg <- cor_matrix_spearman[kegg_index,receptor_val]
row.names(cor_matrix_spearman_kegg) <- kegg_df$id
average_ranking <- ranking_cal(cor_matrix_spearman_kegg,kegg_val,2)

cor_matrix_spearman_all <- cor_matrix_spearman[,receptor_val]
row.names(cor_matrix_spearman_all)[kegg_index] <- kegg_df$id
average_ranking <- ranking_cal(cor_matrix_spearman_all,kegg_val,3)

