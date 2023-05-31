library(readxl)
library(dplyr)
library(tidyr)
library(progress)

# Based on the Average Ranking result in which the plot is extremely left skewed,
# we decided to develop an algorithm applied onto cor_matrix to determine whether correlated
# Here, we introduce two types of algorithms:
# 1. A more straightforward one -- Iteration of the cut-off value
# 2. A more relative one -- Iteration of top value

## Randomly divide the curation dataset into two groups, testing data (2/3) and validation data (1/3)
curation_pathway <- read_excel('./curation_pathway.xlsx',sheet = "Curation_Pathway_for_R",col_names = TRUE)
curation_pathway <- curation_pathway %>% as_tibble() %>% separate_rows(receptors,sep = ',')  #reshape the dataset
curation_pathway_filtered <- filter(curation_pathway,receptors %in% ligand_receptor)
curation_pathway_filtered$index <- 1:nrow(curation_pathway_filtered)
# write.csv(curation_pathway_filtered,'./curation_pathway_filtered.csv')
validation_pathway <- curation_pathway_filtered %>% group_by(source) %>% sample_frac(size = 1/3) #leave alone for final validation
testing_pathway<- anti_join(curation_pathway_filtered,validation_pathway,by = 'index')


## Apply two algorithms to calculate the performance matrix

# Define a function which transfer the correaltion coefficient matrix into correlated or not matrix
correlated_or_not_mtx <- function(cor_matrix, cut_off_value){
  cor_matrix[cor_matrix >= cut_off_value] <- 'TRUE'
  cor_matrix[cor_matrix < cut_off_value] <- 'FALSE'
  return(cor_matrix)
}

correlated_or_not_mtx <- function(cor_matrix, algorithm, value,view){
  if(algorithm == 'cut_off'){
    cut_off_value <- value
    cor_matrix[cor_matrix >= cut_off_value] <- 'TRUE'
    cor_matrix[cor_matrix < cut_off_value] <- 'FALSE'
    return(cor_matrix)
  }
  else if(algorithm == 'top'){
    top_value <- value
    cor_or_not_mtx <- matrix(NA, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix),
                             dimnames = list(rownames(cor_matrix),colnames(cor_matrix)))
    if(view == 'receptor'){
      for (i in 1:ncol(cor_matrix)) {
        col <- cor_matrix[, i]
        sorted_col <- sort(col, decreasing = TRUE)
        threshold <- sorted_col[top_value]
        cor_or_not_mtx[, i] <- col >= threshold
      }
    }
    else if(view == 'pathway'){
      for (i in 1:nrow(cor_matrix)) {
        row <- cor_matrix[i,]
        sorted_row <- sort(row, decreasing = TRUE)
        threshold <- sorted_row[top_value]
        cor_or_not_mtx[i,] <- row >= threshold
      }
    }
    return(cor_or_not_mtx)
  }
}
# Define a function which calculate and output performance matrix
perf_mtx <- function(cor_or_not_mtx, view, curation_pathway){
  if(view == "receptor"){
    perf_matrix <- matrix(0,nrow = 3,ncol = ncol(cor_or_not_mtx),dimnames = list(c('sensitivity','specificity','FDR'),colnames(cor_or_not_mtx)))
    for(i in 1:ncol(cor_or_not_mtx)){
      correlated_pathway <- rownames(cor_or_not_mtx)[cor_or_not_mtx[,i] == "TRUE"]
      ref_pathway <- filter(curation_pathway,receptors == colnames(cor_or_not_mtx)[i])$gsea_symbol
      tp <- length(intersect(correlated_pathway,ref_pathway))
      fp <- length(correlated_pathway) -tp
      tn <- nrow(cor_or_not_mtx) - length(ref_pathway) - fp
      sensitivity <- tp/length(ref_pathway)
      specificity <- tn/(nrow(cor_or_not_mtx) - length(ref_pathway))
      FDR <- fp/length(correlated_pathway)
      perf_matrix[,i] <- c(sensitivity,specificity,FDR)
    }
  }
  else if(view == 'pathway'){
    perf_matrix <- matrix(0,nrow = nrow(cor_or_not_mtx),ncol = 3,dimnames = list(rownames(cor_or_not_mtx),c('sensitivity','specificity','FDR')))
    for(i in 1:nrow(cor_or_not_mtx)){
      correlated_receptor <- colnames(cor_or_not_mtx)[cor_or_not_mtx[i,] == "TRUE"]
      ref_receptor <- filter(curation_pathway,gsea_symbol == rownames(cor_or_not_mtx)[i])$receptors
      tp <- length(intersect(correlated_receptor,ref_receptor))
      fp <- length(correlated_receptor) -tp
      tn <- nrow(cor_or_not_mtx) - length(ref_receptor) - fp
      sensitivity <- tp/length(ref_receptor)
      specificity <- tn/(nrow(cor_or_not_mtx) - length(ref_receptor))
      FDR <- fp/length(correlated_receptor)
      perf_matrix[i,] <- c(sensitivity,specificity,FDR)
  }
  }
  return(perf_matrix)
}

# Test two functions
# t <- correlated_or_not_mtx(cor_matrix_spearman,0.5)
# a <- perf_mtx(t,'pathway',curation_pathway_filtered)
# 1. Algorithm 1: Iteration of the cut-off value
perf_matrix_list <- list()
pb <- progress_bar$new(total = 15)
for (cut_off in seq(0.2, 0.9, 0.05)) {
  cor_or_not_mtx <- correlated_or_not_mtx(cor_matrix_spearman,'cut_off',cut_off,'receptor')
  perf_matrix <- perf_mtx(cor_or_not_mtx,'receptor',curation_pathway_filtered)
  matrix_name <- paste('perf_matrix',cut_off,sep = '_')
  perf_matrix_list[[matrix_name]] <- perf_matrix
  pb$tick()
}



