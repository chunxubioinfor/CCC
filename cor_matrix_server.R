## The Rscript run on the server using 4 cores to calculate the correlation matrix ##
## Assign the method 'spearman' and 'perason' to function ##
## But the test showed more time cost? I don't know why? ##

library(OmnipathR)
library(dplyr)
library(stringr)
library(parallel)
library(snowfall,lib = '/home/projects/kvs_ccc/R_packages/')
setwd('/home/projects/kvs_ccc/')

# Generation of correlation matrix
## Import the expression matrix
gene_expression_matrix <- readRDS('./output/gene_expression_matrix.rds')
gene_expression_matrix <- as.matrix(gene_expression_matrix)
gene_set_expression_matrix <- readRDS('./output/gene_set_expression_matrix.rds')
icn <- readRDS('./data/icn.rds')
ligand_receptor <- icn %>% dplyr::filter(category_intercell_source == 'ligand',
                                         category_intercell_target == 'receptor') %>%
  dplyr::select(category_intercell_source,source,source_genesymbol,
                category_intercell_target,target,target_genesymbol)
ligand <- unique(ligand_receptor$source_genesymbol)
print(paste('The number of the ligand is',length(ligand)))
receptor <- unique(ligand_receptor$target_genesymbol)
print(paste('The number of the receptor is',length(receptor)))

## Define a function to calculate the geomatrix mean
## Condidering the zero issues
complex_expression <- function(subunit){
  exp <- gene_expression_matrix[subunit,]
  gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
    if(any(x < 0, na.rm = TRUE)){
      return(NaN)
    }
    if(zero.propagate){
      if(any(x == 0, na.rm = TRUE)){
        return(0)
      }
      exp(mean(log(x), na.rm = na.rm))
    } else {
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
  }
  complex_exp <- apply(exp,2,gm_mean,na.rm=TRUE,zero.propagate = TRUE)
  return(complex_exp)
}

## Define a function to calculate scores for the ligand and receptor complex
## A list-like objects of subunits of the complex as the input
cor_p_matrix <- function(lr,gene_expression=gene_expression_matrix,gene_set_expression=gene_set_expression_matrix,method = 'spearman'){
  lr <- unlist(lr)
  complex_expression <- function(subunit){
    exp <- gene_expression_matrix[subunit,]
    gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
      if(any(x < 0, na.rm = TRUE)){
        return(NaN)
      }
      if(zero.propagate){
        if(any(x == 0, na.rm = TRUE)){
          return(0)
        }
        exp(mean(log(x), na.rm = na.rm))
      } else {
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }
    }
    complex_exp <- apply(exp,2,gm_mean,na.rm=TRUE,zero.propagate = TRUE)
    return(complex_exp)
  }	
  cor_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(lr),
                       dimnames = list(rownames(gene_set_expression),lr))
  p_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(lr),
                     dimnames = list(rownames(gene_set_expression),lr))
  loop <- 0
  pb <- txtProgressBar(style=3)
  for (i in 1:ncol(cor_matrix)){
    lr_symbol <- colnames(cor_matrix)[i]
    print(lr_symbol)
    if (str_detect(lr_symbol,'_')){
      subunit <- unlist(strsplit(lr_symbol,split = '_'))
      if (all(subunit %in% rownames(gene_expression))){
        lr_score <- complex_expression(subunit)
      } else {
        loop <- loop + nrow(cor_matrix)  #skip the non-covered columns
        next
      }
    }
    else {
      if (lr_symbol %in% rownames(gene_expression)){
        lr_score <- gene_expression[lr_symbol,]
      } else {
        loop <- loop + nrow(cor_matrix)
        next
      }
    }
    for (j in 1:nrow(cor_matrix)){
      gene_set_symbol <- rownames(cor_matrix)[j]
      gene_set_score <- gene_set_expression[gene_set_symbol,]
      if (method == 'spearman'){
        cor_test <- cor.test(gene_set_score,lr_score,method = 'spearman',exact=FALSE)
      } else if (method == 'pearson') {
        cor_test <- cor.test(gene_set_score,lr_score,method = 'pearson')
      }
      cor <- cor_test$estimate
      p <- cor_test$p.value
      cor_matrix[j,i] <- cor
      p_matrix[j,i] <- p
      loop <- loop + 1
      setTxtProgressBar(pb, loop/(nrow(cor_matrix)*ncol(cor_matrix)))
    }
  }
  result <- list(cor_matrix,p_matrix)
  saveRDS(result,paste('./output/',Sys.time(),'_matrix.rds',sep = ''))
  return(result)
  close(pb)
}

chunk_number <- 6
lr_group <- split(c(ligand,receptor),cut(seq_along(c(ligand,receptor)),chunk_number,labels = FALSE))
cl <- makeCluster(4)
clusterEvalQ(cl,{library(stringr)
  library(dplyr)})
clusterExport(cl,c('gene_expression_matrix','gene_set_expression_matrix'))
print('Preparation is done!')
start_time <- Sys.time()
cor_result <- parLapply(cl,lr_group,cor_p_matrix)
saveRDS(cor_result,'./output/cor_p_matrix.rds')
print('The result list containing both correaltion coefficient and p.value matrix has been saved in ./output/cor_p_matrix.rds')
end_time <- Sys.time()
run_time <- end_time - star_time
print(paste('The time cost is',run_time))
stopCluster(cl)