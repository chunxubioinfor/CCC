## The Rscript run on my Macbook to calculate the correlation matrix ##
## Subset a small matrix to estimate the time cost ##
## But it is so slow, so I'm thinking about trying something new ##
## Btw the time cost: 1.3h for only 10 samples (10/16000+) and 30 gene sets (30/3152) ##


library(Hmisc)
library(OmnipathR)
library(dplyr)
library(stringr)

setwd('/home/projects/kvs_ccc/')
# Generation of correlation matrix
## Import the expression matrix
gene_expression <- readRDS('./output/gene_expression_matrix.rds')
gene_set_expression <- readRDS('./output/gene_set_expression_matrix.rds')
icn <- import_intercell_network(high_confidence = TRUE)
ligand_receptor <- icn %>% dplyr::filter(category_intercell_source == 'ligand',
                                         category_intercell_target == 'receptor') %>%
dplyr::select(category_intercell_source,source,source_genesymbol,
              category_intercell_target,target,target_genesymbol)
ligand <- unique(ligand_receptor$source_genesymbol)
receptor <- unique(ligand_receptor$target_genesymbol)

## Establish the spearman correlation matrix
cor_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(ligand) + length(receptor),
                     dimnames = list(rownames(gene_set_expression),c(ligand,receptor)))
p_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(ligand) + length(receptor),
                   dimnames = list(rownames(gene_set_expression),c(ligand,receptor)))
pb <- txtProgressBar(style=3)
print('Begin!')
star_time <- Sys.time()
loop <- 0
for (i in 1:nrow(cor_matrix)){
  for (j in 1:ncol(cor_matrix)){
    lr_symbol <- colnames(cor_matrix)[j]
    gene_set_symbol <- rownames(cor_matrix)[i]
    if (str_detect(lr_symbol,'_')){
      subunit <- unlist(strsplit(lr_symbol,split = '_'))
      if (all(subunit %in% rownames(gene_expression))){
        lr_score <- complex_expression(subunit)
      } else {
        next
      }
    }
    else {
      if (lr_symbol %in% rownames(gene_expression)){
        lr_score <-gene_expression[lr_symbol,]
      } else {
        next
      }
    }
    gene_set_score <- gene_set_expression[rownames(cor_matrix)[i],]
    cor_test <- cor.test(gene_set_score,lr_score,method = 'spearman',exact=FALSE)
    cor <- cor_test$estimate
    p <- cor_test$p.value
    cor_matrix[i,j] <- cor
    p_matrix[i,j] <- p
    loop <- loop + 1
    setTxtProgressBar(pb, loop/(30*ncol(cor_matrix)))
  }
  print('A row is done!')
}

print('All is done!')
end_time <- Sys.time()
close(pb)
run_time <- end_time - star_time

### Calculation of geomatrix mean
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

### For the ligand and receptor complex
complex_expression <- function(subunit){
  exp <- gene_expression[subunit,]
  complex_exp <- apply(exp,2,gm_mean,na.rm=TRUE,zero.propagate = TRUE)
  return(complex_exp)
}


