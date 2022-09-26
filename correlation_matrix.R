## The Rscript run on my Macbook to calculate the correlation matrix ##
## Subset a small matrix to estimate the time cost ##

library(Hmisc)
library(OmnipathR)
library(dplyr)

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

## Establish
cor_matrix <- matrix(nrow = 3152,ncol = 1078,dimnames = list(rownames(gene_set_expression),c(ligand,receptor)))
pb <- txtProgressBar(style=3)
star_time <- Sys.time()
for (i in 1:nrow(cor_matrix)){
  for (j in 1:ncol(cor_matrix)){
    gene_set_score <- gene_set_expression[rownames(cor_matrix)[i],]
    lr_score <- gene_expression[colnames(cor_matrix)[j],]
    cor <- rcorr(gene_set_score,lr_score,type = 'spearman')
    setTxtProgressBar(pb, i/(nrow(cor_matrix)+ncol(cor_matrix)))
  }
}

end_time <- Sys.time()
close(pb)
run_time <- end_time - star_time
