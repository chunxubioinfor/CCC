## The Rscript run on the server using 16 cores to calculate the correlation matrix ##
## Assign the method 'spearman' and 'perason' to function ##

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
print(paste('The number of the ligand is',length(ligand)))
receptor <- unique(ligand_receptor$target_genesymbol)
print(paste('The number of the receptor is',length(receptor)))

## Define a function to calculate the geomatrix mean
## Condidering the zero issues
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

## Define a function to calculate scores for the ligand and receptor complex
## A list-like objects of subunits of the complex as the input
complex_expression <- function(subunit){
  exp <- gene_expression[subunit,]
  complex_exp <- apply(exp,2,gm_mean,na.rm=TRUE,zero.propagate = TRUE)
  return(complex_exp)
}

## Initialise the matrix for correlation coefficient and p.value
cor_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(ligand) + length(receptor),
                     dimnames = list(rownames(gene_set_expression),c(ligand,receptor)))
p_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(ligand) + length(receptor),
                   dimnames = list(rownames(gene_set_expression),c(ligand,receptor)))

pb <- txtProgressBar(style=3)
