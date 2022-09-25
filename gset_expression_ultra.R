## Only gene expression to gene set expression ##

# install and import packages
library(BiocManager)
library(IsoformSwitchAnalyzeR)
library(GSVA)
library(msigdbr)
library(rhdf5)
library(OmnipathR)
library(tidyverse)
library(snow)
library(dynamicTreeCut,lib='/home/projects/kvs_ccc/R_packages/')
library(fastcluster,lib='/home/projects/kvs_ccc/R_packages/')
library(WGCNA,lib='/home/projects/kvs_ccc/R_packages/')
library(clusterProfiler,lib='/home/projects/kvs_ccc/R_packages/')
options(stringsAsFactors = F) 

setwd('/home/projects/kvs_ccc/')
gene_expression <- readRDS('./data/gene_expression_matrix.rds')
print(paste('The dimension of the gene expression matrix is',dim(gene_expression)))
# Quantification from genes to gene-sets level
## Import gene-sets from MSigDB (a combination of three gene sets -- CP + GO slim + Hallmarks)
cp_gene_sets <- clusterProfiler::read.gmt('./data/c2.cp.v2022.1.Hs.symbols.gmt')
h_gene_sets <- clusterProfiler::read.gmt('./data/h.all.v2022.1.Hs.symbols.gmt')
go_gene_sets_df <- OmnipathR::go_annot_slim(organism = 'human',
                                            slim = 'agr',
                                            aspects = c('C','F','P'),
                                            cache = TRUE)
go_gene_sets <- dplyr::select(go_gene_sets_df,go_id,db_object_symbol)
colnames(go_gene_sets) <- c('term','gene')
gene_sets <- rbind(cp_gene_sets,h_gene_sets,go_gene_sets)
gset_list <- split(gene_sets$gene, gene_sets$term)
print('The preparation of gene sets matrix generation is done!')

## Bring the matrix to gene-sets level
start_time_gset <- Sys.time()

gene_expression_matrix <- as.matrix(gene_expression)
gene_set_expression <- gsva(gene_expression_matrix,
                            gset.idx.list = gset_list,
                            method = 'ssgsea',
                            verbose = T,
                            ssgsea.norm = F,
                            parallel.sz = 16)
print('The generation of gene sets expression matrix has been done!')
print(paste('The dimension of the gene expression matrix is',dim(gene_set_expression),sep=' '))
write_rds(gene_set_expression,'./output/gene_set_expression_matrix.rds',compress = 'gz')
print('The gene expression matrix has been saved in /home/projects/kvs_ccc/output/gene_set_expression_matrix.rds.gz')

end_time_gset <- Sys.time()
print(paste('Done! The overall time cost is',end_time_gset - start_time_gset))
