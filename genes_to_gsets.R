library(GSVA)
library(GSEABase)
library(msigdbr)
library(dplyr)

# Import gene sets
KEGG_df_all <-  msigdbr(species = "Homo sapiens",
                        category = "C2",
                        subcategory = "CP:KEGG")
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name)

# 
gene_sets_expression <- gsva(gene_expression,
                             gset.idx.list = kegg_list,
                             method = 'ssgsea',
                             kcdf = 'Poisson',
                             verbose = T)
