## This is a R script used to test the feasiblity and estimate the running time ##
## We use 1000 sampels (sample names are stored in sample_1k.csv) and apply 1 core of the PC ##

# install and import packages
library(BiocManager)
library(IsoformSwitchAnalyzeR)
library(GSVA)
library(GSVAdata)
library(msigdbr)
library(rhdf5)
library(biomaRt)
library(GSA)
library(clusterProfiler)
library(OmnipathR)

options(stringsAsFactors = F) 

start_time <- Sys.time()
setwd('~/Desktop/CCI_DK/ccc/')
dir.create('./test_output')
ARCHS4_transcript_file = "~/Desktop/CCI_DK/ccc/human_transcript_v11_tpm.h5"
sample_1k <- read.csv('~/Desktop/CCI_DK/ccc/sample_1k.csv',col.names = 'sample_name')$sample_name


# Extract 1k samples from ARCHS4_transcript_file
human_transcript_v11 <- H5Fopen(ARCHS4_transcript_file)
h5dump(human_transcript_v11,load=FALSE)    # Dump the content of an HDF5 file
samples <- human_transcript_v11$meta$samples$geo_accession
transcripts <- human_transcript_v11$meta$transcripts$transcripts
sample_locations = which(samples %in% sample_1k)
transcript_expression <- t(h5read(human_transcript_v11, "data/expression", index=list(sample_locations, 1:length(transcripts))))
H5close()
rownames(transcript_expression) <- transcripts
colnames(transcript_expression) <- samples[sample_locations]
write.csv(transcript_expression,'./test_output/transcript_expression.csv')
print('The generation of transcripts expression matrix has been done!')
saveRDS(transcript_expression,'./test_output/transcript_expression_matrix.rds')


# Quantification from transcripts to genes level
## Import GTF file
aSwitchList <- importGTF(pathToGTF='./Homo_sapiens.GRCh38.87.chr_patch_hapl_scaff.gtf')

## Create dataframe with associations
isoGene <- unique(aSwitchList$isoformFeatures[c('isoform_id','gene_name')])
colnames(isoGene) <- c('isoform_id','gene_id')

start_time_gene <- Sys.time()
gene_expression <- IsoformSwitchAnalyzeR::isoformToGeneExp(transcript_expression,
                                                           isoformGeneAnnotation = isoGene,
                                                           quiet = FALSE)
end_time_gene <- Sys.time()
print(start_time_gene - end_time_gene)
#write.csv(gene_expression,'./test_output/gene_expression.csv')
print('The generation of genes expression matrix has been done!')
#saveRDS(gene_expression,'./test_output/gene_expression_matrix.rds')


# Quantification from genes to gene-sets level
## Import gene-sets from MSigDB (a combination of three gene sets -- CP + GO(BP+MF+CC) + Hallmarks)
cp_gene_sets <- clusterProfiler::read.gmt('./c2.cp.v2022.1.Hs.symbols.gmt')
h_gene_sets <- clusterProfiler::read.gmt('./h.all.v2022.1.Hs.symbols.gmt')
#go_gene_sets <- clusterProfiler::read.gmt('c5.go.v2022.1.Hs.symbols.gmt')
go_gene_sets_df <- OmnipathR::go_annot_slim(organism = 'human',
                                         slim = 'agr',
                                         aspects = c('C','F','P'),
                                         cache = TRUE)
go_gene_sets <- dplyr::select(go_gene_sets_df,go_id,db_object_symbol)
colnames(go_gene_sets) <- c('term','gene')
gene_sets <- rbind(cp_gene_sets,h_gene_sets,go_gene_sets)
gset_list <- split(gene_sets$gene, gene_sets$term)

## Bring the matrix to genesets level
start_time_gset <- Sys.time()
gene_expression_matrix <- as.matrix(gene_expression)
gene_set_expression <- gsva(gene_expression_matrix,
                            gset.idx.list = gset_list,
                            method = 'ssgsea',
                            kcdf = 'Poisson',
                            verbose = T,
                            ssgsea.norm = F,
                            parallel.sz = 16)
#write.csv(gene_set_expression,'./gset_expression.csv')
print('The generation of gene sets expression matrix has been done!')

end_time_gset <- Sys.time()
print(end_time_gset - start_time_gset)


