library(IsoformSwitchAnalyzeR)

# Example codes 
salmonQuant <- importIsoformExpression(system.file("extdata/", package="IsoformSwitchAnalyzeR"))
geneRepCount <- isoformToGeneExp(
  isoformRepExpression  = salmonQuant$counts,
  isoformGeneAnnotation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR")
)

# 
gene_expression <- IsoformSwitchAnalyzeR::isoformToGeneExp(isoform_expression,
                                        isoformGeneAnnotation = '/home/databases/archs4/v11/Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf',
                                        )
