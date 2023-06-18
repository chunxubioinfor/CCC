## Provide some examples to demonstrate the biological significance of this method
## Pick critical receptors whose corresponding pathways exhibit significant biological information

# load the performance matrix
perf_mtx_R <- readRDS('./perf_mtx_R.rds')
perf_mtx_co_0.5 <- as.data.frame(t(perf_mtx_R[[1]][["perf_matrix_0.5"]]))

# select some well-performed receptors
well_R <- filter(perf_mtx_co_0.5,sensitivity >= 0.8,specificity >= 0.8,FDR <= 0.5)


# Pay attention to some common signaling pathways