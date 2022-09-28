## The Rscript run on my Macbook to calculate the correlation matrix ##
## Make some alterations to the order of the for and if loop ##
## Dramatically save time and space ##
## 14min for 10 samples and all gene sets ##

# In cooperation with the script file correlation_matrix.R

cor_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(ligand) + length(receptor),
                     dimnames = list(rownames(gene_set_expression),c(ligand,receptor)))
p_matrix <- matrix(nrow = nrow(gene_set_expression),ncol = length(ligand) + length(receptor),
                   dimnames = list(rownames(gene_set_expression),c(ligand,receptor)))
pb <- txtProgressBar(style=3)
print('Begin!')
star_time <- Sys.time()
loop <- 0

for (i in 1:ncol(cor_matrix)){
  lr_symbol <- colnames(cor_matrix)[i]
  if (str_detect(lr_symbol,'_')){
    subunit <- unlist(strsplit(lr_symbol,split = '_'))
    if (all(subunit %in% rownames(gene_expression))){
      lr_score <- complex_expression(subunit)
    } else {
      loop <- loop + nrow(cor_matrix)
      next
    }
  }
  else {
    if (lr_symbol %in% rownames(gene_expression)){
      lr_score <-gene_expression[lr_symbol,]
    } else {
      loop <- loop + nrow(cor_matrix)
      next
    }
  }
  for (j in 1:nrow(cor_matrix)){
    gene_set_symbol <- rownames(cor_matrix)[j]
    gene_set_score <- gene_set_expression[gene_set_symbol,]
    cor_test <- cor.test(gene_set_score,lr_score,method = 'spearman',exact=FALSE)
    cor <- cor_test$estimate
    p <- cor_test$p.value
    cor_matrix[j,i] <- cor
    p_matrix[j,i] <- p
    loop <- loop + 1
    setTxtProgressBar(pb, loop/(nrow(cor_matrix)*ncol(cor_matrix)))
  }
}

print('All is done!')
end_time <- Sys.time()
close(pb)
run_time <- end_time - star_time
print(run_time)
