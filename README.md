# CCC
This is a graduate project by Chunxu Han, who study Bioinformatics in DTU. 


### File

| **File**                        | **Description**                                                                                                                               |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| *gene_sets_exp_matrix_breast.R* | Test R script run on Macbook, so just cover lung samples.                                                                                     |
| *gene_sets_exp_matrix.R*        | Test R script run on server, cover all samples.                                                                                               |
| *gset_test_1k_1core.R*          | The R script used to test the feasiblity and estimate the running time,costs 3.8 hours on my own Macbook.                                     |
| *gset_expression_origin.R*      | The original R script which cover all high-performance samples, but too much time-consuming.                                                  |
| *gset_expression_pro.R*         | The enhanced version of 'original' R script, in which the matrix transposition part was revised.                                              |
| *gset_expression_ultra.R*       | The bypass version, which directly read the RDS file of gene expression matrix and only perform GSEA process. A time cost of around two days. |


