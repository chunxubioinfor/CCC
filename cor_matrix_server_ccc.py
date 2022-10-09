# The python scripts which apply a new correlation calculation method #

import numpy as np
import pandas as pd
from ccc.coef import ccc
import os

os.chdir('/Users/kimhan/Desktop/CCI_DK/ccc')
## Import the expression matrix
gene_expression_df = pd.read_csv('./gene_expression_matrix.csv', index_col=0)
print(gene_expression_df.shape)
gene_set_expression_df = pd.read_csv('./gene_set_expression_matrix.csv', index_col=0)
print(gene_set_expression_df.shape)
icn_df = pd.read_csv('./icn.csv',index_col=0)

ligand = icn_df.loc[:,'source_genesymbol'].drop_duplicates().values.tolist()
print(len(ligand))
receptor = icn_df.loc[:,'target_genesymbol'].drop_duplicates().values.tolist()
print(len(receptor))

