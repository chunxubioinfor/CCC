# The python scripts which apply a new correlation calculation method #

import numpy as np
import pandas as pd
from ccc.coef import ccc
import os
import math
from scipy import stats

os.chdir('/Users/kimhan/Desktop/CCI_DK/ccc')
## Import the expression matrix
gene_expression_df = pd.read_csv('./gene_expression_matrix.csv', index_col=0)
print(gene_expression_df.shape)
gene_set_expression_df = pd.read_csv('./gene_set_expression_matrix.csv', index_col=0)
print(gene_set_expression_df.shape)
icn_df = pd.read_csv('./icn.csv', index_col=0)

ligand = icn_df.loc[:, 'source_genesymbol'].drop_duplicates().values.tolist()
print(len(ligand))
receptor = icn_df.loc[:, 'target_genesymbol'].drop_duplicates().values.tolist()
print(len(receptor))

test = [1, 2, 3, 0]
def gm_mean(x, zero_propagate=False):
    if (np.array(x) < 0).any():
        return None
    if zero_propagate is True:
        if 0 in x:
            return 0
        g_mean = math.exp(np.mean(np.log(x)))
        return g_mean
    else:
        g_mean = math.exp(sum(np.log([i for i in x if i > 0]))/len(x))
        return g_mean


