# The python scripts which apply a new correlation calculation method #

import numpy as np
import pandas as pd
from ccc.coef import ccc
import os
import math
from scipy import stats

os.chdir('/Users/kimhan/Desktop/CCI_DK/ccc')
# Import the expression matrix
gene_expression_df = pd.read_csv('./gene_expression_matrix.csv', index_col=0)
print(gene_expression_df.shape)
gene_set_expression_df = pd.read_csv('./gene_set_expression_matrix.csv', index_col=0)
print(gene_set_expression_df.shape)
icn_df = pd.read_csv('./icn.csv', index_col=0)

ligand = icn_df.loc[:, 'source_genesymbol'].drop_duplicates().values.tolist()
print(len(ligand))
receptor = icn_df.loc[:, 'target_genesymbol'].drop_duplicates().values.tolist()
print(len(receptor))
ligand_receptor = ligand + receptor

test = [1, 2, 3, 0]


def gm_mean(x, zero_propagate=True):
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


def complex_expression(unit_list):
    exp = gene_expression_df.loc[unit_list, :]
    complex_exp = exp.apply(gm_mean, axis=0).tolist()
    return complex_exp


cor_df = pd.DataFrame(columns= ligand_receptor,index= gene_set_expression_df.index.values)
p_df = pd.DataFrame(columns= ligand_receptor,index= gene_set_expression_df.index.values)
loop = 0
for lr_symbol in ligand_receptor:
    if '_' in lr_symbol:
        unit_list = lr_symbol.split(sep='_')
        if set(unit_list).issubset(set(gene_expression_df.index.values)):
            lr_score = complex_expression(unit_list)
        else:
            loop = loop + len(cor_df.index)
            continue
    else:
        if lr_symbol in gene_expression_df.index:
            lr_score = gene_expression_df.loc[lr_symbol]
        else:
            loop = loop + len(cor_df.index)
            continue
    for gene_set_symbol in gene_set_expression_df.index.values:
        gene_set_score = gene_set_expression_df.loc[gene_set_symbol]
        ccc(lr_score,gene_set_score)