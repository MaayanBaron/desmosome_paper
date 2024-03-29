import warnings
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
import scipy

warnings.filterwarnings('ignore')

#bootstrap functions
def get_sample(df_diff):
    x = np.random.choice(df_diff.values,size=df_diff.size,replace=False)
    return x

def bootstrap_stat(df_diff,coh,num_of_bs):
    df_scramble = pd.DataFrame(columns=['mean','median'],index=range(0,num_of_bs))
    x = [get_sample(df_diff.iloc[:,0]) for i in range(num_of_bs)]
    x = pd.DataFrame(x).T
    x.index =  df_diff.index
    df_scramble['mean'] = x.loc[des_genes].mean()
    df_scramble['median'] = x.loc[des_genes].median()
    return df_scramble

import os
os.chdir('/cellar/users/m2baron/work2020/MEL_desmosomes/')

#define des genes
des_genes = ['DSP','DSC1','DSC2','DSC3','DSG1','DSG2',
               'DSG3','DSG4','PKP1','PKP2', 'PKP3','PKP4','JUP']

#calculate p-value for each cohort based on MutSig output
#1. import MutSig expected mutation number
path = '/cellar/users/m2baron/work2020/MEL_desmosomes/mutsigcv1.4_pancanceratlas_expect.csv'
df_expected = pd.read_csv(path,index_col=0,sep='\t')
#fix columns name to match
new_col = []
for col in df_expected.columns:
    new_col.append(col.split('_')[0])#+'_tcga')
df_expected.columns = new_col
#2. import observed TCGA mutation number
path = '/cellar/users/m2baron/work2020/MEL_desmosomes/mutsigcv1.4_pancanceratlas_observe.csv'
df_observed = pd.read_csv(path,index_col=0,sep='\t')
#fix columns name to match
new_col = []
for col in df_observed.columns:
    new_col.append(col.split('_')[0])#+'_tcga')
df_observed.columns = new_col
#3. calculate the diffrence and sort genes
df_diff = df_observed-df_expected

df_pval = pd.DataFrame(index=df_diff.columns,columns=['pval'])
for coh in df_diff.columns:
    print(coh)
    df_scramble = bootstrap_stat(df_diff=df_diff,coh=coh,num_of_bs=1000)
    df_pval.loc[coh,'pval'] = (np.mean(df_diff[coh].loc[des_genes])>df_scramble['median']).mean()

#adding conway datatset for melanoma:
path = '/cellar/users/m2baron/work2020/MEL_desmosomes/Conway_et_al/Conway_expected_mutsig.txt'
df_conway_exp = pd.read_csv(path, sep='\t',header=None,index_col=0)
path = '/cellar/users/m2baron/work2020/MEL_desmosomes/Conway_et_al/mutsigCV_output.txt'
df_conway_obs = pd.read_csv(path, sep='\t',index_col=0)

diff = df_conway_obs['n_nonsilent']-df_conway_exp.iloc[:,0]
df_conway_diff = pd.DataFrame(data=diff,columns=['conway'])

df_scramble = bootstrap_stat(df_diff=df_conway_diff,coh='conway',num_of_bs=1000)
df_pval.loc['conway','pval'] = (np.mean(df_conway_diff['conway'].loc[des_genes])<df_scramble['median']).mean()

#adding chang datatset for melanoma:
path = '/cellar/users/m2baron/work2020/MEL_desmosomes/CutaneousSquamous/chang_expected_mutsig.txt'
df_chang_exp = pd.read_csv(path, sep='\t',header=None,index_col=0)
path = '/cellar/users/m2baron/work2020/MEL_desmosomes/CutaneousSquamous/mutsig_chang.txt.sig_genes.txt'
df_chang_obs = pd.read_csv(path, sep='\t',index_col=0)

diff = df_chang_obs['n_nonsilent']-df_chang_exp.iloc[:,0]
df_chang_diff = pd.DataFrame(data=diff,columns=['chang'])

df_scramble = bootstrap_stat(df_diff=df_chang_diff,coh='chang',num_of_bs=1000)
df_pval.loc['chang','pval'] = (np.mean(df_chang_diff['chang'].loc[des_genes])<df_scramble['median']).mean()

df_pval.to_csv('bootstrap_pvals_100k.tsv',sep='\t')
