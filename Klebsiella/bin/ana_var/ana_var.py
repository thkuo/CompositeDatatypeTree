# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This script uses principal component analysis (PCA)
# to detect highly divergent samples among the
# 100 Klebsiella pneumoniae samples

import gzip
import os
import re
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from adjustText import adjust_text

# load the vcf
vcf_gz_f='../../results/multisample.vcf.gz'
vcf_lines= []
with gzip.open(vcf_gz_f, mode= 'rb') as fh:
    vcf_lines= [l.decode('utf-8').strip().split('\t') for l in fh.readlines() if (not re.match('##', l.decode('utf-8')))]
    vcf_df= pd.DataFrame(vcf_lines, dtype= str)
    vcf_df.rename(columns= vcf_df.iloc[0], inplace= True)
    # only retain the variant data per strain
    vcf_df= vcf_df.iloc[1:,9:]

# convert the variants to binary states
vcf_bin_df= vcf_df.applymap(lambda x: 0 if re.search('\w', x) is None else 1)
print(vcf_df.iloc[:3, :3])
print(vcf_bin_df.iloc[:3, :3])
# per column
assert all(vcf_bin_df.sum(axis= 1).apply(lambda s: s>0))
print(vcf_bin_df.sum(axis= 1).quantile([0, .25, .5, .75, 1.0]))
# per row
assert all(vcf_bin_df.sum(axis= 0).apply(lambda s: s>0))
print(vcf_bin_df.sum(axis= 0).quantile([0, .25, .5, .75, 1.0]))
print('Binarization done')

pca= PCA(n_components= 2)
vcf_bin_mat= vcf_bin_df.T.to_numpy()
samples= vcf_bin_df.columns.values.tolist()
# remember to put the samples in rows
vcf_bin_transformed_mat=pca.fit_transform(vcf_bin_mat)
plt.scatter(vcf_bin_transformed_mat[:, 0],  
            vcf_bin_transformed_mat[:, 1],
            facecolor="#0D65D9", lw=0, alpha= 0.2)
plt.axis('scaled' )
plt.style.use('classic')
texts = []
sample_filter= vcf_bin_df.sum(axis= 0).apply(lambda x: x<=1000)
for i, txt in enumerate(samples):
    if not sample_filter[txt]:
        # unusually many SNPs
        texts.append(plt.text(vcf_bin_transformed_mat[i, 0],
                              vcf_bin_transformed_mat[i, 1],
                              txt))
adjust_text(texts, only_move={'texts':'xy'}, 
           arrowprops=dict(arrowstyle='-', color='grey'))
plt.title('PCA of SNPs presence or absence')
plt.xlabel('PCA1')
plt.ylabel('PCA2')

fig_name='../../results/ana_var/SNPs_pca.pdf'
if not os.path.isdir(os.path.dirname(fig_name)):
    os.mkdir(os.path.dirname(fig_name))
plt.savefig(fig_name)

