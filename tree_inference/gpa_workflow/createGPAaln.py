import pandas as pd
import textwrap

gpa_f= snakemake.input['roary_raw_gpa']
output_aln_f= snakemake.output['gpa_aln']

## read roary result and identify single-copy genes
df=pd.read_csv(gpa_f, sep= ',',
        header= 0, index_col= 0, quotechar= '"', low_memory=False)

## filter and convert the states
print('read gene orthologous families')
sub_df=df.iloc[:, 13:df.shape[1]]
sub_df= sub_df.applymap(lambda x: '1' if (not pd.isna(x)) else '0')
sub_df= sub_df.transpose() # strains in rows

## print 
print('create fasta of gpa alignment')
out_fh=open(output_aln_f, 'w')
for s in sub_df.index.values.tolist():
    out_fh.write('>'+s+'\n') # fasta header
    seq_str='\n'.join(textwrap.wrap(''.join(sub_df.loc[s,:]), 
        width=60))+'\n'
    out_fh.write(seq_str)
out_fh.close()
