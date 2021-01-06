import re

nuc_tr= config['nuc_tr']
nuc_aln= config['nuc_aln']
gpa_aln= config['gpa_aln'] if 'gpa_aln' in config else '-'
raxml_s_model= config['raxml_substitution_model']
raxml_d_model= config['raxml_distribution_model']
# estimate the invariant proportion or not
raxml_i_model= ''
if 'raxml_invariant' in config :
    assert config['raxml_invariant'] in ['', 'I'] 
    raxml_i_model= config['raxml_invariant'] 
# estimate the base frequencies or not
raxml_f_model= ''
if 'raxml_freq' in config :
    assert config['raxml_freq'] in ['', 'X'] 
    raxml_i_model= config['raxml_freq'] 
raxml_model= raxml_s_model+raxml_d_model+raxml_i_model+raxml_f_model
#' asc correction
if not (re.search('ASC_', raxml_model) is None) :
    raxml_model= raxml_model+ ' --asc-corr lewis'
col_nuc_tr= config['col_nuc_tr']
#' the cutoff of log-likelihood score
col_cutoff_pr= 0.75
if 'col_cutoff_pr' in config:
    assert config['col_cutoff_pr'] > 0 & config['col_cutoff_pr'] < 1
    col_cutoff_pr= config['col_cutoff_pr']

#' external rules
include: 'compute_psll_gpa.smk'

rule all: 
    input:
        col_nuc_tr

rule collapse_nuc_tr_with_psll_sums:
    input:
        nuc_pslls= '{result_dir}/psll/raxml/RAxML_perSiteLLs.nuc_PSLL',
        rear_trs= '{result_dir}/psll/rear_trs.nwk'
    output:
        nuc_psll_sum_per_br='{result_dir}/psll/nuc_psllChangeSumPerBranch.tsv',
        col_nuc_tr='{result_dir}/nuc_col_by_psll.nwk' 
    params: 
        col_cutoff_pr= col_cutoff_pr
    threads: 16
    script: '../scripts/psll_analysis.R'

rule compute_nni_trees:
    input:
        nuc_tr=nuc_tr
    output:
        rear_trs= '{result_dir}/psll/rear_trs.nwk'
    script: '../scripts/make_nni_trees.R'
        
rule compute_psll:
    input:
        nuc_aln= nuc_aln,
        rear_trs= '{result_dir}/psll/rear_trs.nwk'
    output:
        nuc_pslls= '{result_dir}/psll/raxml/RAxML_perSiteLLs.nuc_PSLL'
    params: 
        raxml_id='nuc_PSLL',
        raxml_wd= '{result_dir}/psll/raxml',
        raxml_model= raxml_model
    threads: 16
    conda: '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
        fi
        $RAXML_BIN -f g \
  -T {threads} -s {input.nuc_aln} \
  -m {params.raxml_model} \
  -z {input.rear_trs} \
  -n {params.raxml_id} \
  -w {params.raxml_wd} 
        """
