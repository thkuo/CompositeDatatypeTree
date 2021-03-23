nuc_tr= config['nuc_tr']
nuc_aln= config['nuc_aln']
raxml_model= config['raxml_model']
col_nuc_tr= config['col_nuc_tr']
cutoff_perc= config['cutoff_perc']
rule all: 
    input:
        col_nuc_tr

rule collapse_nuc_tr_with_psll_sums:
    input:
        nuc_pslls= '{result_dir}/psll/raxml/RAxML_perSiteLLs.nuc_PSLL',
        rear_trs= '{result_dir}/psll/rear_trs.nwk'
    output:
        nuc_psll_sum_per_br='{result_dir}/psll/nuc_psllChangeSumPerBranch.tsv',
        col_nuc_tr='{result_dir}/col_nuc.nwk' 
    params:
        cutoff_perc= cutoff_perc
    threads: 16
    script: 'scripts/psll_analysis.R'

rule compute_nni_trees:
    input:
        nuc_tr=nuc_tr
    output:
        rear_trs= '{result_dir}/psll/rear_trs.nwk'
    script: 'scripts/make_nni_trees.R'
        
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
