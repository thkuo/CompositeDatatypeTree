nuc_tr= config['nuc_tr']
nuc_aln= config['nuc_aln']
raxml_model= config['raxml_model']
col_nuc_tr= config['col_nuc_tr']
cutoff_perc= config['cutoff_perc']
rule all: 
    input:
        col_nuc_tr

rule collapse_nuc_tr_with_ll_sums:
    input:
        nuc_lls= '{result_dir}/ll/raxml/RAxML_info.nuc_LL',
        rear_trs= '{result_dir}/ll/rear_trs.nwk'
    output:
        nuc_ll_sum_per_br='{result_dir}/ll/nuc_llChangeSumPerBranch.tsv',
        col_nuc_tr='{result_dir}/nuc_col_by_ll.nwk' 
    params:
        cutoff_perc= cutoff_perc
    threads: 16
    script: 'scripts/ll_analysis.R'

rule compute_nni_trees:
    input:
        nuc_tr=nuc_tr
    output:
        rear_trs= '{result_dir}/ll/rear_trs.nwk'
    script: 'scripts/make_nni_trees.R'
        
rule compute_ll:
    input:
        nuc_aln= nuc_aln,
        rear_trs= '{result_dir}/ll/rear_trs.nwk'
    output:
        nuc_lls= '{result_dir}/ll/raxml/RAxML_info.nuc_LL'
    params: 
        raxml_id='nuc_LL',
        raxml_wd= '{result_dir}/ll/raxml',
        raxml_model= raxml_model
    threads: 16
    conda: '../shared_envs_yaml/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
        fi
        $RAXML_BIN -f n \
  -T {threads} -s {input.nuc_aln} \
  -m {params.raxml_model} \
  -z {input.rear_trs} \
  -n {params.raxml_id} \
  -w {params.raxml_wd} 
        """
