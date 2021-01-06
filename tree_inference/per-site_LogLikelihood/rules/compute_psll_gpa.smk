
rule compute_psll_gpa:
    input:
        gpa_aln= gpa_aln,
        rear_trs= '{result_dir}/psll/rear_trs.nwk'
    output:
        gpa_pslls= '{result_dir}/psll/raxml/RAxML_perSiteLLs.gpa_PSLL'
    params: 
        raxml_id='gpa_PSLL',
        raxml_wd= '{result_dir}/psll/raxml',
        raxml_model= 'BIN{}'.format(raxml_d_model)
    threads: 16
    conda: '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
        fi
        $RAXML_BIN -f g \
  -T {threads} -s {input.gpa_aln} \
  -m {params.raxml_model} \
  -z {input.rear_trs} \
  -n {params.raxml_id} \
  -w {params.raxml_wd} 
        """
