import os

gpa_aln= config['gpa_aln']
nuc_aln= config['nuc_aln']
conc_aln= config['conc_aln']
raxml_prtn= config['raxml_prtn']
col_constraint_tr= config['col_constraint_tr']
#best_ml_tree= config['best_ml_tree'] 
#bs_trees= config['bs_trees']
#bs_best_ml_tree= config['bs_best_ml_tree']

#cw_withBS_tree= ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
#'TreePaper/nuc_tree_comparison.v2/results/mapping/'
#'conc.raxml/RAxML_bestTree.conc')
#col_cw_withBS_tree= os.path.join(results_dir,
#    'RAxML_bipartitions.nuc_uw-gpa.withBS.collapse')
#max_core_n= config['max_cores']
max_core_n=len(os.sched_getaffinity(0)) 
max_per_part_core_n= int(max_core_n/5)
#raxml_bin= config['raxml_bin'] 
model= config['model']
#rule all:
#    input:
#        bs_best_ml_tree,
#        best_ml_tree

#rule collapse_conc_tree:
#    input: cw_withBS_tree
#    output: col_cw_withBS_tree
#    shell:
#        '''
#        ./collapse_short.interface.R \
# {input} {output}
#        '''

rule for_cw_map_bs_values_to_tree:
    input:  
        cw_all_bs_trees='{result_dir}/bootstrap/RAxML_bootstrap.all',
        conc_best_tr='{result_dir}/RAxML_bestTree.{suffix}'
    output:
        cw_withBS_tree= '{result_dir}/RAxML_bipartitions.{suffix}.bs'
    params:
#        raxml_bin=raxml_bin,
        tree_full_suffix='{suffix}.bs',
        raxml_model= model
    threads:1
    conda:'/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
        fi
        $RAXML_BIN \
-f b \
-t {input.conc_best_tr} \
-z {input.cw_all_bs_trees} \
-w {wildcards.result_dir} \
-m {params.raxml_model} \
-n {params.tree_full_suffix}
        """ 

rule for_cw_collect_bs_trees:
    input:
        bs_trees=lambda wildcards:[
        os.path.join(wildcards.result_dir,'bootstrap',
        'RAxML_bootstrap.{}.20'.format(part_num)) for part_num in range(1, 5+1)]
    output:
        cw_all_bs_trees='{result_dir}/bootstrap/RAxML_bootstrap.all'
    shell:
        """
        cat {input.bs_trees} > {output.cw_all_bs_trees}
        """

rule bootstrap:
    input: 
        tr= col_constraint_tr,
        raxml_input_aln=conc_aln,
        raxml_input_partitions=raxml_prtn
    output:
        bs_tree='{result_dir}/bootstrap/RAxML_bootstrap.{part_num}.{per_part_tr_num}'
    params:
        bs_dir= '{result_dir}/bootstrap',
#        raxml_bin=raxml_bin,
        raxml_model= model,
        raxml_starting_num= 1,
    conda:'/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    threads: max_per_part_core_n
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'     
        fi
        $RAXML_BIN \
 -T {threads} \
 -g {input.tr} \
 -q {input.raxml_input_partitions} \
 -m {params.raxml_model} \
 -s {input.raxml_input_aln} \
 -n {wildcards.part_num}.{wildcards.per_part_tr_num} \
 -p {wildcards.part_num} -x {wildcards.part_num} \
 -# {wildcards.per_part_tr_num} \
 -w {params.bs_dir}
        """

rule compute_conc_tree:
    input:
        tr= col_constraint_tr,
        raxml_input_aln=conc_aln,
        raxml_input_partitions=raxml_prtn
    output:
        conc_best_tr='{result_dir}/RAxML_bestTree.{suffix}'
    params:
#        raxml_bin=raxml_bin,
        raxml_model= model,
        raxml_starting_num= 1
    threads: max_core_n
    conda:'/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'     
        fi
        $RAXML_BIN \
  -T {threads} \
  -g {input.tr} \
  -q {input.raxml_input_partitions} \
  -m {params.raxml_model} \
  -s {input.raxml_input_aln} \
  -p 1 \
  -# {params.raxml_starting_num} \
  -w {wildcards.result_dir} \
  -n {wildcards.suffix}
        """

rule prepare_concatenated_data:
    input:
        nuc_aln= nuc_aln,
        gpa_aln= gpa_aln
    output:
        raxml_input_aln=conc_aln,
        raxml_input_partitions=raxml_prtn
    params:
        conc_data_script= './create_raxml_conc_inputs.py'
    shell:
        '''
        {params.conc_data_script} --nuc_aln {input.nuc_aln} \
 --gpa_aln {input.gpa_aln} --out_aln {output.raxml_input_aln} \
 --out_prtn {output.raxml_input_partitions}
        '''
        
    
