'''
__purpose__=Tree inference for nucleotide data
__author__=Tzu-Hao Kuo
__description__=from vcf file of multiple samples to nucleotide tree
'''
import re
from multiprocessing import cpu_count
ref_fa=config['ref_fa']
multisample_vcf=config['multisample_vcf']
raxml_model= config['raxml_model']
strains= config['strains'].strip().split(',')
bootstrap_cores=5
rule bs_values_mapped_to_tree:
    input:  
        cw_all_bs_trees='{result_dir}/raxml/bootstrap/RAxML_bootstrap.{suffix}.all',
        nuc_best_tr='{result_dir}/raxml/RAxML_bestTree.{suffix}'
    output:
        cw_withBS_tree= '{result_dir}/raxml/RAxML_bipartitions.{suffix}.bs'
    params:
        raxml_wd= '{result_dir}/raxml',
        tree_full_suffix='{suffix}.bs',
        raxml_model= raxml_model
    threads:1
    conda: '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
        fi
        $RAXML_BIN \
-f b \
-t {input.nuc_best_tr} \
-z {input.cw_all_bs_trees} \
-w {params.raxml_wd} \
-m {params.raxml_model} \
-n {params.tree_full_suffix}
        """ 

rule collect_bs_trees:
    input:
        bs_trees=lambda wildcards:[
        os.path.join(wildcards.result_dir, 'raxml','bootstrap',
        'RAxML_bootstrap.{}.{}.20'.format(wildcards.suffix, part_num)) for part_num in range(1, 5+1)]
    output:
        cw_all_bs_trees='{result_dir}/raxml/bootstrap/RAxML_bootstrap.{suffix}.all'
    shell:
        """
        cat {input.bs_trees} > {output.cw_all_bs_trees}
        """

rule bootstrap:
    input:
        one_big_var_aln= '{result_dir}/alignment/{suffix}.var.aln'
    output:
        bs_tree='{result_dir}/raxml/bootstrap/RAxML_bootstrap.{suffix}.{part_num}.{per_part_tr_num}'
    params:
        bs_dir= '{result_dir}/raxml/bootstrap',
        bs_num= 20,
        tree_id= '{suffix}.{part_num}.{per_part_tr_num}',
        raxml_model= raxml_model,
        raxml_starting_num= 1,
    conda: '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    threads: bootstrap_cores 
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
        fi
        $RAXML_BIN \
  -T {threads} -N {params.bs_num}\
  -m {params.raxml_model} \
  -s {input.one_big_var_aln} \
  -n {params.tree_id} \
  -w {params.bs_dir} \
  -# {wildcards.per_part_tr_num} \
  -p {wildcards.part_num} -x {wildcards.part_num}
        """

rule nuc_best_tree:
    input:
        one_big_var_aln= '{result_dir}/alignment/{suffix}.var.aln'
    output:
        nuc_best_tr='{result_dir}/raxml/RAxML_bestTree.{suffix}'
    params:
        tr_id='{suffix}',
        raxml_seed= '1',
        raxml_starting_num= 1,
        raxml_model= raxml_model,
        raxml_wd=lambda wildcards: os.path.join(wildcards.result_dir, 'raxml')
    conda: '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    threads:
        lambda cores: max(1, cpu_count() - 1)
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
        fi
        $RAXML_BIN \
  -T {threads} \
  -m {params.raxml_model} \
  -s {input.one_big_var_aln} \
  -n {params.tr_id} \
  -# {params.raxml_starting_num} \
  -w {params.raxml_wd} \
  -p {params.raxml_seed}
        """

rule mapping_trim_invariant:
    '''
    Trim the conserved sites and retain those where 
    the dominant state comprises less than (the number of samples - 1)
    '''
    input:
        one_big_aln= '{result_dir}/alignment/{suffix}.full.aln'
    output:
        one_big_var_aln= '{result_dir}/alignment/{suffix}.var.aln'
    threads: 1
    shell:
        '''
        ./removeInvariant.py --in {input} \
--out {output} --cn 2
        '''

rule cons_seqs_to_aln:
    '''
    Create the nucleotide alignment 
    '''
    input:
        per_sample_seq= expand('{result_dir}/alignment/strains/{sample}.fa', 
            result_dir= '{result_dir}', sample= strains)
    output:
        one_big_aln= '{result_dir}/alignment/{suffix}.full.aln'
    threads: 1
    shell:
        '''
        cat {input.per_sample_seq} > {output.one_big_aln}
        '''

rule vcf_to_seq:
    '''
    Introduce the alternative alleles to the reference
    '''
    input:
        multisample_snps_vcf= multisample_vcf,
        multisample_snps_vcf_ix= multisample_vcf+'.tbi',
        ref_fa=ref_fa
    output:
        strain_cons_seq= temp('{result_dir}/alignment/strains/{sample}.fa')
    conda: '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/bcftools_env.yml'
    threads: 1
    shell:
        '''
        bcftools consensus \
  --fasta-ref {input.ref_fa} \
  --sample {wildcards.sample}  {input.multisample_snps_vcf}| 
  sed -e "/^>/s/^>.\+/>{wildcards.sample}/" > {output.strain_cons_seq}
        '''
    
