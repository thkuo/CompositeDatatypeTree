# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Tree inference for nucleotide data
import re
from multiprocessing import cpu_count
def parse_strains(vcf_gz):
    import gzip
    vcf= gzip.open(vcf_gz)
    current= '' 
    next_l= vcf.readline().decode('ascii')
    while re.search('^#', next_l):
        current= next_l
        next_l= vcf.readline().decode('ascii')
    strains= [x.strip() for x in current.split('\t')[9:]]
    return(strains)
    
ref_fa=config['ref_fa']
multisample_vcf=config['multisample_vcf']
strains= parse_strains(multisample_vcf)
raxml_model= config['raxml_model']

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
    conda: '../shared_envs_yaml/raxml_env.yml'
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
    conda: '../shared_envs_yaml/raxml_env.yml'
    threads: 8 
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
    conda: '../shared_envs_yaml/raxml_env.yml'
    threads: 16 
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
    params:
        cutoff_num= 2
    run:
        from Bio import SeqIO
        from Bio import AlignIO
        aln_f= input['one_big_aln']
        out_f= output['one_big_var_aln']
        cutoff= int(params['cutoff_num'])
        case= False
        aln= AlignIO.read(aln_f, 'fasta')
        sample_num= len(aln)
        col_num= aln.get_alignment_length()

        new_aln= aln[:, 0:1] ## removed afterward
        new_aln.sort()
        print([s.id for s in new_aln])
        for n in range(col_num):
            col_str= aln[:, n]
            if not case:
                col_str= col_str.lower()
            dominant_char_num= max(
                [col_str.count(uc) for uc in set(col_str)])
            if (sample_num-dominant_char_num) >= cutoff:
                col_aln= aln[:, n:n+1]
                col_aln.sort()
                new_aln=new_aln+col_aln
        new_aln= new_aln[:, 1:]

        with open(out_f, 'w') as out_fh:
            AlignIO.write(new_aln, out_fh, "fasta")

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
        per_sample_seq= '{result_dir}/alignment/strains/{sample}.fa'
    conda: '../shared_envs_yaml/bcftools_env.yml'
    threads: 1
    shell:
        '''
        bcftools consensus \
  --fasta-ref {input.ref_fa} \
  --sample {wildcards.sample}  {input.multisample_snps_vcf}| 
  sed -e "/^>/s/^>.\+/>{wildcards.sample}/" > {output.per_sample_seq}
        '''
    
