# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Tree inference for nucleotide data
import re
import os
from multiprocessing import cpu_count
import gzip


def parse_strains(vcf_gz):
    vcf= gzip.open(vcf_gz)
    current= ''
    next_l= vcf.readline().decode('ascii')
    while re.search('^#', next_l):
        current= next_l
        next_l= vcf.readline().decode('ascii')
    strains= [x.strip() for x in current.split('\t')[9:]]
    return(strains)


ref_fa=config['ref_fa']
if os.path.isfile(config['multisample_vcf']):
    # if the vcf files have been merged
    multisample_vcf=config['multisample_vcf']
    strains= parse_strains(multisample_vcf)
else:
    # if the vcf files are put in a folder
    strain_vcf_dir = os.path.join(
        os.path.dirname(os.path.dirname(config['multisample_vcf'])),
        'filtered.single_strain_vcf')
    strains = [re.sub('.vcf.gz', '', f)
        for f in os.listdir(strain_vcf_dir)
        if re.search('.vcf.gz$', f) is not None]
    assert len(strains) == len(set(strains)), 'Repeated strain names'
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
    priority: 1
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

rule merge_filtered_aln:
    input:
        expand('{{result_dir}}/alignment/{{suffix}}.full.aln_{p}-filtered',
                p=list(range(10)))
    output:
        one_big_var_aln= '{result_dir}/alignment/{suffix}.var.aln'
    run:
        from Bio import AlignIO
        aln = AlignIO.read(input[0], 'fasta')
        aln.sort()
        for f in input[1:]:
            sub_filtered_aln = AlignIO.read(f, 'fasta')
            sub_filtered_aln.sort()
            aln += sub_filtered_aln
        with open(output.one_big_var_aln, 'w') as out_fh:
            AlignIO.write(aln, out_fh, 'fasta')


rule sub_filtered_aln:
    input:
        sub_aln_f='{result_dir}/alignment/{suffix}.full.aln_{p}'
    output:
        sub_filtered_aln_f=temp('{result_dir}/alignment/{suffix}.full.aln_{p}-filtered')
    params:
        cutoff_num=2,
        case_sensitive=False
    threads: 1
    run:
        from Bio import AlignIO
        from tqdm import trange
        aln = AlignIO.read(input.sub_aln_f, 'fasta')
        cutoff = int(params.cutoff_num)
        case = params.case_sensitive
        sample_num = len(aln)
        col_num = aln.get_alignment_length()

        new_aln = aln[:, 0:1]  # removed afterward
        new_aln.sort()
        for n in trange(col_num):
            col_str = aln[:, n]
            if not case:
                col_str = col_str.lower()
            dominant_char_num = max(
                [col_str.count(uc) for uc in set(col_str)])
            if (sample_num - dominant_char_num) >= cutoff:
                col_aln = aln[:, n:n+1]
                col_aln.sort()
                new_aln = new_aln+col_aln

        new_aln = new_aln[:, 1:]
        with open(output.sub_filtered_aln_f, 'w') as sub_filtered_aln_fh:
            AlignIO.write(new_aln, sub_filtered_aln_fh, 'fasta')


rule divide_aln:
    input:
        one_big_aln= '{result_dir}/alignment/{suffix}.full.aln'
    output:
        sub_aln_files=temp(expand(
            '{{result_dir}}/alignment/{{suffix}}.full.aln_{p}',
             p=list(range(10))))
    run:
        import math
        from Bio import AlignIO
        aln = AlignIO.read(input.one_big_aln, 'fasta')
        print(aln.get_alignment_length())
        num_partitions = 10
        chunk_size = float(aln.get_alignment_length()) / num_partitions

        processed_size = 0
        for n in range(num_partitions):
            sub_aln = aln[:,
                          math.floor(n*chunk_size):
                          math.floor((n+1)*chunk_size)]
            sub_aln_f = '{}_{}'.format(input.one_big_aln, str(n))
            with open(sub_aln_f, 'w') as sub_aln_fh:
                AlignIO.write(sub_aln, sub_aln_fh, 'fasta')
            processed_size += sub_aln.get_alignment_length()
        assert processed_size == aln.get_alignment_length()


rule cons_seqs_to_aln:
    # Create the nucleotide alignment 
    input:
        per_sample_seq= expand('{result_dir}/alignment/strains/{sample}.fa', 
            result_dir= '{result_dir}', sample= strains)
    output:
        one_big_aln= '{result_dir}/alignment/{suffix}.full.aln'
    threads: 1
    run:
        from Bio import SeqIO
        seq_records = [SeqIO.read(f, 'fasta') for f in input.per_sample_seq]
        with open(output.one_big_aln, 'w') as out_fh:
            SeqIO.write(seq_records, out_fh, 'fasta')

rule vcf_to_seq:
    # Introduce the alternative alleles to the reference
    input:
        strain_snps_vcf= '{result_dir}/filtered.single_strain_vcf/{sample}.vcf.gz',
        strain_snps_vcf_ix= '{result_dir}/filtered.single_strain_vcf/{sample}.vcf.gz.tbi',
        ref_fa=ref_fa
    output:
        per_sample_seq= '{result_dir}/alignment/strains/{sample}.fa'
    conda: '../shared_envs_yaml/bcftools_env.yml'
    threads: 1
    shell:
        '''
        bcftools consensus \
  --fasta-ref {input.ref_fa} \
  --sample {wildcards.sample}  {input.strain_snps_vcf}|
  sed -e "/^>/s/^>.\+/>{wildcards.sample}/" > {output.per_sample_seq}
        '''
