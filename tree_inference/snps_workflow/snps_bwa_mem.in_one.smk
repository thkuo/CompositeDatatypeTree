# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

# SNP calling with paired-end sequencing reads of DNA
import pandas as pd
from snakemake.utils import validate
from multiprocessing import cpu_count

list_f= config['list_f']
dna_reads= {}
with open(list_f, 'r') as list_fh:
    for l in list_fh:
        d=l.strip().split('\t')
        dna_reads[d[0]]= d[1].split(',')

strains= list(dna_reads.keys())
ref_fasta=config['ref_fasta']
adaptor_f= config['adaptor']
new_reads_dir= config['new_reads_dir']

rule merge_filtered_vcf:
    '''
    Merge the per-strain vcf files 
    '''
    input:
        all_filtered_vcf_gz= lambda wildcards: [
            ('{}/filtered.single_strain_vcf/{}.vcf.gz'
            ).format(wildcards.results_d, strain) for strain in strains],
        all_filtered_vcf_gz_index=lambda wildcards: [
            ('{}/filtered.single_strain_vcf/{}.vcf.gz.tbi'
            ).format(wildcards.results_d, strain) for strain in strains]
    output:
        multisample_filtered_vcf_gz= '{results_d}/merged_vcf/multisample.snp.vcf.gz',
        multisample_filtered_vcf_gz_index= '{results_d}/merged_vcf/multisample.snp.vcf.gz.tbi'
    threads: 8 
    conda: '../shared_envs_yaml/bcftools_env.yml'
    params:
        tabix_bin= 'tabix',
        bgzip_bin= 'bgzip'
    shell:
        """
        bcftools merge --threads {threads} -m none -O z -o \
{output.multisample_filtered_vcf_gz} {input.all_filtered_vcf_gz}
        {params.tabix_bin} -p vcf {output.multisample_filtered_vcf_gz}
        """

rule filter_each_vcf:
    '''
    This block (1) filters for the snps and (2) turn the subset of vcf into snp table
    '''
    input:
        sample_vcf_gz='{results_d}/single_strain_vcf/{strain}.vcf.gz',
        sample_vcf_gz_index='{results_d}/single_strain_vcf/{strain}.vcf.gz.tbi'
    output:
        filtered_sample_vcf_gz='{results_d}/filtered.single_strain_vcf/{strain}.vcf.gz',
        filtered_sample_vcf_gz_index='{results_d}/filtered.single_strain_vcf/{strain}.vcf.gz.tbi'
    threads: 4
    conda: '../shared_envs_yaml/bcftools_env.yml'
    params:
        tabix_bin= 'tabix',
        bgzip_bin= 'bgzip'
    shell:
        """
        ## retain only the snps
        bcftools filter \
  --threads={threads} \
  --include='TYPE="snp" && MIN(FORMAT/DP)>=10 && FORMAT/AO > FORMAT/RO' \
  -O z -o {output.filtered_sample_vcf_gz} {input.sample_vcf_gz}
        {params.tabix_bin} -p vcf {output.filtered_sample_vcf_gz}

        """

rule var_calling:
    input: 
        sorted_bam='{results_d}/{strain}.sorted.bam',
        sorted_bam_index='{results_d}/{strain}.sorted.bam.bai', 
        reffile=ref_fasta,
        reffile_bwa_index=ref_fasta+".fai",
    output:
        raw_vcf_gz= '{results_d}/single_strain_vcf/{strain}.vcf.gz',
        raw_vcf_gz_index= '{results_d}/single_strain_vcf/{strain}.vcf.gz.tbi'
    params:
        tabix_bin= 'tabix',
        bgzip_bin= 'bgzip',
        frbayes_params= '-p 1 --min-alternate-fraction 0.5',
        frbayes_bin= 'freebayes'
    threads: 1
    conda:'../shared_envs_yaml/freebayes_1_3_env.yml'
    shell:
        '''
        freebayes {params.frbayes_params} -f {input.reffile} \
  {input.sorted_bam}| {params.bgzip_bin} -c > {output.raw_vcf_gz} 
        {params.tabix_bin} -p vcf {output.raw_vcf_gz}
        '''

rule mapping:
    '''
    mapping the reads and compute the .sam file
    '''
    input:
        infile1= lambda wildcards: os.path.join(
        new_reads_dir, '{}.cleaned.1.fq.gz'.format(wildcards.strain)),
        infile2= lambda wildcards: os.path.join(
        new_reads_dir, '{}.cleaned.2.fq.gz'.format(wildcards.strain)),
        reffile=ref_fasta,
        ref_index_stampy=ref_fasta+'.stidx',
        ref_index_bwa=ref_fasta+'.bwt',
    output:
        sorted_bam='{results_d}/{strain}.sorted.bam',
        sorted_bam_index='{results_d}/{strain}.sorted.bam.bai'
    threads:8
    params:
        REF_PREFIX=ref_fasta
    conda: '../shared_envs_yaml/snps_tab_mapping.yml'
    shell:
        """
        bwa mem -v 2 -M \
          -t {threads} \
          -R $(echo "@RG\tID:snps\tSM:snps") \
          {input.reffile} {input.infile1} {input.infile2} |
        bamtools sort -out {output.sorted_bam}
        samtools index {output.sorted_bam}
        """

rule redirect_and_preprocess_reads:
    input: 
        infile1=lambda wildcards: dna_reads[wildcards.strain][0],
        infile2=lambda wildcards: dna_reads[wildcards.strain][1]
    output:
        log_f= os.path.join(new_reads_dir, '{strain}.log'),
        f1= os.path.join(new_reads_dir, '{strain}.cleaned.1.fq.gz'),
        f2= os.path.join(new_reads_dir, '{strain}.cleaned.2.fq.gz')
    params:
        adaptor_f= adaptor_f,
        tmp_f1= lambda wildcards: os.path.join(
            new_reads_dir, '{}.cleaned.1.fq'.format(wildcards.strain)),
        tmp_f2= lambda wildcards: os.path.join(
            new_reads_dir, '{}.cleaned.2.fq'.format(wildcards.strain))
    threads: 1
    shell:
        '''
        if [ -e "{params.adaptor_f}" ]
        then
            fastq-mcf -l 50 -q 20 {params.adaptor_f} {input.infile1} {input.infile2} \
  -o {params.tmp_f1} -o {params.tmp_f2} > {output.log_f}
            gzip -9 {params.tmp_f1}
            gzip -9 {params.tmp_f2}
        else
            echo 'Reads not trimmed'
            echo 'No trimming' > {output.log_f}
            echo $(readlink {input.infile1}) >> {output.log_f}
            echo $(readlink {input.infile2}) >> {output.log_f}
            cp {input.infile1} {output.f1}
            cp {input.infile2} {output.f2}
        fi
        '''

rule stampy_index_ref:
    input:
        reffile=ref_fasta
    output:
        ref_fasta+'.fai',
        ref_fasta+'.bwt',
        ref_fasta+'.stidx',
        ref_fasta+'.sthash'
    conda: '../shared_envs_yaml/snps_tab_mapping.yml'
    threads: 1
    shell:
        '''
        stampy.py -G {input.reffile} {input.reffile}
	stampy.py -g {input.reffile} -H {input.reffile}
        bwa index -a bwtsw {input.reffile}
        '''
