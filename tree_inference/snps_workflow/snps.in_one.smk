'''
__purpose__=SNP calling with paired-end sequencing reads of DNA
__author__=Tzu-Hao Kuo
__description__=adjusted from seq2geno
'''
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
#snps_table=config['snps_table']
adaptor_f= config['adaptor']
new_reads_dir= config['new_reads_dir']

rule filter_vcf:
    '''
    This rule (1) filters for the snps and (2) turn the subset of vcf into snp table
    '''
    input:
        multisample_raw_vcf_gz= '{results_d}/merged_vcf/multisample.vcf.gz',
        multisample_raw_vcf_gz_index= '{results_d}/merged_vcf/multisample.vcf.gz.tbi'
    output:
        multisample_snp_vcf= temp('{results_d}/merged_vcf/multisample.snp.vcf'),
        multisample_snp_vcfgz= '{results_d}/merged_vcf/multisample.snp.vcf.gz',
        multisample_snp_vcfgz_index= '{results_d}/merged_vcf/multisample.snp.vcf.gz.tbi'
    threads: 8
    conda: '../shared_envs_yaml/bcftools_env.yml'
    params:
        filter_pat= 'TYPE="snp" && FORMAT/DP>=10 && FORMAT/AO > FORMAT/RO',
        tabix_bin= 'tabix',
        bgzip_bin= 'bgzip'
    shell:
        """
        ## retain only the snps
        bcftools norm --multiallelics=-any --check-ref=x --threads={threads} \
  {input.multisample_raw_vcf_gz}| \
  bcftools filter --include='{params.filter_pat}' \
  --threads={threads} -O v \
  -o  {output.multisample_snp_vcf}
        {params.bgzip_bin} -c {output.multisample_snp_vcf} > {output.multisample_snp_vcfgz}
        {params.tabix_bin} -p vcf {output.multisample_snp_vcfgz}
        """

rule merge_vcf:
    '''
    Merge the per-strain vcf files 
    '''
    input:
        all_raw_vcf_gz= lambda wildcards: [
            ('{}/single_strain_vcf/{}.vcf.gz'
            ).format(wildcards.results_d, strain) for strain in strains],
        all_raw_vcf_gz_index=lambda wildcards: [
            ('{}/single_strain_vcf/{}.vcf.gz.tbi'
            ).format(wildcards.results_d, strain) for strain in strains]
    output:
        multisample_raw_vcf_gz= '{results_d}/merged_vcf/multisample.vcf.gz',
        multisample_raw_vcf_gz_index= '{results_d}/merged_vcf/multisample.vcf.gz.tbi'
    threads: 8
    conda: '../shared_envs_yaml/bcftools_env.yml'
    params:
        tabix_bin= 'tabix',
        bgzip_bin= 'bgzip'
    shell:
        """
        bcftools merge --threads {threads} -m none -O z -o \
{output.multisample_raw_vcf_gz} {input.all_raw_vcf_gz}
        {params.tabix_bin} -p vcf {output.multisample_raw_vcf_gz}
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

rule sam2bam:
    input:
        sam='{results_d}/{strain}.sam',
    output:
        bam=temp('{results_d}/{strain}.bam'),
        sorted_bam='{results_d}/{strain}.sorted.bam',
        sorted_bam_index='{results_d}/{strain}.sorted.bam.bai'
    threads: 1
    conda: '../shared_envs_yaml/snps_tab_mapping.yml'
    shell:
        """
        # sam to bam 
        samtools view -@ {threads} -bS {input.sam} > {output.bam}

        # sort
        bamtools sort -in {output.bam} -out {output.sorted_bam}
        samtools index {output.sorted_bam}
        """

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
        p_bwa_sai1=temp('{results_d}/{strain}.bwa1.sai'),
        p_bwa_sai2=temp('{results_d}/{strain}.bwa2.sai'),
        bwa_sam= temp('{results_d}/{strain}.bwa.sai'),
        bwa_bam= temp('{results_d}/{strain}.pre_bwa.bam'),
        sam=temp('{results_d}/{strain}.sam')
    threads:1
    params:
        REF_PREFIX=ref_fasta,
        BWA_OPT='-q10'
    conda: '../shared_envs_yaml/snps_tab_mapping.yml'
    shell:
        """
        bwa aln {params.BWA_OPT} -t{threads} {input.reffile} {input.infile1} \
  > {output.p_bwa_sai1}

        bwa aln {params.BWA_OPT} -t{threads} {input.reffile} {input.infile2} \
  > {output.p_bwa_sai2} 

        bwa sampe -r '@RG\\tID:{wildcards.strain}\\tSM:{wildcards.strain}' \
  {input.reffile} {output.p_bwa_sai1} {output.p_bwa_sai2} \
  {input.infile1} {input.infile2} > {output.bwa_sam}
        samtools view -@ {threads} -Sb {output.bwa_sam} >  {output.bwa_bam} 

        ## ng_stampy_remapping:
        stampy.py \
  --readgroup=ID:{wildcards.strain},SM:{wildcards.strain}\
  -g {params.REF_PREFIX} -h {params.REF_PREFIX} \
  -t{threads}  --bamkeepgoodreads -M  {output.bwa_bam}\
  > {output.sam}
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
