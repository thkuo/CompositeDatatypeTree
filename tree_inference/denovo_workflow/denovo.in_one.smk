import os
import pandas as pd

list_f=config['list_f']
dna_reads= {}
with open(list_f, 'r') as list_fh:
    for l in list_fh:
        d=l.strip().split('\t')
        dna_reads[d[0]]= d[1].split(',')
out_prokka_dir=config['out_prokka_dir']
out_roary_dir=config['out_roary_dir']
out_spades_dir= config['out_spades_dir']
adaptor_f= config['adaptor']
new_reads_dir= config['new_reads_dir']

rule roary:
    input:
        gff_files= expand(os.path.join(out_prokka_dir, '{strain}', '{strain}.gff'),
strain= list(dna_reads.keys()))
    output:
        gpa_csv=os.path.join('{roary_dir}', 'gene_presence_absence.csv'),
        gpa_rtab=os.path.join('{roary_dir}', 'gene_presence_absence.Rtab'),
        prot_tab=os.path.join('{roary_dir}', 'clustered_proteins')
    conda: '../shared_envs_yaml/perl5_22_env.yml'
    params:
        check_add_perl_env_script= 'install_perl_mods.sh',
        check_add_software_script= 'set_roary_env.sh',
        roary_bin= 'roary'
    threads: 16
    shell:
        '''
        set +u
#        {params.check_add_software_script}
        ROARY_HOME=$(dirname $(dirname $(which roary)))
        # required perl modules
        {params.check_add_perl_env_script}

	export PATH=\
$ROARY_HOME/build/fasttree:\
$ROARY_HOME/build/mcl-14-137/src/alien/oxygen/src:\
$ROARY_HOME/build/mcl-14-137/src/shmcl:\
$ROARY_HOME/build/ncbi-blast-2.4.0+/bin:\
$ROARY_HOME/build/prank-msa-master/src:\
$ROARY_HOME/build/cd-hit-v4.6.6-2016-0711:\
$ROARY_HOME/build/bedtools2/bin:\
$ROARY_HOME/build/parallel-20160722/src:\
$PATH
	export PERL5LIB=$ROARY_HOME/lib:\
$ROARY_HOME/build/bedtools2/lib:$PERL5LIB
#        PERL5LIB=$ROARY_HOME/lib:\
#$ROARY_HOME/build/bedtools2/lib:\
#$PERL5LIB
#	ln -rsf $ROARY_HOME/lib/Bio $CONDA_PREFIX/lib/perl5/5.22.0/
	which perl
        echo $PERL5LIB
        echo $PERLLIB
        rm -r {wildcards.roary_dir}
        {params.roary_bin} -f {wildcards.roary_dir} \
-v {input.gff_files} -p 30 -g 100000 -z
        set -u
        ''' 
#        {params.roary_bin} -f {wildcards.roary_dir} \
#-e -n -v {input.gff_files} -r -p 30 -g 100000 -z

rule create_gff:
    input: os.path.join(out_spades_dir,'{strain}', 'contigs.fasta')
    output: 
        os.path.join(out_prokka_dir, '{strain}', '{strain}.gff'), 
        os.path.join(out_prokka_dir, '{strain}', '{strain}.ffn')
    threads: 1
    conda: '../shared_envs_yaml/prokka_env.yml'
    shell:
        '''
#        echo $PERL5LIB
        which prokka
        prokka --locustag {wildcards.strain} \
--prefix  {wildcards.strain} \
--force  --cpus {threads} --metagenome --compliant \
--outdir prokka/{wildcards.strain} {input}
        '''

rule spades_create_assembly:
    input: 
        READS= lambda wildcards: [
            os.path.join(new_reads_dir,'{}.cleaned.{}.fq.gz'.format(
            wildcards.strain, str(n))) for n in [1,2]]
    output: os.path.join(out_spades_dir,'{strain}', 'contigs.fasta')
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: (2**attempt) * 10000
    params:
        spades_outdir= os.path.join(out_spades_dir, '{strain}'),
        SPADES_OPT='--careful',
        SPADES_BIN='spades.py'
    conda: '../shared_envs_yaml/spades_3_10_env.yml'
    script:'run_spades.py'

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
            ln  -s {input.infile1} {output.f1}
            ln  -s {input.infile2} {output.f2}
        fi
        '''

