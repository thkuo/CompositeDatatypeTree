import math
tr_f = config['tr_f']
total_size = config['total_size']
var_cutoff = config['var_cutoff']
partitions_num = int(config['partitions_num'])
# divide the gpa alignments
sizes = [total_size / partitions_num] * partitions_num
for n in range(partitions_num - 1):
    # except for the last
    sizes[n] = math.floor(sizes[n])
# the last implements the others
sizes[partitions_num - 1] = total_size - sum(sizes[:(partitions_num - 2)])


rule all:
    input:
        gpa_aln='{wd}/gpa.var.aln'


rule sim_gpa:
    input:
        tr_f=tr_f
    output:
        temp('{wd}/gpa.var.aln-{p}')
    params:
        p=lambda wildcards: int(wildcards.p),
        gain_r=1.0,
        loss_r=0.12,
        sim_num=lambda wildcards: int(sizes[int(wildcards.p)]),
        a0=0.8,
        cutoff=var_cutoff
    script: 'sim_gpa_snakemake.R'


rule integrate_sub_gpa_aln:
    input:
        sub_gpa_aln_files=expand('{{wd}}/gpa.var.aln-{p}', p=range(partitions_num))
    output:
        gpa_aln='{wd}/gpa.var.aln'
    run:
        from Bio import AlignIO
        gpa = AlignIO.read(input.sub_gpa_aln_files[0], 'fasta')
        gpa.sort()
        for n in range(1, len(input.sub_gpa_aln_files)):
            sub_gpa = AlignIO.read(input.sub_gpa_aln_files[n], 'fasta')
            sub_gpa.sort()
            gpa = gpa + sub_gpa
        AlignIO.write(gpa, output.gpa_aln, 'fasta')
