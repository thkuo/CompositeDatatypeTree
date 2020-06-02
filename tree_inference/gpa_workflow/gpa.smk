roary_gpa_csv= config['roary_gpa']
strains= config['strains'].strip().split(',')

rule trim_invariant:
    input:
        gpa_aln='{results_dir}/gpa/gpa.aln'
    output:
        gpa_var_aln='{results_dir}/gpa/gpa.var.aln'
    threads: 1
    shell:
        '''
        ./removeInvariant.py --in {input} \
--out {output} --cn 2
        '''

rule create_gpa_aln:
    input:
        roary_raw_gpa=roary_gpa_csv
    output:
        gpa_aln='{results_dir}/gpa/gpa.aln'
    params: 
        strains= strains
    script: 'createGPAaln.py'
