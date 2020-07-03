'''
Use all strains in the file
'''
import os
roary_gpa_csv= config['roary_gpa']
#strains= config['strains'].strip().split(',')

rule trim_invariant:
    input:
        gpa_aln='{results_dir}/gpa/gpa.aln'
    output:
        gpa_var_aln='{results_dir}/gpa/gpa.var.aln'
    threads: 1
    params:
        cutoff_num= 2
    run:
        
        from Bio import SeqIO
        from Bio import AlignIO
        aln_f= input['gpa_aln']
        out_f= output['gpa_var_aln']
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

rule create_gpa_aln:
    input:
        roary_raw_gpa=roary_gpa_csv
    output:
        gpa_aln='{results_dir}/gpa/gpa.aln'
    script:  'createGPAaln.py'
