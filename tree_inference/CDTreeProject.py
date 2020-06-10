## This class determines the directory of a project

class CDTreeProject:
    def __init__(self, list_f= '', project_dir= '', ref= '', adaptor= ''):
        import os
        self.list_f= list_f
        self.project_dir= project_dir
        self.ref= ref
        self.adaptor= '-' if adaptor == '' else adaptor
        ## shared data
        self.new_reads_dir= os.path.join(project_dir, 'dna_reads_processed')
        ## mapping
        self.vcf_out=os.path.join(project_dir,
                                  'mapping','merged_vcf/multisample.snp.vcf.gz')
        ## denovo
        self.roary_out=os.path.join(project_dir, 'roary',
                                    'gene_presence_absence.csv')
        ## gpa aln
        self.roary_out=os.path.join(project_dir, 'gpa', 'gpa.var.aln')
        ## nuc tree
        self.nuc_tr_out=os.path.join(project_dir,
                                     'mapping','raxml/RAxML_bipartitions.nuc.bs')
        ## cd tree
        self.cd_tr_out= os.path.join(project_dir,
                                     'cd','RAxML_bipartitions.conc.bs')



