## This class determines the directory of a project
'''
NOTE:
    raxml substitution model
'''

class CDTreeProject:
    def __init__(self, list_f= '', project_dir= '', ref= '', adaptor= '',
                 nuc_subs_model= 'GTR', rate_model= 'GAMMA',
                 br_cutoff= 5e-4, bs_cutoff=60, outgroup= []):
        import os
        self.list_f= list_f
        self.project_dir= project_dir
        self.ref= ref
        if not os.path.isfile(self.list_f):
            raise FileNotFoundError('input reads list not found')
        if not os.path.isfile(self.ref):
            raise FileNotFoundError('input reference not found')
        if os.makedirs(self.project_dir) and os.listdir(self.project_dir):
            print('WARNING: {} not empty'.format(self.project_dir) )
        elif not os.path.isdir(self.project_dir):
            os.makedirs(self.project_dir)


        self.adaptor= '-' if adaptor == '' else adaptor
        ## shared data
        self.new_reads_dir= os.path.join(project_dir, 'dna_reads_processed')
        ## mapping
        self.vcf_out=os.path.join(project_dir,
                                  'mapping', 'merged_vcf',
                                  'multisample.snp.vcf.gz')
        ## denovo
        self.roary_out=os.path.join(project_dir, 'roary',
                                    'gene_presence_absence.csv')
        ## gpa aln
        self.roary_out=os.path.join(project_dir, 'gpa', 'gpa.var.aln')
        ## nuc tree
        self.nuc_tr_out=os.path.join(project_dir,
                                     'mapping','raxml',
                                     'RAxML_bipartitions.nuc.bs')
        self.nuc_subs_model= nuc_subs_model 
        self.rate_model= rate_model
        self.nuc_model=''.join(self.nuc_subs_model, self.rate_model)
        self.col_nuc_tr_out='{}_col'.format(self.nuc_tr_out)
        ## cd tree
        self.cd_tr_out= os.path.join(project_dir,
                                     'cd','RAxML_bipartitions.conc.bs')
        # only the rate model will be extracted by RAxML
        self.cd_model='BIN{}'.format(self.rate_model)
        self.col_cd_tr_out= '{}_col'.format(self.cd_tr_out)

        ## tree post processing
        self.outgroup= outgroup
        self.br_cutoff= br_cutoff
        self.bs_cutoff= bs_cutoff

