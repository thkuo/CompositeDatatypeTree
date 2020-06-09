'''
Load the five workflows one by one

NOTE:
    integrate the config arguments (many -> one config file)
'''
import snakemake

'''
# reads -> snps
def mapping():
    print('Function: mapping...')
    just_dryrun= True
    config_f= ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
               'WhichTree_Sim.v7/bin.v5/snps_workflow/config.yml')
    target_f=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
              'WhichTree_Sim.v7/results.v5/mapping/merged_vcf/multisample.snp.vcf.gz')
    snakemake.snakemake(
        dryrun= just_dryrun,
        snakefile='snps_workflow/snps.in_one.smk',
        configfile= config_f,
        conda_prefix=('/net/metagenomics/data/from_moni/old.tzuhao/'
                      'TreePaper/shared_envs'),
        use_conda=True,
        workdir='snps_workflow/',
        cores=16,
        targets= [target_f])

def nuc_tr():
    print('Function: computing nucleotide tree...')
    just_dryrun= True
    config_f= ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
               'WhichTree_Sim.v7/bin.v5/nuc_workflow/nuc_config.yml')
    target_f=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
              'WhichTree_Sim.v7/results.v5/mapping/raxml/RAxML_bipartitions.nuc.bs')
    snakemake.snakemake(
        dryrun= just_dryrun,
        snakefile='nuc_workflow/nuc_tr.smk',
        configfile= config_f,
        conda_prefix=('/net/metagenomics/data/from_moni/old.tzuhao/'
                      'TreePaper/shared_envs'),
        use_conda=True,
        workdir='nuc_workflow/',
        cores=16,
        targets= [target_f])

def cd_tr():
    print('Function: computing composite datatype tree...')
    just_dryrun= True
    config_f=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
              'TreePaper/WhichTree_Sim.v7/bin.v5/conc_workflow/conc.config.yml')
    target_f= ''
    snakemake.snakemake(
        dryrun= just_dryrun,
        snakefile='conc_workflow/conc_tree.smk',
        configfile= config_f,
        conda_prefix=('/net/metagenomics/data/from_moni/old.tzuhao/'
                      'TreePaper/shared_envs'),
        use_conda=True,
        workdir='conc_workflow/',
        cores=16,
        targets= [])

def gpa_aln():
    print('Function: making gpa alignment...')
    just_dryrun= True
    config_f=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/gpa_workflow/config.yml')
    target_f=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/gpa/gpa.var.aln')
    snakemake.snakemake(
        dryrun= just_dryrun,
        snakefile='gpa_workflow/gpa.smk',
        configfile= config_f,
        conda_prefix=('/net/metagenomics/data/from_moni/old.tzuhao/'
                      'TreePaper/shared_envs'),
        use_conda=True,
        workdir='gpa_workflow',
        cores=16,
        targets= [target_f])

def denovo():
    print('Function: computing gpa...')
    just_dryrun= True
    config_f=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/seq2geno/denovo/denovo_config.yml')
    target_f=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/seq2geno/denovo/roary/gene_presence_absence.csv')
    snakemake.snakemake(
        dryrun= just_dryrun,
        snakefile='./denovo_workflow/denovo.in_one.smk',
        configfile= config_f,
        conda_prefix=('/net/metagenomics/data/from_moni/old.tzuhao/'
                      'TreePaper/shared_envs'),
        use_conda=True,
        workdir='denovo_workflow',
        cores=16,
        targets= [target_f])
'''

class ngs_workflow():
    def __init__(self, config_f, target_f, snakefile, conda_prefix, workdir):
        self.config_f= config_f
        self.target_f= target_f
        self.snakefile= snakefile
        self.conda_prefix= conda_prefix
        self.workdir= workdir
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        import re 
        snakemake.snakemake(
            dryrun= just_dryrun,
            snakefile=self.snakefile,
            configfile= self.config_f,
            conda_prefix=self.conda_prefix,
            use_conda=True,
            workdir=self.workdir,
            cores=cpu_num,
            targets= ([] if re.search('\w', self.target_f) is None else
                      [self.target_f]))

class denovo(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: computing gpa\n...initiating')
        super().__init__(
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
             'TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/'
             'seq2geno/denovo/denovo_config.yml'), 
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
             'TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/'
             'seq2geno/denovo/roary/gene_presence_absence.csv'),
            './denovo_workflow/denovo.in_one.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'denovo_workflow')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class gpa_aln(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: making gpa alignment\n...initiating')
        super().__init__(
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
             'TreePaper/WhichTree_Sim.v7/bin.v5/gpa_workflow/config.yml'), 
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
             'TreePaper/WhichTree_Sim.v7/results.v5/gpa/gpa.var.aln'),
            './gpa_workflow/gpa.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'gpa_workflow')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class cd_tr(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: making gpa alignment\n...initiating')
        super().__init__(
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
              'TreePaper/WhichTree_Sim.v7/bin.v5/conc_workflow/conc.config.yml'), 
            '',
            './conc_workflow/conc_tree.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'conc_workflow/')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class nuc_tr(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: computing nucleotide tree\n...initiating')
        super().__init__(
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
               'WhichTree_Sim.v7/bin.v5/nuc_workflow/nuc_config.yml'), 
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
              'WhichTree_Sim.v7/results.v5/mapping/raxml/RAxML_bipartitions.nuc.bs'),
            './nuc_workflow/nuc_tr.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'nuc_workflow/')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class mapping(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: mapping\n...initiating')
        super().__init__(
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
               'WhichTree_Sim.v7/bin.v5/snps_workflow/config.yml'), 
            ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
              'WhichTree_Sim.v7/results.v5/mapping/merged_vcf/multisample.snp.vcf.gz'),
            './snps_workflow/snps.in_one.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'snps_workflow/')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

if __name__=='__main__':
    def determine_workflow(x, config_f, target_f):
        import sys
        target_func= ''
        if x == 'denovo':
            #target_func= denovo
            return(denovo())
        elif x == 'gpa':
            #target_func= gpa_aln
            return(gpa_aln())
        elif x == 'cd_tr':
            #target_func= cd_tr 
            return(cd_tr())
        elif x == 'nuc_tr':
            #target_func= nuc_tr
            return(nuc_tr())
        elif x == 'mapping':
            #target_func= mapping
            return(mapping())
        else:
            sys.exit('Unknown function')
        return(target_func)
    ## let the user opt the function
    import argparse
    parser = argparse.ArgumentParser(description=('''
                         Tree inference using composite datatype of microbial representations
                                                 '''))
    parser.add_argument('--c',dest= 'config_f', type=str, required= True,
            help='the config file for the workflow')
    parser.add_argument('--o',dest= 'target_f', type=str, required= True,
            help='the target file from the workflow')
    parser.add_argument('--cpu',dest= 'cpu', type=str, default= 1,
            help='cpu number')
    parser.add_argument('f', type=str,
            help='the subworkflow to launch',
            choices=['mapping', 'nuc_tr', 'cd_tr', 'gpa', 'denovo'])
    args= parser.parse_args()
    target_func= determine_workflow(args.f, args.config_f, args.target_f)
    target_func.run_workflow(cpu_num= args.cpu)

