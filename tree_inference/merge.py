'''
Load the five workflows one by one

NOTE:
    integrate the config arguments (many -> one config file)
'''
import snakemake
import os
class ngs_workflow():
    def __init__(self, config_params, target_f, snakefile, conda_prefix, workdir):
        self.config_params= config_params
        self.target_f= target_f
        self.snakefile= snakefile
        self.conda_prefix= conda_prefix
        self.workdir= workdir
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        import re 
        snakemake.snakemake(
            dryrun= just_dryrun,
            snakefile=self.snakefile,
            config= self.config_params,
            conda_prefix=self.conda_prefix,
            use_conda=True,
            workdir=self.workdir,
            cores=cpu_num,
            targets= ([] if re.search('\w', self.target_f) is None else
                      [self.target_f]))

class denovo(ngs_workflow):
    def __init__(self, list_f= '', project_dir= ''):
        print('Function: computing denovo assemblies and clustering orthologues\n...initiating')
        config={
            'list_f':list_f,
            'adaptor': '-',
            'new_reads_dir': (
                '/net/metagenomics/data/from_moni/old.tzuhao/'
                'TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/'
                'seq2geno/reads/dna'),
            'out_prokka_dir': os.path.join(project_dir,'prokka'),
            'out_roary_dir': os.path.join(project_dir,'roary'),
            'out_spades_dir': os.path.join(project_dir,'spades')
        }
        target_f=os.path.join(project_dir, 'roary', 'gene_presence_absence.csv')
        super().__init__(config, target_f,
            './denovo_workflow/denovo.in_one.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'denovo_workflow')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class gpa_aln(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: making gpa alignment\n...initiating')
        super().__init__(config_f, target_f,
            './gpa_workflow/gpa.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'gpa_workflow')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class cd_tr(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: making gpa alignment\n...initiating')
        super().__init__(config_f, target_f,
            './conc_workflow/conc_tree.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'conc_workflow/')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class nuc_tr(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: computing nucleotide tree\n...initiating')
        super().__init__(config_f, target_f,
            './nuc_workflow/nuc_tr.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'nuc_workflow/')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class mapping(ngs_workflow):
    def __init__(self, config_f= '', target_f= ''):
        print('Function: mapping\n...initiating')
        super().__init__(config_f, target_f,
            './snps_workflow/snps.in_one.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'snps_workflow/')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

#list_f=(
#    '/net/metagenomics/data/from_moni/'
#    'old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/'
#    'run_seq2geno/dna_list')
#project_dir= (
#    '/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
#    'TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/'
#    'seq2geno/denovo')

#denovo_obj= denovo(list_f, project_dir)
#denovo_obj.run_workflow(just_dryrun= True)

if __name__=='__main__':
    #def determine_workflow(x, config_f, target_f):
    def determine_workflow(x, list_f, project_dir):
        import sys
        target_func= ''
        if x == 'denovo':
            #target_func= denovo
            #return(denovo(config_f, target_f))
            return(denovo(list_f, project_dir))
#        elif x == 'gpa':
#            #target_func= gpa_aln
#            return(gpa_aln(config_f, target_f))
#        elif x == 'cd_tr':
#            #target_func= cd_tr 
#            return(cd_tr(config_f, target_f))
#        elif x == 'nuc_tr':
#            #target_func= nuc_tr
#            return(nuc_tr(config_f, target_f))
#        elif x == 'mapping':
#            #target_func= mapping
#            return(mapping(config_f, target_f))
        else:
            sys.exit('Unknown function')
        return(target_func)

    ## let the user opt the function
    import argparse
    parser = argparse.ArgumentParser(description=('''
                         Tree inference using composite datatype of microbial representations
                                                 '''))
#    parser.add_argument('--c',dest= 'config_f', type=str, required= True,
#            help='the config file for the workflow')
#    parser.add_argument('--o',dest= 'target_f', type=str, required= True,
#            help='the target file from the workflow')
    parser.add_argument('--p',dest= 'project_dir', type=str, required= True,
            help='the directory for project')
    parser.add_argument('--l',dest= 'list_f', type=str, required= True,
            help='the list of sequencing reads')
    parser.add_argument('--cpu',dest= 'cpu', type=str, default= 1,
            help='cpu number')
    parser.add_argument('--dry',dest= 'dryrun', action= 'store_true',
            help='cpu number')
    parser.add_argument('f', type=str,
            help='the subworkflow to launch',
            choices=['mapping', 'nuc_tr', 'cd_tr', 'gpa', 'denovo'])
    args= parser.parse_args()
    #target_func= determine_workflow(args.f, args.config_f, args.target_f)
    target_func= determine_workflow(args.f, args.list_f, args.project_dir)
    target_func.run_workflow(cpu_num= args.cpu, just_dryrun= args.dryrun)

