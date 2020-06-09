'''
Load the five workflows one by one

NOTE:
    integrate the config arguments (many -> one config file)
'''
import snakemake

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
        super().__init__(config_f, target_f,
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

if __name__=='__main__':
    def determine_workflow(x, config_f, target_f):
        import sys
        target_func= ''
        if x == 'denovo':
            #target_func= denovo
            return(denovo(config_f, target_f))
        elif x == 'gpa':
            #target_func= gpa_aln
            return(gpa_aln(config_f, target_f))
        elif x == 'cd_tr':
            #target_func= cd_tr 
            return(cd_tr(config_f, target_f))
        elif x == 'nuc_tr':
            #target_func= nuc_tr
            return(nuc_tr(config_f, target_f))
        elif x == 'mapping':
            #target_func= mapping
            return(mapping(config_f, target_f))
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
    parser.add_argument('--dry',dest= 'dryrun', action= 'store_true',
            help='cpu number')
    parser.add_argument('f', type=str,
            help='the subworkflow to launch',
            choices=['mapping', 'nuc_tr', 'cd_tr', 'gpa', 'denovo'])
    args= parser.parse_args()
    target_func= determine_workflow(args.f, args.config_f, args.target_f)
    target_func.run_workflow(cpu_num= args.cpu, just_dryrun= args.dryrun)

