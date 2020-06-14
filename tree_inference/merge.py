'''
Interface to the five workflows 

Required variables:
    home: to determine the working directory 
    env: for the external packages
        default: ./shared_envs
        configurable using CDTREE_SHARED_ENV: 
            export \
            CDTREE_SHARED_ENV=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs

NOTE:
    The user determines: reads, reference, project directory, full list of strains(???), and outgroup(???) 
    This software determines output: 
        vcf (mapping)
        nuc tree
        gpa matrix (roary)
        gpa alignment
        cd tree
'''
import snakemake
import os
class ngs_workflow():
    def __init__(self, config_params, target_f, snakefile, workdir):
        self.config_params= config_params
        self.target_f= target_f
        self.snakefile= snakefile
        self.workdir= workdir
    def determine_smk(self):
        return(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            self.snakefile))

    def determine_conda_env_dir(self):
        return(os.environ['CDTREE_SHARED_ENV'] if 
                'CDTREE_SHARED_ENV' in os.environ else
                os.path.join(os.path.dirname(os.path.realpath(__file__)),'shared_envs'))
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        import re 
        print(self.determine_conda_env_dir())
        snakemake.snakemake(
            dryrun= just_dryrun,
            #snakefile=self.snakefile,
            snakefile=self.determine_smk(),
            config= self.config_params,
            conda_prefix=self.determine_conda_env_dir(),
            use_conda=True,
            workdir=self.workdir,
            cores=cpu_num,
            targets= ([] if re.search('\w', self.target_f) is None else
                      [self.target_f]))

class tr_inf_workflow(ngs_workflow):
    def collapse_tree(self):
        tr_f= str(self.cd_tr_out)
        out_tr_f= str(self.col_cd_tr_out)
        br_cutoff= str(self.br_cutoff)
        bs_cutoff= str(self.bs_cutoff)
        outgroup= ' '.join(self.outgroup)
        import subprocess
        home_dir= os.path.dirname(os.path.realpath(__file__))
        script_f= os.path.join(home_dir, 'collapse_br', 'collapse_tree.R' )
        subprocess.run([script_f, '--i', tr_f, '--o', out_tr_f,
                        '--br', br_cutoff, '--bs', bs_cutoff, 
                        '--og', outgroup])

#class cd_tr(ngs_workflow):
class cd_tr(tr_inf_workflow):
    def __init__(self, proj):
        print('Function: making gpa alignment\n...initiating')
        config={
            'nuc_aln': str(proj.nuc_aln),
            'gpa_aln': str(proj.gpa_aln),
            'conc_aln': str(proj.conc_aln),
            'raxml_prtn': str(proj.raxml_prtn),
            'col_constraint_tr': str(proj.col_nuc_tr_out),
            'model': str(proj.cd_model)
        }
        target_f= str(proj.cd_tr_out)
        workdir= str(proj.project_dir)
        super().__init__(config, target_f,
            './conc_workflow/conc_tree.smk',
            workdir
                        )
        print(self.__dict__)
#    def run_workflow(self, cpu_num= 1, just_dryrun= True):
#        super().run_workflow(cpu_num, just_dryrun)

class gpa_aln(ngs_workflow):
    def __init__(self, proj):
        print('Function: making gpa alignment\n...initiating')
        config={
            'roary_gpa': str(proj.roary_out)
        }
        target_f= str(proj.gpa_aln)
        workdir= str(proj.project_dir)
        super().__init__(config, target_f,
            './gpa_workflow/gpa.smk',
            workdir
                        )
        print(self.__dict__)
#    def run_workflow(self, cpu_num= 1, just_dryrun= True):
#        super().run_workflow(cpu_num, just_dryrun)

class denovo(ngs_workflow):
    def __init__(self, proj):
        print('Function: computing denovo assemblies and clustering orthologues\n...initiating')
        config={
            'list_f':str(proj.list_f),
            'adaptor': str(proj.adaptor),
            'new_reads_dir': str(proj.new_reads_dir),
            'out_prokka_dir': str(proj.prokka_dir),
            'out_roary_dir': str(proj.roary_dir),
            'out_spades_dir':str(proj.spades_dir)
        }
        target_f= str(proj.roary_out)
        workdir= str(proj.project_dir)
        super().__init__(config, target_f,
            './denovo_workflow/denovo.in_one.smk',
            workdir
                        )
        print(self.__dict__)
#    def run_workflow(self, cpu_num= 1, just_dryrun= True):
#        super().run_workflow(cpu_num, just_dryrun)

#class nuc_tr(ngs_workflow):
class nuc_tr(tr_inf_workflow):
    def __init__(self, proj):
        print('Function: computing nucleotide tree\n...initiating')
        config={
            'raxml_model':str(proj.nuc_model),
            'multisample_vcf': str(proj.multisample_vcf),
            'ref_fa': str(proj.ref)
        }
        target_f=str(proj.nuc_tr_out)
        workdir= str(proj.project_dir)
        super().__init__(config, target_f,
            './nuc_workflow/nuc_tr.smk',
            workdir
                        )
        print(self.__dict__)
#    def run_workflow(self, cpu_num= 1, just_dryrun= True):
#        super().run_workflow(cpu_num, just_dryrun)
#    def collapse_tree(self):
#        tr_f= self.cd_tr_out
#        out_tr_f= self.col_cd_tr_out
#        import subprocess
#        home_dir= os.path.dirname(os.path.realpath(__file__))
#        script_f= os.path.join(home_dir, 'collapse_br', 'collapse_tree.R' )
#        subprocess.run([script_f, '--i', tr_f, '--o', out_tr_f,
#                        '--br', str(self.br_cutoff), '--bs', str(self.bs_cutoff), 
#                        '--og', ' '.join(self.outgroup)])

class mapping(ngs_workflow):
    def __init__(self, proj):
        print('Function: mapping\n...initiating')
        config={
            'list_f':str(proj.list_f),
            'adaptor': str(proj.adaptor),
            'new_reads_dir': str(proj.new_reads_dir),
            'ref_fasta': str(proj.ref)
        }
        target_f=str(proj.multisample_vcf)
        workdir= str(proj.project_dir)
        super().__init__(config, target_f,
            './snps_workflow/snps.in_one.smk',
            workdir
                        )
        print(self.__dict__)
#    def run_workflow(self, cpu_num= 1, just_dryrun= True):
#        super().run_workflow(cpu_num, just_dryrun)


#list_f=(
#    '/net/metagenomics/data/from_moni/'
#    'old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/'
#    'run_seq2geno/dna_list')
#project_dir= (
#    '/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
#    'TreePaper/WhichTree_Sim.v7/results.v5/')
#ref= ('/net/metagenomics/data/from_moni/old.tzuhao/'
#      'TreePaper/WhichTree_Sim/data/reference/ATCC_700669.fasta')

#from CDTreeProject import CDTreeProject
#cd_proj= CDTreeProject(list_f, project_dir, ref)
#
### 
#print(cd_proj.__dict__)
#target_func= mapping(cd_proj)
#target_func.run_workflow(cpu_num= 1, just_dryrun= True)

#
if __name__=='__main__':
    def determine_workflow(x, p):
        import sys
        target_func= ''
        if x == 'denovo':
            return(denovo(p))
        elif x == 'gpa':
            return(gpa_aln(p))
        elif x == 'cd_tr':
            return(cd_tr(p))
        elif x == 'nuc_tr':
            return(nuc_tr(p))
        elif x == 'mapping':
            return(mapping(p))
        else:
            sys.exit('Unknown function')
        return(target_func)

    # read arguments
    from CDTreeArgParser import CDTreeArgParser
    parser= CDTreeArgParser()
    args= parser.parse()
    #print(args.__dict__)
    ## control the output filenames
    from CDTreeProject import CDTreeProject
    cd_proj= CDTreeProject(
        args.list_f, args.project_dir, args.ref)
    cd_proj.user_define(args.config_f)

    from pprint import pprint
    pprint(cd_proj.__dict__)
    ## 
    target_func= determine_workflow(args.f, cd_proj)
    target_func.run_workflow(cpu_num= args.cpu, just_dryrun= args.dryrun)

