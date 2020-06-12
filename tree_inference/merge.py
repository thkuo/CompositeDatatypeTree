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

class cd_tr(ngs_workflow):
    def __init__(self, proj):
        print('Function: making gpa alignment\n...initiating')
        config={
            'nuc_aln': proj.roary_out,
            'gpa_aln': proj.gpa_aln,
            'conc_aln':os.path.join(proj.project_dir, 'cd', 'materials',
                                    'conc.aln'),
            'raxml_prtn':os.path.join(proj.project_dir, 'cd', 'materials',
                                    'conc.prtn'),
            'col_constraint_tr':proj.col_nuc_tr_out,
            'model': proj.cd_model
        }
        target_f=proj.cd_tr_out
        super().__init__(config, target_f,
            './conc_workflow/conc_tree.smk',
            proj.project_dir
                        )
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)
    def collapse_tree(self):
        tr_f= self.cd_tr_out
        out_tr_f= self.col_cd_tr_out
        import subprocess
        home_dir= os.path.dirname(os.path.realpath(__file__))
        script_f= os.path.join(home_dir, 'collapse_br', 'collapse_tree.R' )
        subprocess.run([script_f, '--i', tr_f, '--o', out_tr_f,
                        '--br', str(self.br_cutoff), '--bs', str(self.bs_cutoff), 
                        '--og', ' '.join(self.outgroup)])

class gpa_aln(ngs_workflow):
    def __init__(self, proj):
        print('Function: making gpa alignment\n...initiating')
        config={
            'roary_gpa': proj.roary_out,
            'strains': []
        }
        target_f=proj.gpa_aln
        super().__init__(config, target_f,
            './gpa_workflow/gpa.smk',
            proj.project_dir
                        )
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class denovo(ngs_workflow):
    def __init__(self, proj):
        print('Function: computing denovo assemblies and clustering orthologues\n...initiating')
        config={
            'list_f':proj.list_f,
            'adaptor': proj.adaptor,
            'new_reads_dir': proj.new_reads_dir,
            'out_prokka_dir': os.path.join(proj.project_dir,'prokka'),
            'out_roary_dir': os.path.join(proj.project_dir,'roary'),
            'out_spades_dir': os.path.join(proj.project_dir,'spades')
        }
        target_f=proj.roary_out
        super().__init__(config, target_f,
            './denovo_workflow/denovo.in_one.smk',
            proj.project_dir
                        )
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class nuc_tr(ngs_workflow):
    def __init__(self, proj):
        print('Function: computing nucleotide tree\n...initiating')
        config={
            'raxml_model':proj.nuc_model,
            'adaptor': proj.adaptor,
            'multisample_vcf': proj.vcf_out,
            'ref_fa': proj.ref
        }
        target_f=proj.nuc_tr_out
        super().__init__(config_f, target_f,
            './nuc_workflow/nuc_tr.smk',
            proj.project_dir
                        )
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)
    def collapse_tree(self):
        tr_f= self.cd_tr_out
        out_tr_f= self.col_cd_tr_out
        import subprocess
        home_dir= os.path.dirname(os.path.realpath(__file__))
        script_f= os.path.join(home_dir, 'collapse_br', 'collapse_tree.R' )
        subprocess.run([script_f, '--i', tr_f, '--o', out_tr_f,
                        '--br', str(self.br_cutoff), '--bs', str(self.bs_cutoff), 
                        '--og', ' '.join(self.outgroup)])

class mapping(ngs_workflow):
    def __init__(self, proj):
        print('Function: mapping\n...initiating')
        config={
            'list_f':proj.list_f,
            'adaptor': proj.adaptor,
            'new_reads_dir': proj.new_reads_dir,
            'ref_fasta': proj.ref
        }
        target_f=proj.vcf_out
        super().__init__(config, target_f,
            './snps_workflow/snps.in_one.smk',
            proj.project_dir
                        )
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)


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
    #def determine_workflow(x, config_f, target_f):
    def determine_workflow(x, list_f, project_dir, ref):
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
        elif x == 'mapping':
            #target_func= mapping
            #return(mapping(config_f, target_f))
            return(mapping(list_f, project_dir, ref))
        else:
            sys.exit('Unknown function')
        return(target_func)

    # read arguments
    from CDTreeArgParser import CDTreeArgParser
    parser= CDTreeArgParser()
    args= parser.parse()
    print(args.__dict__)
    ## control the output filenames
    from CDTreeProject import CDTreeProject
    cd_proj= CDTreeProject(
        args.list_f, args.project_dir, args.ref,
        args.adaptor, args.nuc_subs_model, args.rate_model,
        args.br_cutoff, args.bs_cutoff, args.outgroup)

    ## 
#    target_func= determine_workflow(args.f, cd_proj)
#    target_func.run_workflow(cpu_num= args.cpu, just_dryrun= args.dryrun)

