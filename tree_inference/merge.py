#!/usr/bin/env python3
'''
Interface to each workflow
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
        env_dir= os.path.join(os.path.dirname(os.path.realpath(__file__)),'shared_envs')
        if 'CDTREE_SHARED_ENV' in os.environ:
            env_dir=os.environ['CDTREE_SHARED_ENV'] 
        return(env_dir)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        import re 
        snakemake.snakemake(
            dryrun= just_dryrun,
            printshellcmds= True, 
            snakefile=self.determine_smk(),
            config= self.config_params,
            conda_prefix=self.determine_conda_env_dir(),
            use_conda=True,
            workdir=self.workdir,
            cores=cpu_num,
            targets= ([] if re.search('\w', self.target_f) is None else
                      [self.target_f]))

class tr_inf_workflow(ngs_workflow):
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

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
#        self.br_cutoff= proj.br_cutoff
#        self.bs_cutoff= proj.bs_cutoff
        self.outgroup= proj.outgroup
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

class col_tr(tr_inf_workflow):
    def __init__(self, proj):
        print('Function: collapsing the nucleotide tree using log-likelihood '
              'score\n...initiating')
        config={
            'nuc_tr': str(proj.nuc_tr_out),
            'nuc_aln': str(proj.nuc_aln),
            'raxml_model':str(proj.nuc_model),
            'col_nuc_tr': str(proj.col_nuc_tr_out),
            'cutoff_perc': proj.cutoff_perc
        }
        target_f=str(proj.col_nuc_tr_out)
        workdir=str(proj.project_dir)
        super().__init__(config, target_f,
            './per-site_LogLikelihood/compute_psll.smk',
            workdir
                        )


class nuc_tr(tr_inf_workflow):
    def __init__(self, proj):
        print('Function: computing nucleotide tree\n...initiating')
        config={
            'raxml_model':str(proj.nuc_model),
            'multisample_vcf': str(proj.multisample_vcf),
            'ref_fa': str(proj.ref)
        }
        target_f=str(proj.nuc_tr_out)
        workdir=str(proj.project_dir)
#        self.br_cutoff= proj.br_cutoff
#        self.bs_cutoff= proj.bs_cutoff
        self.outgroup= proj.outgroup
        super().__init__(config, target_f,
            './nuc_workflow/nuc_tr.smk',
            workdir
                        )

class mapping(ngs_workflow):
    def __init__(self, proj, sub_func):
        print('Function: mapping\n...initiating')
        config={
            'list_f':str(proj.list_f),
            'adaptor': str(proj.adaptor),
            'new_reads_dir': str(proj.new_reads_dir),
            'ref_fasta': str(proj.ref)
        }
        target_f=str(proj.multisample_vcf)
        workdir= str(proj.project_dir)
        smk_file= '' 
        if sub_func == 'mapping':
            smk_file= './snps_workflow/snps.in_one.smk' 
        elif sub_func == 'fast_mapping':
            smk_file= './snps_workflow/snps_bwa_mem.in_one.smk'
        super().__init__(config, target_f,
            smk_file, workdir )
        print(self.__dict__)
#
if __name__=='__main__':
    def determine_workflow(x, p):
        import sys
        target_funcs= []
        if x == 'denovo':
            target_funcs.append(denovo(p))
        elif x == 'gpa':
            target_funcs.append(gpa_aln(p))
        elif x == 'cd_tr':
            target_funcs.append(cd_tr(p))
        elif x == 'nuc_tr':
            target_funcs.append(nuc_tr(p))
        elif x == 'col_tr':
            target_funcs.append(col_tr(p))
        elif x in ['mapping', 'fast_mapping']:
            target_funcs.append(mapping(p, x))
        elif x == 'all':
            target_funcs= [mapping(p), nuc_tr(p), denovo(p), gpa_aln(p),
                           cd_tr(p)]
        else:
            sys.exit('Unknown function')
        return(target_funcs)

    # read arguments
    from CDTreeArgParser import CDTreeArgParser
    parser= CDTreeArgParser()
    args= parser.parse()
    #print(args.__dict__)
    ## determine parameters
    from CDTreeProject import CDTreeProject
    cd_proj= CDTreeProject(
        args.list_f, args.project_dir, args.ref)
    if (not args.config_f is None):
        cd_proj.user_define(args.config_f)

    from pprint import pprint
    print('==parameters==')
    pprint(cd_proj.__dict__)

    ## 
    target_funcs= determine_workflow(args.f, cd_proj)
    for target_func in target_funcs:
        target_func.run_workflow(cpu_num= args.cpu, just_dryrun= args.dryrun)


