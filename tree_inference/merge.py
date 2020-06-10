'''
Interface to the five workflows 

NOTE:
    The user determines: reads, reference, outgroup(???), and project directory
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
            'col_constraint_tr':os.path.join(proj.project_dir, 'cd', 
            'best_ml_tree':os.path.join(proj.project_dir, 'cd', 
            'bs_trees':os.path.join(proj.project_dir, 'cd', 
            'bs_best_ml_tree':os.path.join(proj.project_dir, 'cd', 
            'model': 'BINGAMMA'
        }
        target_f=proj.cd_tr
        super().__init__(config, target_f,
            './conc_workflow/conc_tree.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'conc_workflow/')
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
            'strains': 
        }
        target_f=proj.gpa_aln
        super().__init__(config, target_f,
            './gpa_workflow/gpa.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'gpa_workflow')
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
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'denovo_workflow')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)

class nuc_tr(ngs_workflow):
    def __init__(self, proj):
        print('Function: computing nucleotide tree\n...initiating')
        config={
            'raxml_model':'GTRGAMMA ',
            'adaptor': proj.adaptor,
            'multisample_vcf': proj.vcf_out,
            'ref_fa': proj.ref
        }
        target_f=proj.nuc_tr_out
        super().__init__(config_f, target_f,
            './nuc_workflow/nuc_tr.smk',
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'nuc_workflow/')
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
            '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs',
            'snps_workflow/')
        print(self.__dict__)
    def run_workflow(self, cpu_num= 1, just_dryrun= True):
        super().run_workflow(cpu_num, just_dryrun)


list_f=(
    '/net/metagenomics/data/from_moni/'
    'old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/'
    'run_seq2geno/dna_list')
project_dir= (
    '/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
    'TreePaper/WhichTree_Sim.v7/results.v5/')
ref= ('/net/metagenomics/data/from_moni/old.tzuhao/'
      'TreePaper/WhichTree_Sim/data/reference/ATCC_700669.fasta')

from CDTreeProject import CDTreeProject
cd_proj= CDTreeProject(list_f, project_dir, ref)

## 
print(cd_proj.__dict__)
target_func= mapping(cd_proj)
target_func.run_workflow(cpu_num= 1, just_dryrun= True)

#
#if __name__=='__main__':
#    #def determine_workflow(x, config_f, target_f):
#    def determine_workflow(x, list_f, project_dir, ref):
#        import sys
#        target_func= ''
#        if x == 'denovo':
#            #target_func= denovo
#            #return(denovo(config_f, target_f))
#            return(denovo(list_f, project_dir))
##        elif x == 'gpa':
##            #target_func= gpa_aln
##            return(gpa_aln(config_f, target_f))
##        elif x == 'cd_tr':
##            #target_func= cd_tr 
##            return(cd_tr(config_f, target_f))
##        elif x == 'nuc_tr':
##            #target_func= nuc_tr
##            return(nuc_tr(config_f, target_f))
#        elif x == 'mapping':
#            #target_func= mapping
#            #return(mapping(config_f, target_f))
#            return(mapping(list_f, project_dir, ref))
#        else:
#            sys.exit('Unknown function')
#        return(target_func)
#
#    ## let the user opt the function
#    import argparse
#    parser = argparse.ArgumentParser(description=('''
#                         Tree inference using composite datatype of microbial representations
#                                                 '''))
#    parser.add_argument('--p',dest= 'project_dir', type=str, required= True,
#            help='the directory for project')
#    parser.add_argument('--l',dest= 'list_f', type=str, required= True,
#            help='the list of sequencing reads')
#    parser.add_argument('--r',dest= 'ref', type=str, required= True,
#            help='the reference genome for mapping sequencing reads')
#    parser.add_argument('--cpu',dest= 'cpu', type=str, default= 1,
#            help='cpu number')
#    parser.add_argument('--dry',dest= 'dryrun', action= 'store_true',
#            help='cpu number')
#    parser.add_argument('f', type=str,
#            help='the subworkflow to launch',
#            choices=['mapping', 'nuc_tr', 'cd_tr', 'gpa', 'denovo'])
#    args= parser.parse_args()
#    ## control the output filenames
#    from CDTreeProject import CDTreeProject
#    cd_proj= CDTreeProject(args.list_f, args.project_dir, args.ref)
#
#    ## 
#    target_func= determine_workflow(args.f, cd_proj)
#    target_func.run_workflow(cpu_num= args.cpu, just_dryrun= args.dryrun)

