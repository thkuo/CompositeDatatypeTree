#!/usr/bin/env python3
'''
Interface to each workflow
'''
import snakemake
import os
import logging
import pprint
import re
from CDTreeArgParser import CDTreeArgParser
from CDTreeProject import CDTreeProject

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s %(levelname)s: %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p')
ch.setFormatter(formatter)
logger.addHandler(ch)

class ngs_workflow():
    def __init__(self, config_params, target_f, snakefile, workdir):
        self.config_params = config_params
        self.target_f = target_f
        self.snakefile = snakefile
        self.workdir = workdir

    def determine_smk(self):
        return(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            self.snakefile))

    def determine_conda_env_dir(self):
        env_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               'shared_envs')
        if 'CDTREE_SHARED_ENV' in os.environ:
            env_dir = os.environ['CDTREE_SHARED_ENV']
        return(env_dir)

    def run_workflow(self, cpu_num=1, just_dryrun=True):
        # ensure unlocked folder
        snakemake.snakemake(
            unlock=True,
            snakefile=self.determine_smk(),
            config=self.config_params,
            conda_prefix=self.determine_conda_env_dir(),
            use_conda=True,
            workdir=self.workdir,
            cores=cpu_num,
            targets=([] if re.search('\w', self.target_f) is None else
                     [self.target_f]))
        # run the workflow
        snakemake.snakemake(
            dryrun=just_dryrun,
            printshellcmds=True,
            restart_times=3,
            force_incomplete=True,
            snakefile=self.determine_smk(),
            config=self.config_params,
            conda_prefix=self.determine_conda_env_dir(),
            use_conda=True,
            workdir=self.workdir,
            cores=cpu_num,
            targets=([] if re.search('\w', self.target_f) is None else
                     [self.target_f]))


class tr_inf_workflow(ngs_workflow):
    def run_workflow(self, cpu_num=1, just_dryrun=True):
        super().run_workflow(cpu_num, just_dryrun)


class cd_tr(tr_inf_workflow):
    def __init__(self, proj):
        logging.info('Set up composite datatype inference')
        config = {
            'nuc_aln': str(proj.nuc_aln),
            'gpa_aln': str(proj.gpa_aln),
            'conc_aln': str(proj.conc_aln),
            'raxml_prtn': str(proj.raxml_prtn),
            'col_constraint_tr': str(proj.col_constraint_tr),
            'raxml_model': str(proj.cd_model)
        }
        target_f = str(proj.cd_tr_out)
        workdir = str(proj.project_dir)
        self.outgroup = proj.outgroup
        super().__init__(config, target_f,
                         './conc_workflow/conc_tree.smk',
                         workdir)
        logger.info(pprint.pformat(self.__dict__))

class gpa_aln(ngs_workflow):
    def __init__(self, proj):
        logging.info('Set up the generation of gpa alignment')
        config = {
            'roary_gpa': str(proj.roary_out)
        }
        target_f = str(proj.gpa_aln)
        workdir = str(proj.project_dir)
        super().__init__(config, target_f,
                         './gpa_workflow/gpa.smk',
                         workdir)
        logger.info(pprint.pformat(self.__dict__))

class denovo(ngs_workflow):
    def __init__(self, proj):
        logging.info('''Set up denovo assemblies and orthologous
                     clustering''')
        config = {
            'list_f': str(proj.list_f),
            'adaptor': str(proj.adaptor),
            'new_reads_dir': str(proj.new_reads_dir),
            'out_prokka_dir': str(proj.prokka_dir),
            'out_roary_dir': str(proj.roary_out),
            'out_spades_dir': str(proj.spades_dir)}
        target_f = str(proj.roary_out)
        workdir = str(proj.project_dir)
        super().__init__(config, target_f,
                         './denovo_workflow/denovo.in_one.smk',
                         workdir)
        logger.info(pprint.pformat(self.__dict__))


class col_tr(tr_inf_workflow):
    def __init__(self, proj):
        logging.info('Set up branch collapsing for the nucleotide tree '
                     'using the log-likelihood method')
        config = {
            'nuc_tr': str(proj.nuc_tr_out),
            'nuc_aln': str(proj.nuc_aln),
            'raxml_model': str(proj.nuc_model),
            'col_nuc_tr': str(proj.col_nuc_tr_out),
            'cutoff_perc': proj.cutoff_perc
        }
        target_f = str(proj.col_nuc_tr_out)
        workdir = str(proj.project_dir)
        super().__init__(config, target_f,
                         './LogLikelihood/compute_ll.smk',
                         workdir)
        logger.info(pprint.pformat(self.__dict__))


class nuc_tr(tr_inf_workflow):
    def __init__(self, proj):
        logging.info('Set up nucleotide tree inference')
        config = {
            'raxml_model': str(proj.nuc_model),
            'multisample_vcf': str(proj.multisample_vcf),
            'ref_fa': str(proj.ref)
        }
        target_f = str(proj.nuc_tr_out)
        workdir = str(proj.project_dir)
        self.outgroup = proj.outgroup
        super().__init__(config, target_f,
                         './nuc_workflow/nuc_tr.smk',
                         workdir)
        logger.info(pprint.pformat(self.__dict__))

class mapping(ngs_workflow):
    def __init__(self, proj, sub_func):
        logging.info('Set up read mapping ({})'.format(sub_func))
        config = {
            'list_f': str(proj.list_f),
            'adaptor': str(proj.adaptor),
            'new_reads_dir': str(proj.new_reads_dir),
            'ref_fasta': str(proj.ref)
        }
        target_f = str(proj.multisample_vcf)
        workdir = str(proj.project_dir)
        smk_file = ''
        if sub_func == 'mapping':
            smk_file = './snps_workflow/snps.in_one.smk'
        elif sub_func == 'fast_mapping':
            smk_file = './snps_workflow/snps_bwa_mem.in_one.smk'

        super().__init__(config, target_f,
                         smk_file, workdir)
        logger.info(pprint.pformat(self.__dict__))


def determine_workflow(opted_funcs, p):
    import sys
    target_funcs = []
    # include all as opted
    if 'all' in opted_funcs:
        opted_funcs = ['mapping', 'nuc_tr',
                       'col_tr', 'denovo', 'gpa', 'cd_tr']
    funcs_order = ['mapping', 'fast_mapping', 'nuc_tr',
                   'col_tr', 'denovo', 'gpa', 'cd_tr']
    funcs_classes = {'mapping':mapping,
                     'fast_mapping':mapping,
                     'nuc_tr':nuc_tr,
                     'col_tr':col_tr,
                     'denovo':denovo,
                     'gpa':gpa_aln,
                     'cd_tr':cd_tr}
    for x in funcs_order:
        if x in opted_funcs:
            logger.info("{} is opted".format(x))
            try:
                target_funcs.append(funcs_classes[x](p))
            except TypeError:
                target_funcs.append(funcs_classes[x](p, x))
            except Exception as e:
                logger.error("Unexpected error: {}".format(e))
                sys.exit()
    return(target_funcs)


if __name__ == '__main__':
    # read arguments
    parser = CDTreeArgParser()
    args = parser.parse()
    cd_proj = CDTreeProject(
        args.list_f, args.project_dir, args.ref)
    if (args.config_f is not None):
        cd_proj.user_define(args.config_f)

    logging.info('Parse the parameters \n{}'.format(
        pprint.pformat(cd_proj.__dict__)))

    target_funcs = determine_workflow(args.f, cd_proj)
    for target_func in target_funcs:
        target_func.run_workflow(cpu_num=args.cpu,
                                 just_dryrun=args.dryrun)
