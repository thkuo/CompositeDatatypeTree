# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

# A class that contains all settings of a project

class CDTreeProject:
    def __init__(self, list_f= '', project_dir= '', ref= '', adaptor= '-'):
        import os
        from pathlib import Path

        self.list_f= Path(os.path.abspath(list_f))
        self.project_dir= Path(os.path.abspath(project_dir))
        self.ref= Path(os.path.abspath(ref))
        self.adaptor= Path(os.path.abspath(adaptor))
        if not self.list_f.is_file():
            raise FileNotFoundError('input reads list not found')
        if not self.ref.is_file():
            raise FileNotFoundError('input reference not found')
        if self.project_dir.is_dir() and len(os.listdir(self.project_dir))>0:
            print('WARNING: {} not empty'.format(self.project_dir) )
        elif not self.project_dir.is_dir():
            os.makedirs(self.project_dir)
        if str(os.path.basename(self.adaptor)) != '-' and not Path(self.adaptor).is_file():
            raise FileNotFoundError(
                'input adaptor file not found: {}'.format(self.adaptor))

        ##=== Default intermediate files/params
        ## shared data
        self.new_reads_dir= Path(os.path.join(self.project_dir,
                                              'dna_reads_processed'))
        ## mapping
        self.multisample_vcf= Path(os.path.join(self.project_dir,
                                  'mapping', 'merged_vcf',
                                  'multisample.snp.vcf.gz'))
        ## denovo
        self.prokka_dir= Path(os.path.join(self.project_dir,'prokka'))
        self.roary_dir= Path(os.path.join(self.project_dir,'roary'))
        self.spades_dir= Path(os.path.join(self.project_dir,'spades'))
        self.roary_out=Path(os.path.join(self.roary_dir,
                                    'gene_presence_absence.csv'))
        ## gpa aln
        self.gpa_aln=Path(os.path.join(self.project_dir, 'gpa', 'gpa.var.aln'))
        ## nuc tree
        self.nuc_aln= Path(os.path.join(self.project_dir,
                                     'mapping', 'alignment', 'nuc.var.aln'))
        self.nuc_subs_model= str('GTR')
        self.rate_model= str('GAMMA')
        self.nuc_model=str(''.join([self.nuc_subs_model, self.rate_model]))
        self.nuc_tr_out= Path(os.path.join(self.project_dir,
                                     'mapping','raxml',
                                     'RAxML_bipartitions.nuc.bs'))
        ## collapsed nuc tree
        self.col_nuc_tr_out=Path(os.path.join(self.project_dir,
                                              'LogLikelihood',
                                              'col_nuc.nwk'))
        self.cutoff_perc= 0.75

        ## cd tree
        self.conc_aln= Path(os.path.join(self.project_dir, 'cd', 'materials',
                      'conc.aln'))
        self.raxml_prtn= Path(os.path.join(self.project_dir, 'cd', 'materials',
                               'conc.prtn'))
        self.cd_tr_out= Path(os.path.join(self.project_dir,
                                     'cd','RAxML_bipartitions.conc.bs'))
        self.col_constraint_tr= Path(self.col_nuc_tr_out)

        # only the rate model will be extracted by RAxML
        self.cd_model=str('ASC_BIN{} --asc-corr lewis'.format(self.rate_model))
        self.col_cd_tr_out= Path('{}_col'.format(self.cd_tr_out))
        ## tree post processing
        self.outgroup= [] 
    def user_define(self, config_f):
        import yaml
        import os
        from pathlib import PosixPath
        from pathlib import Path
        config_dict= yaml.safe_load(open(config_f, 'r'))
        args_dict= self.__dict__
        for k in [k for k in args_dict if k in config_dict]:
            if isinstance(args_dict[k], PosixPath):
                print('{} --> {}'.format(args_dict[k], config_dict[k]))
                args_dict[k]= Path(config_dict[k])
            elif isinstance(args_dict[k], str):
                args_dict[k]= str(config_dict[k])
            elif isinstance(args_dict[k], list):
                args_dict[k]= config_dict[k].split(' ')
            elif isinstance(args_dict[k], int):
                args_dict[k]= int(config_dict[k])
            elif isinstance(args_dict[k], float):
                args_dict[k]= float(config_dict[k])
        self.__dict__= args_dict

