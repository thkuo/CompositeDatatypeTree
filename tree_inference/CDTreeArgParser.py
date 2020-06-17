class CDTreeArgParser:
    def __init__(self):
        ## let the user opt the function
        import argparse
        parser = argparse.ArgumentParser(description=('''
                             Tree inference using composite datatype of microbial representations
                                                     '''))
        parser.add_argument('project_dir', type=str, 
            help='the directory for project')
        parser.add_argument('list_f', type=str, 
            help='the list of sequencing reads')
        parser.add_argument('ref', type=str, 
            help='the reference genome for mapping sequencing reads')
        parser.add_argument('f', type=str,
            help='the function to launch',
            choices=['mapping', 'nuc_tr', 'cd_tr', 'gpa', 'denovo', 'all'])
        parser.add_argument('--config',dest= 'config_f', type=str, 
            help='''
                yaml file to overwrite default parameter settings
            ''')
        parser.add_argument('--cpu',dest= 'cpu', type=int, default= 1,
                help='cpu number; default: %(default)s')
        parser.add_argument('--dry',dest= 'dryrun', action= 'store_true',
                            default= False, 
                help='display the processes and exit')
        self.parser= parser
    def parse(self):
        return(self.parser.parse_args())
