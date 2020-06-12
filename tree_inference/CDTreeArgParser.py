class CDTreeArgParser:
    def __init__(self):
        ## let the user opt the function
        import argparse
        parser = argparse.ArgumentParser(description=('''
                             Tree inference using composite datatype of microbial representations
                                                     '''))
        parser.add_argument('--p', dest= 'project_dir', type=str, 
                required= True,
                help='the directory for project')
        parser.add_argument('--l',dest= 'list_f', type=str, 
                required= True,
                help='the list of sequencing reads')
        parser.add_argument('--r',dest= 'ref', type=str, 
                required= True,
                help='the reference genome for mapping sequencing reads')
        parser.add_argument('f', type=str,
                help='the function to launch',
                choices=['mapping', 'nuc_tr', 'cd_tr', 'gpa', 'denovo'])
        parser.add_argument('--cpu',dest= 'cpu', type=str, default= 1,
                help='cpu number; default: %(default)s')
        parser.add_argument('--dry',dest= 'dryrun', action= 'store_true',
                help='display the processes and exit')
        parser.add_argument('--nuc_subs_m', dest= 'nuc_subs_m', 
                default= 'GTR',
                help=('substitution model for '
                'nucleotide tree inference '
                '(same as those for RAxML); '
                'default: %(default)s'))
        parser.add_argument('--rate_m', dest= 'rate_m', 
                default= 'GAMMA',
                help=('rate variation model for '
                'both nucleotide and composite datatype tree inference '
                '(same as those for RAxML); '
                'default: %(default)s'))
        self.parser= parser
    def parse(self):
        return(self.parser.parse_args())
