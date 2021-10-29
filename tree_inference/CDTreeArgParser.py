# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

class CDTreeArgParser:
    def __init__(self):
        ## let the user opt the function
        import argparse
        import textwrap
        desc_text= '''\
    Tree inference using composite datatype of microbial representations
    ---
    mapping: mapping the paired-end reads to reference
    with bwa-backtrack and Stampy.py; variants detected
    with freebayes (main output: multi-samples vcf)

    fast_mapping: mapping workflow for paired-end reads
    that are longer than 100bp with bwa-mem; the
    others keep the same as those in the mapping
    function

    nuc_tr: creating the nucleotide alignment from the vcf file
    and computing the phylogenetic tree with RAxML (main
    output: nucleotide tree)

    col_tr: collapsing the nucleotide tree with the log-likelihood
    method (main output: multifurcating nucleotide tree)

    denovo: computing de novo assemblies with spades.py;
    coding regions detection with prokka; orthologous
    clustering with roary (main output: orthologous families
    matrix)

    gpa: encoding the orthologous families table into binary
    states alignment of gene presence/absence (main
    output: gpa alignment)

    cd_tr: inferring composite datatype tree with RAxML (main
    output: composite datatype tree
    ---'''
        parser = argparse.ArgumentParser(
            description= textwrap.dedent(desc_text),
            prog='cdtree', usage='%(prog)s [options]',
            formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('project_dir', type=str,
            help='the directory for project')
        parser.add_argument('list_f', type=str,
            help='the list of paired-end DNA sequencing reads')
        parser.add_argument('ref', type=str,
            help='the reference genome for read mapping')
        parser.add_argument('f', type=str,
            help='the workflow(s) to launch',nargs= '*',
            choices=['mapping', 'fast_mapping', 'nuc_tr', 'col_tr', 'denovo', 'gpa', 'cd_tr', 'all'])
        parser.add_argument('--config',dest= 'config_f', type=str,
            help='yaml file to overwrite default parameter settings')
        parser.add_argument('--cpu',dest= 'cpu', type=int, default= 1,
                help='cpu number; default: %(default)s')
        parser.add_argument('--dry',dest= 'dryrun', action= 'store_true',
                            default= False,
                help='display the processes and exit')
        self.parser= parser
    def parse(self):
        return(self.parser.parse_args())
