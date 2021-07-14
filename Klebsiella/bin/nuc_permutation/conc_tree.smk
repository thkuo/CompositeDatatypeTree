# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

import os

nuc_tr= config['nuc_tr']
nuc_aln= config['nuc_aln']
gpa_aln= config['gpa_aln']

tr_num= config['tr_num']
sample_ixs= [str(n) for n in range(1, (tr_num+1))]

raxml_bin= config['raxml_bin'] 
raxml_model=config['raxml_model']
max_cores= config['max_cores']


rule for_cw_map_bs_values_to_tree:
    input:  
        cw_all_bs_trees='{result_dir}/bootstrap_gpa_trees.nwk',
        nuc_tr= nuc_tr
    output:
        cw_withBS_tree='{result_dir}/RAxML_bipartitions.{suffix}'
    params:
        raxml_w_dir= '{result_dir}',
        raxml_bin=raxml_bin,
        tree_full_suffix='{suffix}',
        raxml_model= raxml_model,
    threads:1
    shell:
        """
        {params.raxml_bin} \
-f b \
-t {input.nuc_tr} \
-z {input.cw_all_bs_trees} \
-w {params.raxml_w_dir} \
-m {params.raxml_model} \
-n {params.tree_full_suffix}
        """ 

rule for_cw_collect_bs_trees:
    input:
        bs_trees=expand('{{result_dir}}/set_{sample_ix}/RAxML_bestTree.gpa',
            sample_ix= sample_ixs)
    output:
        cw_all_bs_trees='{result_dir}/bootstrap_gpa_trees.nwk'
    shell:
        """
        cat {input.bs_trees} >> {output.cw_all_bs_trees}
        """

rule compute_gpa_tree:
    input:
        perm_gpa_aln= '{result_dir}/set_{sample_ix}/gpa_perm.aln'
    output:
        conc_best_tr='{result_dir}/set_{sample_ix}/RAxML_bestTree.gpa'
    params:
        raxml_suffix= 'gpa',
        raxml_w_dir= '{result_dir}/set_{sample_ix}',
        raxml_bin=raxml_bin,
        raxml_model= raxml_model,
        raxml_starting_num= 1
    threads: max_cores
    shell:
        """
        {params.raxml_bin} -T {threads} \
 -m {params.raxml_model} \
 -s {input.perm_gpa_aln} \
 -p 1 \
 -# {params.raxml_starting_num} \
 -w {params.raxml_w_dir} \
 -n {params.raxml_suffix}
        """

rule create_perm_gpa_aln:
    input:
        nuc_aln= nuc_aln,
        gpa_aln= gpa_aln
    output:
        perm_gpa_aln= '{result_dir}/set_{sample_ix}/gpa_perm.aln'
    run:
        from Bio import AlignIO
        import random
        gpa= AlignIO.read(input['gpa_aln'], 'fasta')
        gpa_len= gpa.get_alignment_length() 
        nuc= AlignIO.read(input['nuc_aln'], 'fasta')
        sample_size= nuc.get_alignment_length()
        perm_gpa= gpa[:, 0:1]
        random.seed(int(wildcards['sample_ix']))
        for k in random.choices(range(gpa_len), k=sample_size):
            perm_gpa= perm_gpa+ gpa[:, k:(k+1)]
        perm_gpa= perm_gpa[:, 1:]
        with open(output['perm_gpa_aln'], 'w') as out_fh:
            AlignIO.write(perm_gpa, out_fh, 'fasta')
