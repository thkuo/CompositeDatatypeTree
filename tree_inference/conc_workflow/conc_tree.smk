import os

nuc_pslls= config['nuc_pslls']
rear_nuc_tr=config['rear_nuc_trs']
gpa_aln= config['gpa_aln']
nuc_aln= config['nuc_aln']
model= config['model']
output_dir= config['output_dir']
psll_percs= config['psll_percs']
norm_by_smallest= config['norm_weight'] if 'norm_weight' in config else False
weighted_conc= config['weighted'] if 'weighted' in config else False
tr_name= 'conc' if not weighted_conc else 'conc_weighted'

rule all:
    input:
        conc_best_tr=expand('{result_dir}/{psll_perc}/RAxML_bestTree.{tr_name}',
            result_dir= output_dir, psll_perc=psll_percs, tr_name= tr_name)

#rule for_cw_map_bs_values_to_tree:
#    input:  
#        cw_all_bs_trees='{result_dir}/{psll_perc}/bootstrap/RAxML_bootstrap.all',
#        conc_best_tr='{result_dir}/{psll_perc}/RAxML_bestTree.{suffix}'
#    output:
#        cw_withBS_tree= '{result_dir}/{psll_perc}/RAxML_bipartitions.{suffix}.bs'
#    params:
##        raxml_bin=raxml_bin,
#        tree_full_suffix='{suffix}.bs',
#        raxml_model= model
#    threads:1
#    conda:'/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
#    shell:
#        """
#        export RAXML_BIN='raxmlHPC-PTHREADS'
#        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
#            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'
#        fi
#        $RAXML_BIN \
#-f b \
#-t {input.conc_best_tr} \
#-z {input.cw_all_bs_trees} \
#-w {wildcards.result_dir} \
#-m {params.raxml_model} \
#-n {params.tree_full_suffix}
#        """ 
#
#rule for_cw_collect_bs_trees:
#    input:
#        bs_trees=lambda wildcards:[
#        os.path.join(wildcards.result_dir,'bootstrap',
#        'RAxML_bootstrap.{}.20'.format(part_num)) for part_num in range(1, 5+1)]
#    output:
#        cw_all_bs_trees='{result_dir}/{psll_perc}/bootstrap/RAxML_bootstrap.all'
#    shell:
#        """
#        cat {input.bs_trees} > {output.cw_all_bs_trees}
#        """
#
#rule bootstrap:
#    input: 
#        tr='{result_dir}/{psll_perc}/raxml_input/nuc_col.nwk',
#        raxml_input_aln='{result_dir}/{psll_perc}/raxml_input/conc.aln',
#        raxml_input_partitions='{result_dir}/{psll_perc}/raxml_input/conc.prtn'
#    output:
#        bs_tree='{result_dir}/{psll_perc}/bootstrap/RAxML_bootstrap.{part_num}.{per_part_tr_num}'
#    params:
#        bs_dir= '{result_dir}/{psll_perc}/bootstrap',
#        raxml_model= model,
#        raxml_starting_num= 1,
#    conda:'/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
#    threads: 8
#    shell:
#        """
#        export RAXML_BIN='raxmlHPC-PTHREADS'
#        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
#            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'     
#        fi
#        $RAXML_BIN \
#   -T {threads} \
#   -g {input.tr} \
#   -q {input.raxml_input_partitions} \
#   -m {params.raxml_model} \
#   -s {input.raxml_input_aln} \
#   -n {wildcards.part_num}.{wildcards.per_part_tr_num} \
#   -p {wildcards.part_num} -x {wildcards.part_num} \
#   -# {wildcards.per_part_tr_num} \
#   -w {params.bs_dir}
#        """

rule compute_conc_tree:
    input:
        tr='{result_dir}/{psll_perc}/raxml_input/nuc_col.nwk',
        raxml_input_aln='{result_dir}/{psll_perc}/raxml_input/conc.aln',
        raxml_input_partitions='{result_dir}/{psll_perc}/raxml_input/conc.prtn'
    output:
        conc_best_tr=('{result_dir}/{psll_perc}/'
            'RAxML_bestTree.conc')
    params:
#        raxml_bin=raxml_bin,
        tr_name= 'conc',
        raxml_wd= '{result_dir}/{psll_perc}',
        raxml_model= model,
        raxml_starting_num= 1
    threads: 16 
    conda:'/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'     
        fi
        $RAXML_BIN \
  -T {threads} \
  -g {input.tr} \
  -q {input.raxml_input_partitions} \
  -m {params.raxml_model} \
  -s {input.raxml_input_aln} \
  -p 1 \
  -# {params.raxml_starting_num} \
  -w {params.raxml_wd} \
  -n {params.tr_name}
        """


rule compute_conc_tree_weighted:
    input:
        tr='{result_dir}/{psll_perc}/raxml_input/nuc_col.nwk',
        raxml_input_aln='{result_dir}/{psll_perc}/raxml_input/conc.aln',
        raxml_input_partitions='{result_dir}/{psll_perc}/raxml_input/conc.prtn',
        raxml_input_weights='{result_dir}/{psll_perc}/raxml_input/conc.weights'
    output:
        conc_best_tr=('{result_dir}/{psll_perc}/'
            'RAxML_bestTree.conc_weighted')
    params:
        tr_name= 'conc_weighted',
        raxml_wd= '{result_dir}/{psll_perc}',
        raxml_model= model,
        raxml_starting_num= 1
    threads: 16 
    conda:'/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs/raxml_env.yml'
    shell:
        """
        export RAXML_BIN='raxmlHPC-PTHREADS'
        if [ $(grep avx2 /proc/cpuinfo|wc -l) -gt 0 ] ; then
            export RAXML_BIN='raxmlHPC-PTHREADS-AVX2'     
        fi
        $RAXML_BIN \
  -T {threads} \
  -a {input.raxml_input_weights} \
  -g {input.tr} \
  -q {input.raxml_input_partitions} \
  -m {params.raxml_model} \
  -s {input.raxml_input_aln} \
  -p 1 \
  -# {params.raxml_starting_num} \
  -w {params.raxml_wd} \
  -n {params.tr_name}
        """

rule prepare_concatenated_data:
    input:
        nuc_aln= nuc_aln,
        gpa_aln= gpa_aln
    output:
        raxml_input_aln='{result_dir}/{psll_perc}/raxml_input/conc.aln',
        raxml_input_partitions='{result_dir}/{psll_perc}/raxml_input/conc.prtn'
    run:
        from Bio import AlignIO
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        import sys

        def convert_alphabet(aln, alphabet):
            new_aln= aln
            for a in new_aln:
                a.seq.alphabet=alphabet
            return(new_aln)

        n_aln= AlignIO.read(input.nuc_aln, 'fasta')
        g_aln= AlignIO.read(input.gpa_aln, 'fasta')
        n_aln_l= n_aln.get_alignment_length()
        g_aln_l= g_aln.get_alignment_length()
        print('nucleotide sites: {}'.format(str( n_aln_l)))
        print('gpa sites: {}'.format(str( g_aln_l)))
        # test alignments
        try:
            i=set([str(a.id) for a in n_aln]) & set([str(a.id) for a in g_aln])
            u=set([str(a.id) for a in n_aln]) | set([str(a.id) for a in g_aln])
            if len(u-i) > 0:
                raise Exception('The two alignments include different samples')
        except Exception as e:
            sys.exit(str(e))
        # create the alignment
        n_aln.sort()
        n_aln=convert_alphabet(n_aln, generic_dna)
        g_aln.sort()
        g_aln=convert_alphabet(g_aln, generic_dna)
        conc_aln= n_aln+g_aln
        AlignIO.write(conc_aln, output.raxml_input_aln, 'fasta')
        # print the partitions
        with open(output.raxml_input_partitions, 'w') as out_fh:
            out_fh.write('DNA, nuc= 1-{}\n'.format(str(n_aln_l)))
            out_fh.write('BIN, gpa= {}-{}\n'.format(
                str(n_aln_l+1), str(n_aln_l+g_aln_l)
            ))

rule prepare_concatenated_weight:
    input:
        nuc_aln= nuc_aln,
        gpa_aln= gpa_aln
    output:
        raxml_input_weights='{result_dir}/{psll_perc}/raxml_input/conc.weights'
    params:
        norm_by_smallest= norm_by_smallest
    run:
        from Bio import AlignIO
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        import sys
        import argparse

        def convert_alphabet(aln, alphabet):
            new_aln= aln
            for a in new_aln:
                a.seq.alphabet=alphabet
            return(new_aln)

        n_aln= AlignIO.read(input.nuc_aln, 'fasta')
        g_aln= AlignIO.read(input.gpa_aln, 'fasta')
        n_aln_l= n_aln.get_alignment_length()
        g_aln_l= g_aln.get_alignment_length()

        n_weights= [g_aln_l] * n_aln_l
        g_weights= [n_aln_l] * g_aln_l
        all_weights= n_weights + g_weights
        if params.norm_by_smallest:
            import math
            min_weight= min(all_weights)
            all_weights= [math.floor(w/min_weight) for w in all_weights]

        with open(output.raxml_input_weights, 'w') as out_fh:
            out_fh.write(' '.join([str(w) for w in all_weights]))

rule col_tree:
    input:
        nuc_pslls= nuc_pslls,
        rear_trs=rear_nuc_tr
    output:
        col_nuc_tr='{result_dir}/{psll_perc}/raxml_input/nuc_col.nwk'
    threads: 8
    params: 
        ll_perc='{psll_perc}'
    script:'./scripts/col_tree_diff_cutoffs.R'
