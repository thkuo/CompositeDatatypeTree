# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule transphylo:
    input:
        max_credential_tr= '{w_dir}/{tr_name}.d_tree.nex',
        time_tab='{w_dir}/time.tab'
    output:
        ctree_rds='{w_dir}/transphylo/cd_ctree.Rds',
        ttree_rds='{w_dir}/transphylo/cd_ttree.Rds'
    params:
        out_d='{w_dir}/transphylo'
    conda: 'env/r35.yml'
    shell: 
        '''
        Rscript ./transmission_ana.R \
  -t {input.max_credential_tr} \
  -d {input.time_tab} \
  -p {wildcards.tr_name} -o {params.out_d}
        '''

 
rule conclude_beast2_result:
    input:
        mcmc_trees= '{w_dir}/{tr_name}.trees'
    output:
        max_credential_tr= '{w_dir}/{tr_name}.d_tree.nex'
    conda: 'env/beast25_env.yml'
    params: 
        burnin=10
    shell:
        '''
        treeannotator -burnin {params.burnin} \
 {input.mcmc_trees} {output.max_credential_tr}
        '''

rule dating_with_beast2:
    input:
        xml_f= '{w_dir}/beast.xml'
    output:
        mcmc_trees= '{w_dir}/{tr_name}.trees'
    conda: 'env/beast25_env.yml'
    threads: 5
    shell:
        '''
        beast -working -seed 1234 -threads {threads} {input.xml_f}
        '''
