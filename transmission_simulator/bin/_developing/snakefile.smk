
rule all:
    ## given a set of within-host trees
    ## create the cimulation 

rule sim_wtrees:
    input:
        sim_config=
    output:
        wtrees= dynamic('{wtree_dir}/{id}.nwk')
    script: 

rule sim_coding_region:
    input:
        wtree=
        alf_config= 
    output: 
        gene_seqs= 

rule sim_intergenic_region:
    input:
        dawg_config_template=
        intergenic_coor=
        wtree=
    output:
        dawg_config='{in_outdir}/{project}.dawg'
        intergenic_seq='{in_outdir}/{project}/intergenic.fa'
        genome_seq='{in_outdir}'
    threads:
    params:
        alf_db_to_fasta
    conda:
    shell:
        '''
        IN_SIM_OUT="$1"
        GEN_SIM_OUT="$2"
        PROJECT="$3"
        DAWG_CONFIG_FILE="$4"
        INTERGENIC_COOR_FILE="$5"
        CPU_NUM={threads}
        
        #++ use the built-in conda function
        #Using dawg for intergenic regions, with parameter file intergenic.dawg
        #source activate whichtree_env 
        PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2/:$PERL5LIB
        #++
        
        mkdir -p $IN_SIM_OUT/$PROJECT
        dawg --output {output.intergenic_seq} \
        $DAWG_CONFIG_FILE
        $WHICH_TREE_HOME/dawg_to_alf.pl $GEN_SIM_OUT/$PROJECT/speciesMapping.txt \
        $INTERGENIC_COOR_FILE \
        $IN_SIM_OUT/$PROJECT/intergenic.fa
        mkdir -p $IN_SIM_OUT/$PROJECT/MSA/
        mv MSA_i* $IN_SIM_OUT/$PROJECT/MSA/
        #combine with the alf results
        #creating the genome sequences
        parallel -j $CPU_NUM "$WHICH_TREE_HOME/alf_db_to_fasta_splice_intergenic.pl $GEN_SIM_OUT/$PROJECT/DB/{}_dna.fa $IN_SIM_OUT/$PROJECT/MSA | sed 's/>.\+\(DB\/SE[0-9]\+_dna.fa\)/>\1/'> $GEN_SIM_OUT/$PROJECT/DB/{}_genome.fa" \
  ::: `ls $GEN_SIM_OUT/$PROJECT/DB| grep dna.fa| sed 's/_dna.fa//'`
         '''


rule sim_reads_for_tips:
    
    input:
        genome_sim_dir=
        
    output:
        genomes_files_dict_f= 
        reads=  
    thread:
        
    run:
        import subprocess
        subprocess.run(['pirs', 'simulate', '-t', str(thread), 
                        '-l', '100', '-x', '100', '-m', '250', 
                        '-z', '-o', os.path.join(reads_dir, ix),
                        input.genome_seq])
