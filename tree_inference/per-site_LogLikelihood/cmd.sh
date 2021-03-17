source activate cdtree_env
export cpu_num=20
export config_f=config.yml

snakemake -pn \
  --rerun-incomplete --use-conda \
 --conda-prefix /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs \
 -j  $cpu_num \
 --configfile=$config_f \
 --snakefile=./compute_psll.smk 
