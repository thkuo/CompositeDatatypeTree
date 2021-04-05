export out_f=$( realpath ../../results/nuc_rep_by_gpa/RAxML_bipartitions.gpa.bs )
source activate cdtree_env
snakemake -p \
  --restart-times 3 \
  -j 15 \
  --use-conda \
  --configfile=conc.config.yml \
  --snakefile=conc_tree.smk \
  $out_f
