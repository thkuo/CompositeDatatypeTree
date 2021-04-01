
#export d_tree=../../results.5/transmission_ana/cd.best/clade-wise/4/v_best/cd.d_tree.nex
#export time_table_f=../../results.5/transmission_ana/cd.best/clade-wise/4/time.tab
#export tree_name='cd'
#export out_dir=../../results.5/transmission_ana/cd.best/clade-wise/4/v_best/transphylo

if [ $# -eq 0 ]; then
  echo 'Usage: ./run_transphylo.sh [d_tree] [time_table_f] [tree_name] [out_dir]'
  exit
fi

export d_tree=$1
export time_table_f=$2
export tree_name=$3
export out_dir=$4

source activate r35
Rscript ./transmission_ana.R \
  -t $d_tree \
  -d $time_table_f \
  -p $tree_name -o $out_dir 
source deactivate


