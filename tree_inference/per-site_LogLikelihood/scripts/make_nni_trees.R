library(phytools)
library(phangorn)

main_tr_f<- snakemake@input[['nuc_tr']]
print(sprintf('input tree %s', main_tr_f))
main_tr<- read.newick(main_tr_f) 
out_f<-snakemake@output[['rear_trs']]
rear_trs<- nni(main_tr)
write.tree(c(main_tr, rear_trs), out_f)
