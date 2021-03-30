library(seqinr)
library(phytools)

tr_f<- '../../data.v5/outbreak_sim/phyloFromPtree.nwk'
out_f<- '../../Benchmarking/data/gpa.fa'
tr<- read.newick(tr_f)
gain<- 1.0
loss<- 0.12*gain
Q<-matrix(c(-loss,loss, gain, -gain),
	  2,2)
rownames(Q)<-colnames(Q)<-c("1","0")
sim_outs<- sim.history(tree= tr,Q= Q, anc= setNames(c(0.2, 0.8), c('1', '0')),nsim= 2e4)
#summary(sim_outs)
gpa_mat<- sapply(1:length(sim_outs), function(x) sim_outs[[x]]$state)
rownames(gpa_mat)<- names(sim_outs[[1]]$state)
gpa<- as.list(apply(gpa_mat, 1, function(r) paste0(r, collapse='')))
write.fasta(gpa, names(gpa), file.out= out_f, nbchar= 70)

