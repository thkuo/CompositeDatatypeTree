
# CHANGELOG
# v4: simulate only variant sites
# v3: rate heterogeneity like how alf does

library(argparse)
library(seqinr)
library(phytools)
library(ape)


rescale_tr<- function(tr,g.scale= 1, g.shape= 1){
  br.scale<- rgamma(n=1, shape= g.shape, scale= g.scale)
  new_tr<- tr
  new_tr$edge.length <- tr$edge.length * br.scale 
  return(new_tr)
}

run_sim<- function(tree, Q , anc, nsim= 1, cutoff=0){
  is_var<- F
  sim_out<- list()
  while (!is_var){
    sim_out<- sim.history(tree= tree,
			  Q= Q,
			  anc= anc,
			  nsim=nsim, 
			  message=F)
    states<- sim_out$state
    if ((length(unique(states)) > 1) & (max(table(states)) < (length(states)-cutoff))){
      is_var<- T
    }

  }
  return(sim_out)
}


set.seed(snakemake@params[['p']])
tr_f<-snakemake@input[['tr_f']]
out_f<-snakemake@output[[1]]
gain<-snakemake@params[['gain_r']]
loss<-snakemake@params[['loss_r']]
sim_num<-snakemake@params[['sim_num']]
a0<-snakemake@params[['a0']]
p0<- 1-a0
cutoff<-snakemake@params[['cutoff']]

tr<- read.newick(tr_f)
Q<-matrix(c(-loss,loss, gain, -gain),
	  2,2)
rownames(Q)<-colnames(Q)<-c("1","0")

total<- sim_num
sim_outs<- lapply(1:sim_num, function(n) run_sim(tree= rescale_tr(tr),
						 Q= Q,
						 anc= setNames(c(p0, a0), c('1', '0')),
						 nsim= 1,
						 cutoff=cutoff))
gpa_mat<- sapply(1:length(sim_outs), function(x) sim_outs[[x]]$state)
rownames(gpa_mat)<- names(sim_outs[[1]]$state)
gpa<- as.list(apply(gpa_mat, 1, function(r) paste0(r, collapse='')))
write.fasta(gpa, names(gpa), file.out= out_f, nbchar= 70)

