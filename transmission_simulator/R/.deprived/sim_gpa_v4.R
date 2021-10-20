closeAllConnections()
rm(list=ls())

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


set.seed(0)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
    dest="verbose")
parser$add_argument("-t", "--tree", dest= 'tr_f', 
		    help="the tree topology")
parser$add_argument("-gr", "--gain_rate", dest= 'grate', nargs= 1,
		    type="double",
		    default= 1.0, help="the rate of gene gain")
parser$add_argument("-lr", "--loss_rate", dest= 'lrate', nargs= 1,
		    type="double",
		    default= 0.12, help="the rate of gene loss")
parser$add_argument("-c", dest= 'cutoff', nargs= 1,
		    type="integer",
		    default= 2, help="at least <num> samples have a state different from the majority")
parser$add_argument("-l", "--length", dest= 'sim_num', nargs= 1,
		    type="integer",
		    default= 2e4, help="the length of simulated patterns")
parser$add_argument("-a0", "--init_freq_absence", dest= 'a0', nargs= 1,
		    type="double",
		    default= 0.8, help="the probability of initial absence (range: [0.0, 1.0])")
parser$add_argument("-o", "--output_gpa", dest= 'out',
    help="output filename of gene presence or absence")
args <- parser$parse_args()
tr_f<- args$tr_f
out_f<- args$out
gain<- args$grate
loss<- args$lrate 
sim_num<- args$sim_num
a0<- args$a0
p0<- 1-a0
cutoff<- args$cutoff

tr<- read.newick(tr_f)
Q<-matrix(c(-loss,loss, gain, -gain),
	  2,2)
rownames(Q)<-colnames(Q)<-c("1","0")

total<- sim_num
sim_outs<- list()
for (n in 1:sim_num){
  sim_outs[[n]]<- run_sim(tree= rescale_tr(tr),
			  Q= Q,
			  anc= setNames(c(p0, a0), c('1', '0')),
			  nsim= 1,
			  cutoff=cutoff)
  print(n)

}
gpa_mat<- sapply(1:length(sim_outs), function(x) sim_outs[[x]]$state)
rownames(gpa_mat)<- names(sim_outs[[1]]$state)
gpa<- as.list(apply(gpa_mat, 1, function(r) paste0(r, collapse='')))
write.fasta(gpa, names(gpa), file.out= out_f, nbchar= 70)

