closeAllConnections()
rm(list=ls())

library(argparse)
library(seqinr)
library(phytools)

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

tr<- read.newick(tr_f)
Q<-matrix(c(-loss,loss, gain, -gain),
	  2,2)
rownames(Q)<-colnames(Q)<-c("1","0")
sim_outs<- sim.history(tree= tr,Q= Q, anc= setNames(c(p0, a0), c('1', '0')),nsim= sim_num)
#summary(sim_outs)
gpa_mat<- sapply(1:length(sim_outs), function(x) sim_outs[[x]]$state)
rownames(gpa_mat)<- names(sim_outs[[1]]$state)
gpa<- as.list(apply(gpa_mat, 1, function(r) paste0(r, collapse='')))
write.fasta(gpa, names(gpa), file.out= out_f, nbchar= 70)

