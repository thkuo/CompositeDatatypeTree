genome simulator + transmission simulation
- within-host evolution of nucleotide
- within-host evolution through LGT
- incomplete sampling 
- output: genomes, transmission tree, and time-measured tree of sampled cases

## Idea
- Dependencies
    - transphylo: big tree where each node is a small coalescence tree
    - alfsim: 
	- fixed parameters (mutation rates, substitution model...etc)
	- variables: inital genome and the tree
- transphylo+alfsim:
    stage 0. [TransPhylo] Create overall transmission tree and the within-host trees 
    stage 1. Travel along the tree from the root (generation 0)
	1. determine the tree topology
	2. determine the initial genome, which is the starting genome for the origin and a simulated genome for the others
	3. [ALF] simulate the within-host evolution with (1) and (2), and obtain the genomes at tips
	    i. create the config files
	    ii. run alfsim for simulating the coding regions
	    iii. run dawg for simulating the integenic regions
	4. map the simulated genomes, go to the offsprings, and repeat (1)-(3)

## Challenges
- optimize the substitution rate (iteratively?)
- nodes matrix (transphylo) -> time-mesured tree
- how does it differ from alfsim with the whole time-measured tree?
    For the divergence of host environment (e.g. mutation rate fast in this host, gene pool
different in another...)
