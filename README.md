genome simulator + transmission simulation
- within-host evolution of nucleotide
- within-host evolution through LGT
- incomplete sampling 
- output: genomes, transmission tree, and time-measured tree of sampled cases

## Idea
- transphylo: big tree where each node is a small coalescence tree
- alfsim: 
    - fixed parameters (mutation rates, substitution model...etc)
    - variables: inital genome and the tree
- transphylo+alfsim:
Create overall transmission tree and the within-host trees
for the root (G0), 
    1. extract the within-host tree
    2. simulate with the initial genome and the within-host tree
    3. extract the genome of offsprings
    4. go to the offspring (breadth-first traversal) and repeat 1-3 until tips

## Challenges
- optimize the substitution rate (iteratively?)
- nodes matrix (transphylo) -> time-mesured tree
- how does it differ from alfsim with the whole time-measured tree?
    For the divergence of host environment (e.g. mutation rate fast in this host, gene pool
different in another...)
