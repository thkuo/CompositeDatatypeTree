import os 
import subprocess
import re 
from Bio import Phylo


#+++++
# arguments
import argparse

parser = argparse.ArgumentParser(
    description='Simulate sequencing reads using PIRS')
parser.add_argument('-o',dest= 'output_dir', type=str,
                    required= True,
                    help='output directory')
parser.add_argument('-gd', dest='genome_dir', type= str,
                    required= True,
                    help='the directory of simulated genomes')
parser.add_argument('-wd', dest='wtrees_dir', type= str,
                    required= True,
                    help='the directory of wtrees')
parser.add_argument('-n', dest='cores',type= int,
                    default= 1,
                    help='number of threads')
args = parser.parse_args()

#output_dir=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
#            'TreePaper/WhichTree_Sim.v6/data/reads_sim')
output_dir=args.output_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#cpu_num=10
cpu_num= args.cores
#gen_output_dir=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
#            'WhichTree_Sim.v6/data/genome_evol_sim/')
gen_output_dir=args.genome_dir
#wtrees_dir=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
#            'WhichTree_Sim.v6/data/outbreak_sim/wtrees')
wtrees_dir=args.wtrees_dir 
genome_files_mat_f= os.path.join(output_dir, 'genome_paths.tsv')
reads_dir= os.path.join(output_dir, 'fastq')
# determine the starting index
starting_ix= max([int(re.sub('.nwk$', '', f)) for f in os.listdir(wtrees_dir) if
    re.search('.nwk$', f)])

## generation 0
project_name= str(starting_ix)
tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))

## determine who the tips are and where their genomes were made
tr= Phylo.read(tree_f,'newick')
next_projects= [str(tip.name) for tip in tr.get_terminals()]
source_projects= [project_name] * len(tr.get_terminals())
#travel along the trees
genome_files_dict= {}
from datetime import datetime
while (len(next_projects)>0):
    project_name= next_projects.pop(0)
    source_project= source_projects.pop(0) 
#    print(datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
#    print('{}->{}'.format(source_project, project_name))
    tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))
    if project_name == 'pseudo':
        continue
    if not os.path.isfile(tree_f):
        # it is a tip
        # first, know how the previous ALF encode the species names
        species_map_f=os.path.join(gen_output_dir, source_project,
                                   'speciesMapping.txt')
        species= {}
        for l in open(species_map_f, 'r'):
            d= l.strip().split('\t')
            species[d[0]]= d[1]
        alf_id= species[project_name] 
        genome_f= os.path.join(gen_output_dir , source_project,
                                  'DB','{}_genome.fa'.format(alf_id))
        genome_files_dict[project_name]=genome_f
    # update the next steps
    else:
        tr= Phylo.read(tree_f,'newick')
        print([project_name] * len(tr.get_terminals()))
        next_projects= next_projects + [str(tip.name) for tip in tr.get_terminals()]
        source_projects= source_projects + [project_name] * len(tr.get_terminals())

from pprint import pprint
pprint(genome_files_dict)
## simulate the reads
with open(genome_files_mat_f, 'w') as genome_files_mat_fh:
    for ix in genome_files_dict:
        print(ix)
        genome_files_mat_fh.write('{}\t{}\n'.format(ix, genome_files_dict[ix]))
        if not os.path.exists(reads_dir):
            os.mkdir(reads_dir)
        subprocess.run(['pirs', 'simulate', '-t', str(cpu_num), 
                        '-l', '100', '-x', '100', '-m', '250', 
                        '-z', '-o', os.path.join(reads_dir, ix),
                        genome_files_dict[ix]])

