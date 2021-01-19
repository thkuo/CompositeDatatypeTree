import os 
import subprocess
import re 
from Bio import Phylo
import shutil

#+++++
# Core functions
def make_alf_config(config_template, project_name, genome_db_f, tree_f, output_dir):
    # create the ALF config file using the template file
    config_file= os.path.join(output_dir, '{}.drw'.format(project_name))
    with open(config_template, 'r') as config_template_h:
        with open(config_file, 'w') as config_fh:
            for l in config_template_h:
                if re.search('\$ORGANISM_DB\$', l):
                    l= re.sub('\$ORGANISM_DB\$', genome_db_f, l)
                elif re.search('\$TREE_FILE\$', l):
                    l= re.sub('\$TREE_FILE\$', tree_f, l)
                elif re.search('\$PROJECT_NAME\$', l):
                    l= re.sub('\$PROJECT_NAME\$', project_name, l)
                config_fh.write(l)
    return(config_file)

def exec_func(output_dir, alf_config_file, project_name):
    # run ALF
    try:
        t= 0
        success=os.path.isfile(
                os.path.join(output_dir, project_name, 'DB', 'SE002_dna.fa'))
        while (not success) and (t<3):
            if os.path.isdir(os.path.join(output_dir, project_name)):
                shutil.rmtree(os.path.join(output_dir, project_name), ignore_errors=True)
            subprocess.run(['alfsim', '-o', output_dir, alf_config_file])
            success=os.path.isfile(
                os.path.join(output_dir, project_name, 'DB', 'SE002_dna.fa'))
            t=t+1
        if (not success):
            raise FileExistsError
    except FileExistsError as fe:
        print(fe)
        print('Host {}: genomic evolution unsuccessfully'
              'simulated'.format(project_name))

#+++++
# arguments
import argparse

parser = argparse.ArgumentParser(description='Simulate coding regions using ALF')
parser.add_argument('-o',dest= 'output_dir', type=str,
                    required= True,
                    help='output directory')
parser.add_argument('-g', dest='genome_db', 
                    required= True,
                    help='the initial genome')
parser.add_argument('-wd', dest='wtrees_dir', 
                    required= True,
                    help='the directory of wtrees')
parser.add_argument('-c', dest='conf', 
                    required= True,
                    help='config file template of ALF')


args = parser.parse_args()
#output_dir=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
#            'WhichTree_Sim.v6/data/genome_evol_sim/')
output_dir=args.output_dir  
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

## generation 0
#genome_db_f= ('/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
#              'WhichTree_Sim/bin/re_simulate/which_tree/Streptococcus_pneumoniae_ATCC_700669_v1.db')
genome_db_f= args.genome_db 
assert os.path.isfile(genome_db_f)
#wtrees_dir=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
#            'WhichTree_Sim.v6/data/outbreak_sim/wtrees')
wtrees_dir=args.wtrees_dir 
assert os.path.isdir(wtrees_dir)
# determine the starting index
starting_ix= max([int(re.sub('.nwk$', '', f)) for f in os.listdir(wtrees_dir) if
    re.search('.nwk$', f)])
project_name= str(starting_ix)
tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))
#config_template= ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
#            'WhichTree_Sim.v6/data/genome_evol_sim/alf-input_template.drw' )
config_template=args.conf
assert os.path.isfile(config_template)

#+++++
# Simulation starts from the fist generation
# make the config file of ALF
alf_config_file= make_alf_config(config_template, project_name, genome_db_f, tree_f, output_dir)
# run alf
exec_func(output_dir, alf_config_file, project_name)

# determine next tree and genome
tr= Phylo.read(tree_f,'newick')
next_projects= [str(tip.name) for tip in tr.get_terminals()]
source_projects= [project_name] * len(tr.get_terminals())
#travel along the trees
from datetime import datetime
while (len(next_projects)>0):
    project_name= next_projects.pop(0)
    source_project= source_projects.pop(0) 
    print(datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
    print('{}->{}'.format(source_project, project_name))
    tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))
    if not os.path.isfile(tree_f):
        # it is a tip
        continue
    # find the db file
    # first, know how the previous ALF encode the species names
    species_map_f=os.path.join(output_dir, source_project,
                               'speciesMapping.txt')
    species= {}
    for l in open(species_map_f, 'r'):
        d= l.strip().split('\t')
        species[d[0]]= d[1]
    alf_id= species[project_name] 
    genome_db_f= os.path.join(output_dir , source_project,
                              'DB','{}.db'.format(alf_id))
    # make alf config file
    alf_config_file= make_alf_config(config_template, project_name, genome_db_f, tree_f, output_dir)

    # run alf
    exec_func(output_dir, alf_config_file, project_name)

    # update the next steps
    tr= Phylo.read(tree_f,'newick')
    next_projects= next_projects + [str(tip.name) for tip in tr.get_terminals()]
    source_projects= source_projects + [project_name] * len(tr.get_terminals())

