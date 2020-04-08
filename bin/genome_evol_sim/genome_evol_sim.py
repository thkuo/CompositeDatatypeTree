import os 
import subprocess
import re 
from Bio import Phylo

output_dir='/net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator/results/genome_evol_sim/'

## generation 0
genome_db_f= ('/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/'
              'WhichTree_Sim/bin/re_simulate/which_tree/Streptococcus_pneumoniae_ATCC_700669_v1.db')
starting_ix= 57
project_name= str(starting_ix)
wtrees_dir= ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator/'
             'results/outbreak_sim/wtrees/')
tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))
config_template= 'alf-input_template.drw' 

# make the config file of ALF
def make_alf_config(config_template, project_name, genome_db_f, tree_f, output_dir):
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
alf_config_file= make_alf_config(config_template, project_name, genome_db_f, tree_f, output_dir)
# run alf
subprocess.run(['alfsim', '-o', output_dir, alf_config_file])

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
    subprocess.run(['alfsim', '-o', output_dir, alf_config_file])

    # update the next steps
    tr= Phylo.read(tree_f,'newick')
    next_projects= next_projects + [str(tip.name) for tip in tr.get_terminals()]
    source_projects= source_projects + [project_name] * len(tr.get_terminals())

