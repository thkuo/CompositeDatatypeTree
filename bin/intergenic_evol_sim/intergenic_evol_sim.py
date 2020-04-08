import os 
import subprocess
import re 
from Bio import Phylo

cpu_num= 10
output_dir=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
            'transmission_simulator/results/intergenic_evol_sim/')
gen_output_dir=('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
                'transmission_simulator/results/genome_evol_sim')
wtrees_dir= ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/'
             'transmission_simulator/results/outbreak_sim/wtrees/')
config_template= './dawg-input_template.dawg' 
intergenic_coor_f= './intergenic_coordinates.txt'
starting_ix= 57

## generation 0
project_name= str(starting_ix)
tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))

# make the config file of dawg
def make_dawg_config(config_template, project_name, tree_f, output_dir):
    config_file= os.path.join(output_dir, '{}.dawg'.format(project_name))
    # read the tree
    tr_str= ''
    with open(tree_f, 'r') as tree_fh:
        tr_str= tree_fh.readline()
        tr_str= tr_str.strip()
    # fill the template
    with open(config_template, 'r') as config_template_h:
        with open(config_file, 'w') as config_fh:
            for l in config_template_h:
                if re.search('\$TREE\$', l):
                    l= re.sub('\$TREE\$', tr_str, l)
                config_fh.write(l)
    return(config_file)
dawg_config_file= make_dawg_config(config_template, project_name, tree_f,
                                   output_dir)

# run dawg
subprocess.run(['./run_dawg.sh', output_dir, gen_output_dir, project_name,
                dawg_config_file, intergenic_coor_f, str(cpu_num)]) 

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

    # make dawg config file
    dawg_config_file= make_dawg_config(config_template, project_name, tree_f, output_dir)

    # run dawg
    subprocess.run(['./run_dawg.sh', output_dir, gen_output_dir, project_name,
                    dawg_config_file, intergenic_coor_f, str(cpu_num)])

    # update the next steps
    tr= Phylo.read(tree_f,'newick')
    next_projects= next_projects + [str(tip.name) for tip in tr.get_terminals()]
    source_projects= source_projects + [project_name] * len(tr.get_terminals())

