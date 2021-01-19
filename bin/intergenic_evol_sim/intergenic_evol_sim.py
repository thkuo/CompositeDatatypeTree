import os 
import subprocess
import re 
from Bio import Phylo

#+++++
# Core functions
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

def exec_func(output_dir, gen_output_dir, project_name,
              dawg_config_file, intergenic_coor_f, cpu_num):
    # run dawg
    try:
        t= 0
        success=False
        while (not success) and (t<3):
            if os.path.isdir(os.path.join(output_dir, project_name)):
                shutil.rmtree(os.path.join(output_dir, project_name))
            subprocess.run([
                os.path.join(os.path.dirname(
                    os.path.realpath(__file__)),'run_dawg.sh'), 
                output_dir, gen_output_dir, project_name,
                dawg_config_file, intergenic_coor_f, str(cpu_num)]) 
            success=os.path.isfile(
                os.path.join(output_dir, project_name, 'intergenic.fa'))
            t=t+1
        if (not success):
            raise FileExistsError
    except FileExistsError as fe:
        print(fe)
        print('Host {}: evolution of integenic regions failed'.format(
            project_name))

#+++++
# arguments
import argparse 
import os
parser = argparse.ArgumentParser(
    description='Simulate intergenic regions using dawg')
parser.add_argument('-o',dest= 'output_dir', type=str,
                    required= True,
                    help='output directory')
parser.add_argument('-gd', dest='genome_dir', type= str,
                    required= True,
                    help='the directory of simulated genomes')
parser.add_argument('-wd', dest='wtrees_dir', type= str,
                    required= True,
                    help='the directory of wtrees')
parser.add_argument('-c', dest='conf',type= str,
                    required= True,
                    help='config file template of dawg')
parser.add_argument('-co', dest='coor',type= str, 
                    required= True,
                    help='intergenic region coordinates')
parser.add_argument('-n', dest='cores',type= int,
                    default= 1,
                    help='number of threads')

args = parser.parse_args()

os.chdir(os.path.dirname(os.path.abspath(__file__)))
output_dir=args.output_dir 
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

cpu_num= args.cores
gen_output_dir=args.genome_dir
assert os.path.isdir(gen_output_dir)
wtrees_dir=args.wtrees_dir 
assert os.path.isdir(wtrees_dir)
config_template= args.conf 
assert os.path.isfile(config_template)
intergenic_coor_f= args.coor 
# determine the starting index
starting_ix= max([int(re.sub('.nwk$', '', f)) for f in os.listdir(wtrees_dir) if
    re.search('.nwk$', f)])

## generation 0
project_name= str(starting_ix)
tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))

# make the config file of dawg
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

