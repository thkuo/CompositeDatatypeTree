import os 
import subprocess
import re 
from Bio import Phylo
import argparse 

class simulator:
    def __init__(self, config_template, project_name, source_project, tree_f, 
               wtrees_dir, output_dir, gen_output_dir, intergenic_coor_f):
        self.project_name= project_name
        self.source_project= source_project
        self.tree_f= tree_f
        self.output_dir= output_dir
        self.gen_output_dir= gen_output_dir
        self.config_template= config_template
        self.intergenic_coor_f= intergenic_coor_f
        if not os.path.isdir(wtrees_dir):
            sys.exit('{} not existing'.format(wtrees_dir))
        elif not os.path.isfile(tree_f):
            raise FileExistsError('Tip node')

    def run(self, cpu_num=1):
        print(datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
        print('{}->{}'.format(self.source_project, self.project_name))
        # make dawg config file
        dawg_config_file= self.make_dawg_config()
        # run dawg
        self.exec_func(dawg_config_file, cpu_num)

    #+++++
    # Core functions
    def make_dawg_config(self):
        config_file= os.path.join(self.output_dir, '{}.dawg'.format(self.project_name))
        # read the tree
        tr_str= ''
        with open(self.tree_f, 'r') as tree_fh:
            tr_str= tree_fh.readline()
            tr_str= tr_str.strip()
        # fill the template
        with open(self.config_template, 'r') as config_template_h:
            with open(config_file, 'w') as config_fh:
                for l in config_template_h:
                    if re.search('\$TREE\$', l):
                        l= re.sub('\$TREE\$', tr_str, l)
                    config_fh.write(l)
        return(config_file)

    def exec_func(self,dawg_config_file, cpu_num=1):
        # run dawg
        t= 0
        success=False
        while (not success) and (t<3):
            if os.path.isdir(os.path.join(self.output_dir, self.project_name)):
                shutil.rmtree(os.path.join(self.output_dir, self.project_name))
            subprocess.run(['./run_dawg.sh', self.output_dir, 
                            self.gen_output_dir, self.project_name,
                            dawg_config_file, 
                            self.intergenic_coor_f, 
                            str(cpu_num)])
            success=os.path.isfile(
                os.path.join(self.output_dir, self.project_name, 'intergenic.fa'))
            t=t+1
        if (not success):
            raise FileExistsError('Host {}: evolution of integenic regions failed'.format(
            project_name))

