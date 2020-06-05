import os
import subprocess
import re 
from Bio import Phylo
import shutil
import sys 

class simulator:
    def __init__(self, config_template, project_name, source_project, genome_db_f, tree_f,
                 wtrees_dir, output_dir):
        self.config_template= config_template
        self.project_name= str(project_name)
        self.source_project= source_project
        self.genome_db_f= genome_db_f
        self.tree_f= tree_f
        self.output_dir= output_dir
        if not os.path.isdir(wtrees_dir):
            sys.exit('{} not existing'.format(wtrees_dir))
        elif not os.path.isfile(tree_f):
            raise FileExistsError('Tip node')

    def run(self): 
        # find the db file
        # first, know how the previous ALF encode the species names
        species_map_f=os.path.join(self.output_dir, self.source_project,
                                   'speciesMapping.txt')
        species= {}
        for l in open(species_map_f, 'r'):
            d= l.strip().split('\t')
            species[d[0]]= d[1]
        alf_id= species[project_name] 
        genome_db_f= os.path.join(self.output_dir, self.source_project,
                                  'DB','{}.db'.format(self.alf_id))
        # make alf config file
        alf_config_file= self.make_alf_config()
        # run alf
        self.exec_func(alf_config_file)

    #+++++
    # Core functions
    def make_alf_config(self):
        # create the ALF config file using the template file
        config_file= os.path.join(output_dir, '{}.drw'.format(self.project_name))
        with open(self.config_template, 'r') as config_template_h:
            with open(config_file, 'w') as config_fh:
                for l in config_template_h:
                    if re.search('\$ORGANISM_DB\$', l):
                        l= re.sub('\$ORGANISM_DB\$', self.genome_db_f, l)
                    elif re.search('\$TREE_FILE\$', l):
                        l= re.sub('\$TREE_FILE\$', self.tree_f, l)
                    elif re.search('\$PROJECT_NAME\$', l):
                        l= re.sub('\$PROJECT_NAME\$', self.project_name, l)
                    config_fh.write(l)
        return(config_file)

    def exec_func(self, alf_config_file):
        # run ALF
        t= 0
        success=os.path.isfile(
                os.path.join(self.output_dir, self.project_name, 'DB', 'SE002_dna.fa'))
        while (not success) and (t<3):
            if os.path.isdir(os.path.join(self.output_dir, project_name)):
                shutil.rmtree(os.path.join(self.output_dir, project_name), ignore_errors=True)
            subprocess.run(['alfsim', '-o', self.output_dir, alf_config_file])
            success=os.path.isfile(
                os.path.join(self.output_dir, self.project_name, 'DB', 'SE002_dna.fa'))
            t=t+1
        if not success:
            raise FileExistsError('Host {}: genomic evolution unsuccessfully '
                'simulated'.format(self.project_name))

