import os 
import subprocess
import re 
from Bio import Phylo
import shutil
import argparse
import logging
import threading
class SimNode(threading.Thread):
    '''
    ready to work -> work -> ready to give birth -> give birth
    '''
    def __init__(self, source_project, wtrees_dir, config_template, project_name, genome_db_f, tree_f,
                 output_dir, logger, back_board, lock):
        # This function overwrites the __init__ of Thread. 
        # If not calling it, this function won't be initiated like a Thread
        # object
        super(SimNode, self).__init__()
        self.source_project= source_project
        self.wtrees_dir= wtrees_dir
        self.config_template= config_template
        self.project_name= project_name
        self.genome_db_f= genome_db_f
        self.tree_f= tree_f
        self.output_dir= output_dir
        self.logger= logger
        self.is_done= False
        self.back_board= back_board
        self.lock= lock

#    def run(self ):
#        assert self.is_ready():
#        self.exec_func()
#        assert self.is_ready_to_give_birth():
#        children_nodes=self.give_birth()
#        self.lock.acquire()
#        self.logger.info('back board locked by project {}'.format(
#            str(self.project_name)))
#        self.back_board+=children_nodes
#        self.lock.release()
#        self.logger.info('back board released from project {}'.format(
#            str(self.project_name)))

    def is_ready(self):
        try:
            assert os.path.isfile(self.config_template)
            assert os.path.isfile(self.genome_db_f)
            assert os.path.isfile(self.tree_f)
        except AssertionError as e:
            return(False)
        else:
            return(True)

    def make_alf_config(self):
    #def make_alf_config(config_template, project_name, genome_db_f, tree_f, output_dir):
        config_template=self.config_template
        project_name=self.project_name
        genome_db_f=self.genome_db_f
        tree_f=self.tree_f
        output_dir=self.output_dir
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

    #def exec_func(self):
    #def exec_func(output_dir, alf_config_file, project_name, logger):
    def run(self):
        assert self.is_ready()
        output_dir= self.output_dir
        alf_config_file= self.make_alf_config()
        project_name= self.project_name
        logger= self.logger
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
            self.lock.acquire()
            logger.error(fe)
            logger.error('Simulation for host {} failed'.format(project_name))
            self.lock.release()
            sys.exit()

        assert self.is_ready_to_give_birth()
        children_nodes=self.give_birth()
        self.lock.acquire()
        self.logger.info('back board locked by project {}'.format(
            str(self.project_name)))
        self.back_board+=children_nodes
        self.lock.release()
        self.logger.info('back board released from project {}'.format(
            str(self.project_name)))

    def is_ready_to_give_birth(self):
        output_dir= self.output_dir
        child_source_project=self.project_name
        species_map_f=os.path.join(output_dir, child_source_project,
                                   'speciesMapping.txt')
        try:
            assert os.path.isfile(species_map_f)
        except AssertionError as e:
            return(False)
        else:
            return(True)

    def give_birth(self):
        child_source_project=self.project_name
        config_template= self.config_template
        output_dir= self.output_dir
        wtrees_dir= self.wtrees_dir
        children_projects= self.show_children_projects()
        back_board=self.back_board
        lock= self.lock
        logger= self.logger

        children_nodes= []
        for p in children_projects:
            child_project_name=p
            child_tree_f=os.path.join(wtrees_dir, 
                                '{}.nwk'.format(child_project_name))
            if not os.path.isfile(self.tree_f):
                # it is a tip
                logger.info('Tip reached')
                continue
            # find the db file
            # first, know how the species names are encoded
            species_map_f=os.path.join(output_dir, child_source_project,
                                       'speciesMapping.txt')
            species= {}
            for l in open(species_map_f, 'r'):
                d= l.strip().split('\t')
                species[d[0]]= d[1]
            alf_id= species[child_project_name] 
            child_genome_db_f= os.path.join(output_dir , child_source_project,
                                      'DB','{}.db'.format(alf_id))
            child_node= SimNode(source_project=child_source_project,
                                wtrees_dir= wtrees_dir, 
                                config_template= config_template,
                                project_name=child_project_name,
                                genome_db_f=child_genome_db_f, 
                                tree_f= child_tree_f, 
                                output_dir= output_dir,
                                logger= logger, 
                                back_board= back_board, 
                               lock= lock)
            children_nodes.append(child_node)
        print('line 145')
        self.is_done= True
        return(children_nodes)


    def show_children_projects(self):
        tr= Phylo.read(self.tree_f,'newick')
        children= [str(tip.name) for tip in tr.get_terminals()]
        return(children)

