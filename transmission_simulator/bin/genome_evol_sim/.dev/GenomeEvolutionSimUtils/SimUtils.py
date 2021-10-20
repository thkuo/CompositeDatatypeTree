import os 
import subprocess
import re 
from Bio import Phylo
import shutil
import argparse
import logging
import time

class SimNode():
    '''
    ready to work -> work -> ready to give birth -> give birth
    '''
    def __init__(self, source_project, wtrees_dir, config_template,
                 project_name, genome_db_f, tree_f,
                 output_dir, logger, lock):
#                 output_dir, logger, back_board, lock):
        self.source_project= source_project
        self.wtrees_dir= wtrees_dir
        self.config_template= config_template
        self.project_name= project_name
        self.genome_db_f= genome_db_f
        self.tree_f= tree_f
        self.output_dir= output_dir
        self.logger= logger
        self.is_done= False
#        self.back_board= back_board
        self.lock= lock

def is_ready(evolv_obj):
    files_to_test = {
        'config_template':evolv_obj.config_template,
        'db':evolv_obj.genome_db_f,
        'tree_f':evolv_obj.tree_f}
    for f_name in files_to_test:
        assert os.path.isfile(files_to_test[f_name]), (
            '{} ({}) not found'.format(f_name, files_to_test[f_name]))


def exec_func(evolv_obj):
    global all_evol_queue
    # test the materials
    is_ready(evolv_obj)
    # set the paramters
    output_dir= evolv_obj.output_dir
    alf_config_file= make_alf_config(evolv_obj)
    project_name= evolv_obj.project_name
    logger= evolv_obj.logger
    # run ALF
    try:
        t= 0
        max_t=3
        success=os.path.isfile(
                os.path.join(output_dir, project_name, 'DB', 'SE002_dna.fa'))
        while (not success) and (t<max_t):
            if os.path.isdir(os.path.join(output_dir, project_name)):
                shutil.rmtree(os.path.join(output_dir, project_name), ignore_errors=True)
            try:
                subprocess.run(['alfsim', '-o', output_dir, alf_config_file])
            except subprocess.CalledProcessError as e:
                t=t+1
        if t==max_t:
            raise RuntimeError
#            proc_out=subprocess.run(['alfsim', '-o', output_dir, alf_config_file])
#            print(proc_out)
#            success=True if proc_out.returncode == 0 else False
#            assert proc_out is subprocess.CompletedProcess, (
#                'ALF subprocess not returning CompletedProcess object')
#            success=os.path.isfile(
#                os.path.join(output_dir, project_name, 'DB', 'SE002_dna.fa'))
#            t=t+1
#        if (not success):
#            raise RuntimeError
    except RuntimeError as fe:
        logger.error(fe)
        logger.error('Simulation for host {} failed'.format(project_name))
        sys.exit()

    evolv_obj.lock.acquire()
    assert is_ready_to_give_birth(evolv_obj), (
        '{} cannot create children nodes'.format(evolv_obj.project_name))
    children_nodes=give_birth(evolv_obj)
    evolv_obj.logger.info('lock gained by project {}'.format(
        str(evolv_obj.project_name)))
    evolv_obj.logger.info('Job queue: {}'.format(
        str(len(all_evol_queue))))
    all_evol_queue+=children_nodes
#    evolv_obj.back_board+=children_nodes
    evolv_obj.logger.info('lock released from project {}'.format(
        str(evolv_obj.project_name)))
    evolv_obj.logger.info('Updated job queue: {}'.format(
        str(len(all_evol_queue))))
    evolv_obj.lock.release()

def make_alf_config(evolv_obj):
    config_template=evolv_obj.config_template
    project_name=evolv_obj.project_name
    genome_db_f=evolv_obj.genome_db_f
    tree_f=evolv_obj.tree_f
    output_dir=evolv_obj.output_dir
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


def is_ready_to_give_birth(evolv_obj):
    output_dir= evolv_obj.output_dir
    child_source_project=evolv_obj.project_name
    species_map_f=os.path.join(output_dir, child_source_project,
                               'speciesMapping.txt')
    try:
        assert os.path.isfile(species_map_f)
    except AssertionError as e:
        return(False)
    else:
        return(True)

def give_birth(evolv_obj):
    child_source_project=evolv_obj.project_name
    config_template= evolv_obj.config_template
    output_dir= evolv_obj.output_dir
    wtrees_dir= evolv_obj.wtrees_dir
    children_projects= show_children_projects(evolv_obj)
    # back_board=evolv_obj.back_board
    lock= evolv_obj.lock
    logger= evolv_obj.logger

    children_nodes= []
    for p in children_projects:
        child_project_name=p
        child_tree_f=os.path.join(wtrees_dir,
                            '{}.nwk'.format(child_project_name))
        if not os.path.isfile(child_tree_f):
            # it is a tip
            logger.info('Tip or pseudo taxa reached')
            continue
        # find the db file
        # first, know how the species names are encoded
        logger.info('Initiate child project {}'.format(child_project_name))
        species_map_f=os.path.join(output_dir, child_source_project,
                                   'speciesMapping.txt')
        species= {}
        for l in open(species_map_f, 'r'):
            d= l.strip().split('\t')
            species[d[0]]= d[1]
        alf_id= species[child_project_name]
        child_genome_db_f= os.path.join(output_dir, child_source_project,
                                  'DB', '{}.db'.format(alf_id))
        child_node= SimNode(source_project=child_source_project,
                            wtrees_dir= wtrees_dir, 
                            config_template= config_template,
                            project_name=child_project_name,
                            genome_db_f=child_genome_db_f, 
                            tree_f= child_tree_f, 
                            output_dir= output_dir,
                            logger= logger, 
#                            back_board= back_board, 
                            lock= lock)
        children_nodes.append(child_node)
    evolv_obj.is_done= True
    return(children_nodes)


def show_children_projects(evolv_obj):
    tr= Phylo.read(evolv_obj.tree_f,'newick')
    children= [str(tip.name) for tip in tr.get_terminals()]
    return(children)

