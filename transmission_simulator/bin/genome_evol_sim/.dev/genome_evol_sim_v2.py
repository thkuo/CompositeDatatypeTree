'''
Compared to v1, this version aims to avoid back refereing in each child
recursion
'''

import os
import subprocess
import re
from Bio import Phylo
import shutil
import argparse
import logging
import time
from GenomeEvolutionSimUtils.SimUtils import *
from threading import Lock
import threading

def void(logger, thread_ix):
    logger.info('Initiate thread {}'.format(thread_ix))


def make_logger(log_f= ''):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        '%(filename)s %(asctime)s %(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    if log_f == '':
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    else:
        fh= logging.FileHandler(log_f)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(ch)
    return(logger)


def main(args):
    output_dir=args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    logger= make_logger(args.log_f)

    ## generation 0
    genome_db_f = args.genome_db
    wtrees_dir = args.wtrees_dir
    try:
        assert os.path.isfile(genome_db_f)
        assert os.path.isdir(wtrees_dir)
    except AssertionError as e:
        logger.error(e)
        sys.exit()

    # determine the starting index
    starting_ix= max([int(re.sub('.nwk$', '', f))
                      for f in os.listdir(wtrees_dir) if
                      re.search('.nwk$', f)])
    project_name= str(starting_ix)
    tree_f=(os.path.join(wtrees_dir, '{}.nwk'.format(project_name)))
    config_template=args.conf
    try:
        assert os.path.isfile(config_template)
    except AssertionError as e:
        logger.error(e)
        sys.exit()

    # set up multithreading
    cpu_num= int(args.cpu_num)
    threads= [threading.Thread(target= void, args= (logger, n)) for n in range(cpu_num)]
    for n in range(len(threads)):
        threads[n].start()

    # shared by all tasks:
    # this is to avoid racing over the back_board (threads)
    lock= Lock()
    # the back board of every thread
    all_evol_queue= []
    root_task= SimNode(source_project='',
                       wtrees_dir= wtrees_dir,
                       config_template= config_template,
                       project_name= project_name,
                       genome_db_f= genome_db_f,
                       tree_f= tree_f,
                       output_dir= output_dir,
                       logger= logger,
#                       back_board= all_evol_queue,
                       lock= lock)
    all_evol_queue.append(root_task)
    threads[0]= threading.Thread(target= exec_func,
                                 args=(all_evol_queue.pop(),),
                                 daemon= True)
    threads[0].start()
    logger.info(('New task loaded; {} cpus occupied; '
                 '{} tasks in the queue').format(
                     str(sum([t.is_alive() for t in threads])),
                     str(len(all_evol_queue)))
               )
    thread_ix= 0
    # end when there is no more task and every task is done
    while (len(all_evol_queue)>0 or
           any([t.is_alive() for t in threads])):
        # when there is no new task
        if len(all_evol_queue)==0 :
            continue
        # when there is any new task
        while threads[thread_ix].is_alive():
            thread_ix= (thread_ix+1)%cpu_num

        threads[thread_ix]= threading.Thread(target= exec_func,
                                             args=(all_evol_queue.pop(),),
                                             daemon= True)
        threads[thread_ix].start()
        logger.info(('New task loaded; {} cpus occupied; '
                     '{} tasks in the queue').format(
                         str(sum([t.is_alive() for t in threads])),
                         str(len(all_evol_queue)))
                   )
    # ensure that every task is complete before the script exits
    for n in range(len(threads)):
        threads[n].join()

if __name__=='__main__':
    #+++++
    # arguments

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
    parser.add_argument('-l', dest='log_f', 
                        required= False, default= '',
                        help='filename to save log; otherwise, stdout')
    parser.add_argument('-n', dest='cpu_num', 
                        required= False, default= 1,
                        help='max number of cpus to use')

    args = parser.parse_args()
    main(args)
