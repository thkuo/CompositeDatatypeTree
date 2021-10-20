# contain all the required data
class SimNode():
    '''
    ready to work -> work -> ready to give birth -> give birth
    '''
    def __init__(self, source_project, wtrees_dir, config_template,
                 project_name, genome_db_f, tree_f,
                 output_dir, logger, back_board, lock):
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


