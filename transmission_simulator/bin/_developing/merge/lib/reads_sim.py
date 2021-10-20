import os 
import subprocess

class simulator:
    def __init__(self, ix, genome_file,output_dir):
        self.ix= ix # strain id
        self.genome_file= genome_file
        self.output_dir= output_dir

    def run (self, cpu_num= 1):
        reads_dir= os.path.join(self.output_dir, 'fastq')
        if not os.path.exists(reads_dir):
            os.mkdir(reads_dir)
        # run PIRs
        self.exec_func(reads_dir, cpu_num)

    def exec_func(self, reads_dir, cpu_num= 1):
        subprocess.run(['pirs', 'simulate', '-t', str(cpu_num), 
                        '-l', '100', '-x', '100', '-m', '250', 
                        '-z', '-o', os.path.join(reads_dir, self.ix),
                        self.genome_file])

