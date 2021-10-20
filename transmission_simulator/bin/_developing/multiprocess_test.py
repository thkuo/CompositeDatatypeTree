from multiprocessing import Process, Pool, Pipe
import time

def task(x):
    time.sleep(1)
    print('input: {}'.format(str(x)))
    return(x**2)
    

#for n in range(10):
#    task(n)

#parallel version
with Pool(processes= 5) as pool:
    print(pool.map(task, range(10)))
    print([pool.apply(task,args= (x, )) for x in range(10)])
    #print(pool.apply(task,range(10)))
    print([pool.apply_async(task,args= (x, )) for x in range(10)])
    print([pool.apply_async(task,args= (x, )).get() for x in range(10)])
         


#class worker:
#    def __init__(self, materials):
#        # materials: function arguments 
#        self.materials= materials 
#    def do(self, mission):
#        return(mission(self.materials))
#
#class production_line:
#    def __init(self, spaces, materials):
#        self.spaces= spaces
#    def mission_for_every(self, task, x):
#        return(sum(x)) 
#    def achieve(self, mission):
#        mission(self.materials)
#        # materials: function arguments 
##        for m in self.materials:
##            # find a free space
##            free_space_found= False
##            free_space_ix= 0
##            while not free_space_found:
##                free_space_ix=(free_space_ix+1)%len(spaces)
##            spaces[free_space_ix]= worker(mission, m)
#
#ws= [worker([x, 0-x]) for x in range(10)]
#print([ws[x].do(sum) for x in range(10)])
#print([ws[x].do(lambda x: x[0]*x[1]) for x in range(10)])
'''
def task(self, x):
    time.sleep(x)
    print(x)
    if x<10:
        self.factory_todos= self.factory_todos+[x*2, x*3]

class worker:
    def __init__(self, mission, materials, factory_todos):
        self.busy=False 
        self.factory_todos= factory_todos

    def start(self, lock):
        lock.acquire()
        mission(materials)

class factory:
    def __init__(self, task_max_num, init_task_args_q, init_task_q):
        self.init_task_args_q= task_args_q 
        self.init_task_q= task_q
        self.task_max_num= task_max_num 
    def dynamic_parallel(self): 
        ## iteratively update the queue size
        task_args_q= self.init_task_args_q
        task_q= self.init_task_q 
        while len(task_args_q)>0:
            print(task_args_q)
            args= task_args_q.pop(0)
            print(args)
            ## outer layer: control input
            assigned= False
            free_ics= [x for x in range(task_max_num) if not task_q[x].is_alive()]
            while len(free_ics)==0:
                ## inner layer: determine the processor
                free_ics= [x for x in range(task_max_num) if not task_q[x].is_alive()]
            print(free_ics)
            task(args)
        #    task_q[free_ics.pop(0)]= Process(target= task, args= (args,))
        #    task_q[free_ics.pop(0)].start()

        for x in range(len(task_q)):
            if task_q[x].is_alive():
                task_q[x].join()


if __name__=='__main__':
    task_max_num= 2
    task_args_q= [1]
    task_q= [Process(target= task, args= ())  for x in range(task_max_num) ]
    #dynamic_parallel(task_max_num, init_task_args_q= task_args_q, 
    #init_task_q= task_q)
'''
