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
         

'''
determine hosts 
sort hosts
task queue: [[genome,host1], [noncoding, host1], [reads, host1], 
        [genome,host2], [noncoding, host2], [reads, host2], ... 
        [genome,hostx], [noncoding, hostx], [reads, hostx]]
'''
