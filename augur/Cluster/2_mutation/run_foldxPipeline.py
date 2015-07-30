#!/usr/bin/python


import os, socket, sys, time
 
 
ENV_VARS = ( "SLURM_JOB_ID", "SLURM_NODELIST", "SLURM_NODEID",
             "SLURM_NNODES", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK",
             "SLURM_CPUS_ON_NODE", "SLURM_NTASKS_PER_NODE",
             "SLURM_TASK_PID",  "SLURM_PARTITION" )
def timer(f):
    """
    Decorator to report on approximate processor time consumed
    by a function.
    """
    def wrapper(*args, **kwargs):
        start_time = time.clock()
        start_time_wall = time.time()
        r = f(*args, **kwargs)
        print 'Run-time for %s: %s seconds' % \
            ( f.func_name, str(time.clock() - start_time) )
        print 'Elapsed wall clock time for %s: %s seconds' % \
            ( f.func_name, str(time.time() - start_time_wall) )
        return r
    return wrapper
 
def env_show():
    """
    Display information about the current execution environment.
    """
    print "hostname is %s" % socket.gethostname()
    for v in ENV_VARS:
        try:
            print "%s=%s" % ( v, os.environ[v] )
        except KeyError:
            print "%s: " % v
            continue
 
@timer
def main(file_part):
    print("Running...")
    env_show()
    sys.stdout.flush()
    print("Running feed_foldx.py for 2YP7_trimer... output to 2YP7_trimer_ddG_mutations.txt")
    os.system("srun -n 1 -c 1 -t 24:00:00 ../../usr/bin/python feed_foldx.py 2YP7_trimer " + file_part)
    print("Running feed_foldx.py for 1HA0_trimer... output to 1HA0_trimer_ddG_mutations.txt")
    os.system("srun -n 1 -c 1 -t 24:00:00 ../../usr/bin/python feed_foldx.py 1HA0_trimer " + file_part)


 
 
if __name__ == '__main__':
    file_part = sys.argv[1]
    sys.exit(main(file_part))
