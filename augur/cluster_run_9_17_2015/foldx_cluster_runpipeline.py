#!/usr/bin/python


import os, socket, sys, time
from subprocess import Popen, PIPE
 
 
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
 
# check each of the processes that are running. If any are still running, then return true
# otherwise return false
def check_process_completion(running_processes):
	for process in running_processes:
		return_code = process.poll()
		if return_code is None:  # process in not finished
			return True
	return False  # all processes are finished

@timer
def main(number_of_files):
    print("Running...")
    env_show()
    sys.stdout.flush()
    os.system("python split_file.py 0_mutation_trunk.txt " + number_of_files)
    print("Split Files")
    running_processes = []
    running = True
    bash_file = open("foldx_cluster_runcall.sh", 'w')
    bash_file.close()
    bash_file = open("foldx_cluster_runcall.sh", 'a')
    bash_file.write("#!/bin/sh" + "\n")
    for num in range(int(number_of_files)):
        num += 1
        print("NUM: " + str(num))
        directory = os.getcwd() + "/multiple_runs_" + run_number + "/" + str(num) + "_foldx_split/" 
        # add new processes to list of processes
        bash_file.write("srun -n 1 -c 1 -t 24:00:00 -o output_1HA0_" + str(num) + ".txt python " + directory + "feed_foldx.py 1HA0_trimer " + str(num) + " " + run_number + " & \n")
        bash_file.write("srun -n 1 -c 1 -t 24:00:00 -o output_2YP7_" + str(num) + ".txt python " + directory + "feed_foldx.py 1HA0_trimer " + str(num) + " " + run_number + " & \n")
        #running_processes.append(Popen(["./foldx_cluster_runcall.sh", directory, "1HA0", str(num), run_number], stdout=PIPE, stderr=PIPE))
        #running_processes.append(Popen(["./foldx_cluster_runcall.sh", directory, "2YP7", str(num), run_number], stdout=PIPE, stderr=PIPE))
        #bash_command = "./foldx_cluster_run_call.sh " + directory + " " + str(num) + " " + run_number
        #os.system(bash_command)     
        # os.system("python foldx_cluster_runcall.py " + directory + " " + str(num) + " " + run_number)
        '''
        print("Running feed_foldx.py for 1HA0_trimer... output to " + str(num) + "_1HA0_trimer_ddG_mutations.txt")
        os.system("srun -n 1 -c 1 -t 24:00:00 python " + directory + "feed_foldx.py 1HA0_trimer " + str(num) + " " + run_number)

        print("Running feed_foldx.py for 2YP7_trimer... output to " + str(num) + "_2YP7_trimer_ddG_mutations.txt")
        os.system("srun -n 1 -c 1 -t 24:00:00 python " + directory + "feed_foldx.py 2YP7_trimer " + str(num) + " " + run_number)
		'''
    '''
	# Keep checking for completion of all processes
	while running:
	 	time.sleep(5)
	 	running = check_process_completion(running_processes)
	'''
    bash_file.write("wait")
    bash_file.close()
    os.system("./foldx_cluster_runcall.sh")
	
if __name__ == '__main__':
    number_of_files = sys.argv[1]
    run_number = sys.argv[2]
    sys.exit(main(number_of_files))
