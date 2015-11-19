import os
import sys

def main():
    print("Running feed_foldx.py for 1HA0_trimer... output to " + str(split) + "_1HA0_trimer_ddG_mutations.txt")
    os.system("srun -n 1 -c 1 -t 24:00:00 python " + directory + "feed_foldx.py 1HA0_trimer " + str(split) + " " + run_number)

    print("Running feed_foldx.py for 2YP7_trimer... output to " + str(split) + "_2YP7_trimer_ddG_mutations.txt")
    os.system("srun -n 1 -c 1 -t 24:00:00 python " + directory + "feed_foldx.py 2YP7_trimer " + str(split) + " " + run_number)	

directory = sys.argv[1]
split = sys.argv[2]
run_number = sys.argv[3]	
main()
