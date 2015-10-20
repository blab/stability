import os, socket, sys, time
 
'''
Runs each of the split mutation files on both 1HA0 and 2YP7 structures through feed_foldx.py
''' 
 

def main(number_of_files):
    print("Running...")
    print("Using multiple_runs_" + run_number)
    os.system("python split_file.py 0_mutation_trunk.txt " + number_of_files)
    for num in range(int(number_of_files)):
        num += 1
        print("NUM: " + str(num))
        directory = os.getcwd() + "/multiple_runs_" + run_number + "/" + str(num) + "_foldx_split/"
        print("Running feed_foldx.py for 1HA0_trimer... output to " + str(num) + "_1HA0_trimer_ddG_mutations.txt")
        os.system("python " + directory + "feed_foldx.py 1HA0_trimer " + str(num) + " " + run_number)

        print("Running feed_foldx.py for 2YP7_trimer... output to " + str(num) + "_2YP7_trimer_ddG_mutations.txt")
        os.system("python " + directory + "feed_foldx.py 2YP7_trimer " + str(num) + " " + run_number)
	

 
 
if __name__ == '__main__':
    number_of_files = sys.argv[1]
    run_number = sys.argv[2]
    sys.exit(main(number_of_files))
