import os, socket, sys, time
 
'''
Runs each of the split mutation files on both 1HA0 and 2YP7 structures through feed_foldx.py
''' 
 

def main(number_of_files):
    print("Running...")
    os.system("python split_file.py 0_mutation_trunk.txt " + number_of_files)
    for num in range(int(number_of_files)):
        print("NUM: " + num)
        print("Running feed_foldx.py for 1HA0_trimer... output to " + str(num) + "_1HA0_trimer_ddG_mutations.txt")
        os.system("python feed_foldx.py 1HA0_trimer " + str(num))

        print("Running feed_foldx.py for 2YP7_trimer... output to " + str(num) + "_2YP7_trimer_ddG_mutations.txt")
        os.system("python feed_foldx.py 2YP7_trimer " + str(num))
	

 
 
if __name__ == '__main__':
    number_of_files = sys.argv[1]
    sys.exit(main(number_of_files))
