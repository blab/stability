import os
import sys


def main():
    os.system("python src/H3N2_process.py -v" + viruses_per_month + " -y" + years)
    #os.system("python feed_foldx.py " + pdb_name1)
    #os.system("python feed_foldx.py " + pdb_name2)
    #os.system("python analyze_foldx.py " + pdb_name)

viruses_per_month = sys.argv[1]
years = sys.argv[2]
#pdb_name1 = sys.argv[3]
#pdb_name2 = sys.argv[4]
main()





# structure file should already be repaired,  need to specify name of structure
# need to specify which outgroup/structure to use
# need to specify range of sites in structure