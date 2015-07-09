import os
import sys


def main():
    os.system("python src/H3N2_process.py -v1 -y30")
    os.system("python feed_foldx.py " + pdb_name)

pdb_name = sys.argv[1]
main()





# structure file should already be repaired, may need to specify name of structure
# need to specify which outgroup/structure to use
# need to specify range of sites in structure