import os
import sys


def main():
    os.system("python src/H3N2_process.py -v" + viruses_per_month + " -y" + years)
    os.system("python feed_foldx.py " + pdb_name)
    os.system("python analyze_folx.py " + pdb_name)

pdb_name = sys.argv[1]
viruses_per_month = sys.argv[2]
years = sys.argv[3]
main()





# structure file should already be repaired,  need to specify name of structure
# need to specify which outgroup/structure to use
# need to specify range of sites in structure