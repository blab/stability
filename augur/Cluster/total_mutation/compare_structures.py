import sys
import matplotlib.pyplot as plot

pdb_name1 = sys.argv[1]
pdb_name2 = sys.argv[2]


def get_ddG_list(transition_ddG_file):
    ddG_list = []
    for line in transition_ddG_file:
        print(line)
        split_line = line.split("\t")
        if len(split_line) > 1:
            ddG = split_line[1]
            ddG_list.append(ddG)
    return ddG_list

def main():
    transition_ddG1 = pdb_name1 + "_transition_ddG_mutations.txt"
    transition_ddG_file1 = open(transition_ddG1, 'r')
    transition_ddG2 = pdb_name2 + "_transition_ddG_mutations.txt"
    transition_ddG_file2 = open(transition_ddG2, 'r')
    ddG_list1 = get_ddG_list(transition_ddG_file1)
    ddG_list2 = get_ddG_list(transition_ddG_file2)
    print(ddG_list1)
    print(ddG_list2)
    plot.plot([-10, 35], [-10, 35], color='red')
    plot.scatter(ddG_list1, ddG_list2, facecolors='none', color='blue')
    plot.show()

main()
# read in the files line by line, split, get ddG, put in list, then plot one list on x axis, other on y-axis
