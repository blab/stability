
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# makes a run file for the specified pdb file to calculate stability.
def make_run_file(name_extension):
    pdb_file_name = pdb_name + name_extension + ".pdb"
    runfileName = "stability_runfile.txt"
    runfile = open(runfileName, 'w')
    runfile.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<Stability>%s;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;' % (pdb_file_name, "stability.txt"))

# opens the output of the stability command in foldX and gets the dG value for the current structure
def get_stability():
    stability_file_name = "stability.txt"
    stability_file = open(stability_file_name, 'r')
    for line in stability_file:
        if line[:5] == pdb_name + "_":
            stability_line = line.split()
            stability = stability_line[1]
    stability_file.close()
    return float(stability)

def get_trunk_ddG(stability):
    ddG_mutation = "ddG_mutations.txt"
    ddG_mutation_file = open(ddG_mutation, 'r')
    trunk_ddG = "trunk_ddG.txt"
    trunk_ddG_file = open(trunk_ddG, 'w')
    trunk_ddG_file.write(str(stability) + "\n")
    list_stability = [stability]
    trunk_ddG_list = []
    nontrunk_ddG_list = []
    for line in ddG_mutation_file:
        trunk_line = line.split()
        if trunk_line[0] == "True":
            trunk_ddG = float(trunk_line[1])
            trunk_ddG_list.append(trunk_ddG)
            mutations = trunk_line[2]
            dG = trunk_ddG + stability
            trunk_ddG_file.write(str(dG) + "\t" + mutations + "\n")
            list_stability.append(dG)
        else:
            nontrunk_ddG = float(trunk_line[1])
            nontrunk_ddG_list.append(nontrunk_ddG)
    trunk_average = np.average(trunk_ddG_list)
    nontrunk_average = np.average(nontrunk_ddG_list)
    print("trunk average ddG = " + str(trunk_average))
    print("non-trunk average ddG = " + str(nontrunk_average))
    return list_stability

def plot_stability(list_stability):
    plt.plot(list_stability, marker='o')

    plt.show()

# classify the transition between each local parent and child node
def classify(parent_classify, child_classify):
    transition_classification = ""
    if parent_classify == "True":
        if child_classify == "True":
            transition_classification = "TT"
        else:
            transition_classification = "TB"
    else:
        transition_classification = "BB"
    return transition_classification

def add_ddG_classify_list(classification, ddG):
    if classification == "TT":
        TT_ddG_list.append(ddG)
    elif classification == "TB":
        TB_ddG_list.append(ddG)
    else:
        BB_ddG_list.append(ddG)

def print_summary_ddG(file):
    TT_average = np.average(TT_ddG_list)
    TB_average = np.average(TB_ddG_list)
    BB_average = np.average(BB_ddG_list)
    print("ddG averages for different transition types in phylogenetic tree from HA sequences")
    print("Trunk -> Trunk transitions = " + str(TT_average))
    print("Trunk -> Branch transitions = " + str(TB_average))
    print("Branch -> Branch transitions = " + str(BB_average))
    transition_ddG_file.write("ddG averages for different transition types in phylogenetic tree from HA sequences" + "\n")
    transition_ddG_file.write("Trunk -> Trunk transitions = " + str(TT_average) + "\n")
    transition_ddG_file.write("Trunk -> Branch transitions = " + str(TB_average) + "\n")
    transition_ddG_file.write("Branch -> Branch transitions = " + str(BB_average) + "\n")


def main():
    #make_run_file("_formatted_1")
    #os.system("./foldx3b6 -runfile stability_runfile.txt")
    #stability = get_stability()
    #list_stability = get_trunk_ddG(stability)
    #plot_stability(list_stability)
    ddG_mutation = "ddG_mutations.txt"
    ddG_mutation_file = open(ddG_mutation, 'r')
    for line in ddG_mutation_file:
        split_line = line.strip("\n").split("\t")

        parent_classify = split_line[0]
        child_classify = split_line[3]
        transition_classification = classify(parent_classify, child_classify)

        parent_ddG = float(split_line[1])
        child_ddG = float(split_line[4])
        transition_ddG = child_ddG - parent_ddG
        add_ddG_classify_list(transition_classification, transition_ddG)

        parent_mutations = split_line[2]
        # delete these after the first run since fixed in tree_mutations
        parent_mutations = parent_mutations[:len(parent_mutations) - 1]
        parent_mutations_set = set(parent_mutations.split(","))
        child_mutations = split_line[5]
        # delete these after the first run since fixed in tree_mutations
        child_mutations = child_mutations[:len(child_mutations) - 1]
        child_mutations_set = set(child_mutations.split(","))
        transition_mutations = child_mutations_set - parent_mutations_set

        transition_ddG_file.write(transition_classification + "\t" + str(transition_ddG) + "\t" + ','.join(transition_mutations) + "\n")
    print_summary_ddG(ddG_mutation_file)
    ddG_mutation_file.close()
transition_ddG = "transition_ddG_mutations.txt"
transition_ddG_file = open(transition_ddG, 'w')
TT_ddG_list = []
TB_ddG_list = []
BB_ddG_list = []
pdb_name = sys.argv[1]
main()