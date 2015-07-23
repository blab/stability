
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# makes a run file for the specified pdb file to calculate stability.
def make_run_file(name_extension):
    pdb_file_name = pdb_name + name_extension + ".pdb"
    runfileName = "stability_runfile.txt"
    runfile = open(runfileName, 'w')
    runfile.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<Stability>%s;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;' % (pdb_file_name, "Stability.txt"))

# opens the output of the stability command in foldX and gets the dG value for the current structure
def get_stability():
    stability_file_name = "Stability.txt"
    stability_file = open(stability_file_name, 'r')
    for line in stability_file:
        if line.startswith(pdb_name):
            stability_line = line.split()
            stability = stability_line[1]
    stability_file.close()
    return float(stability)

# plot the main trunk pathway with annotated mutations
def plot_trunk_stability(trunk_location, tt_line, mutations):
    plt.plot(trunk_location, tt_line, marker='o', color='blue')
    plt.annotate(mutations, (trunk_location[1], tt_line[1]), color='pink', size='xx-small')

# plot the side branches off the main trunk pathway with annotated mutations
def plot_branch_stability(trunk_location, tb_line, mutations):
    plt.plot([trunk_location, trunk_location + 1], tb_line, color='red', marker='v')
    plt.annotate(mutations, (trunk_location + 1, tb_line[1]), color='green', size='xx-small')

# classify the transition between each local parent and child node
def classify(parent_classify, child_classify):
    transition_classification = ""
    if parent_classify == "True":
        if child_classify == "True":
            transition_classification = "TT"
        else:
            transition_classification = "TB"
    elif child_classify == "False":
        transition_classification = "BB"
    else:
        transition_classification = "Root"
    return transition_classification

# add each ddG value to the corresponding list of ddG values
def add_ddG_classify_list(classification, ddG):
    if classification == "TT":
        TT_ddG_list.append(ddG)
    elif classification == "TB":
        TB_ddG_list.append(ddG)
    elif classification == "BB":
        BB_ddG_list.append(ddG)

# print the averages of the ddG lists for each classification
def print_summary_ddG():
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

#
#def find_unique_mutations():

# tracks the main trunk stability values, used to plot main trunk and sidebranches
def update_stability_pathway(transition_ddG, stability_pathway):
    most_recent_dG = float(stability_pathway[len(stability_pathway) - 1])
    new_stability = most_recent_dG + transition_ddG
    stability_pathway.append(new_stability)

# goes through each line of parent and child mutations and ddG to just give the corresponding transition ddGs and mutations
# organizes for plotting
def get_translation_info(line, stability_pathway, trunk_identifier):
    split_line = line.strip("\n").split("\t")
    parent_classification = split_line[0]
    child_classification = split_line[3]
    transition_classification = classify(parent_classification, child_classification)
    if transition_classification != "Root":
        parent_ddG = float(split_line[1])
        child_ddG = float(split_line[4])
        transition_ddG = child_ddG - parent_ddG
        add_ddG_classify_list(transition_classification, transition_ddG)

        parent_mutations = split_line[2]
        parent_mutations_set = set(parent_mutations.split(","))
        child_mutations = split_line[5]
        child_mutations_set = set(child_mutations.split(","))
        transition_mutations_set = child_mutations_set - parent_mutations_set

        if transition_classification == "TT":
            update_stability_pathway(transition_ddG, stability_pathway)
            index = [len(stability_pathway) - 2, len(stability_pathway) - 1]
            tt_line = [stability_pathway[index[0]], stability_pathway[index[1]]]
            plot_trunk_stability(index, tt_line, ','.join(transition_mutations_set))
            trunk_identifier.append(parent_mutations)

        elif transition_classification == "TB":
            # looks for the side branches that branch off a main trunk node
            if parent_mutations in trunk_identifier:
                trunk_location = trunk_identifier.index(parent_mutations)
                tb_line = [stability_pathway[trunk_location], stability_pathway[trunk_location] + transition_ddG]
                plot_branch_stability(trunk_location, tb_line, ','.join(transition_mutations_set))


        transition_ddG_file.write(transition_classification + "\t" + str(transition_ddG) + "\t" + ','.join(transition_mutations_set) + "\t" + parent_mutations + "\n")

def main():
    make_run_file("_formatted_1")
    os.system("./foldx3b6 -runfile stability_runfile.txt")
    stability = get_stability()
    stability_pathway = [stability]

    ddG_mutation = "ddG_mutations.txt"
    ddG_mutation_file = open(ddG_mutation, 'r')
    trunk_identifier = []
    for line in ddG_mutation_file:
        get_translation_info(line, stability_pathway, trunk_identifier)
    print_summary_ddG()
    ddG_mutation_file.close()
    transition_ddG_file.close()
    plt.show()

pdb_name = sys.argv[1]
transition_ddG = pdb_name + "_transition_ddG_mutations.txt"
transition_ddG_file = open(transition_ddG, 'w')
TT_ddG_list = []
TB_ddG_list = []
BB_ddG_list = []
main()