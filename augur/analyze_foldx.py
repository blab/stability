
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
def get_stability(pdb_name):
    stability_file_name = pdb_name + "_stability.txt"
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
    plt.ylabel("dG (kcal/mol)")
    plt.xlabel("number of mutations")
    if pdb_name == "2YP7_trimer":
        plt.annotate(mutations, (trunk_location[1], 7), color='green', size='xx-small', rotation='vertical', verticalalignment='center')
    elif pdb_name == "1HA0_trimer":
        plt.annotate(mutations, (trunk_location[1], 27), color='green', size='xx-small', rotation='vertical', verticalalignment='center')

# plot the side branches off the main trunk pathway with annotated mutations
def plot_branch_stability(branch_mutation_position, tb_line, mutations):
    plt.plot(branch_mutation_position, tb_line, color='red', marker='v')
    #plt.annotate(mutations, (branch_mutation_position[1], tb_line[1]), color='green', size='xx-small')

# classify the transition between each local parent and child node
# trunk->trunk, trunk->side branch, side branch -> side branch
def classify1(parent_classification, child_classification):
    transition_classification = ""
    if parent_classification == "True":
        if child_classification == "True":
            transition_classification = "TT"
        else:
            transition_classification = "TB"
    elif child_classification == "False":
        transition_classification = "BB"
    else:
        transition_classification = "Root"
    return transition_classification

# classify the transition between each local parent and child node
# trunk->trunk, trunk->branch or branch -> branch, branch -> tips
def classify2(parent_classification, child_classification, children_present):
    transition_classification = ""
    if parent_classification == "True":
        if child_classification == "True":
            transition_classification = "TT"
        else:
            transition_classification = "BB"
    elif parent_classification == "False":
        if children_present == "False":
            transition_classification = "BB"
        else:
            transition_classification = "BC"
    return transition_classification


# tracks the main trunk stability values, used to plot main trunk and sidebranches
def update_stability_pathway(transition_ddG, stability_pathway):
    most_recent_dG = float(stability_pathway[len(stability_pathway) - 1])
    new_stability = most_recent_dG + transition_ddG
    stability_pathway.append(new_stability)

# goes through each line of parent and child mutations and ddG to just give the corresponding transition ddGs and mutations
# organizes for plotting
def get_translation_info(line, stability_pathway, trunk_identifier, mutation_identifier, current_number_mutations):
    split_line = line.strip("\n").split("\t")

    if len(split_line) > 1:
        parent_classification = split_line[0]
        child_classification = split_line[3]
        child_present = split_line[6]
        transition_classification1 = classify1(parent_classification, child_classification)
        transition_classification2 = classify2(parent_classification, child_classification, child_present)
        parent_ddG = float(split_line[1])
        child_ddG = float(split_line[4])
        transition_ddG = child_ddG - parent_ddG

        parent_mutations = split_line[2]
        parent_mutations_set = set(parent_mutations.split(","))
        child_mutations = split_line[5]
        child_mutations_set = set(child_mutations.split(","))
        transition_mutations_set = child_mutations_set - parent_mutations_set
        number_transition_mutations = len(transition_mutations_set)
        new_number_mutations = current_number_mutations + number_transition_mutations
        index = [current_number_mutations, new_number_mutations]
        if number_transition_mutations != 0:
            if transition_classification1 == "TT":
                update_stability_pathway(transition_ddG, stability_pathway)
                index_stability_pathway = [len(stability_pathway) - 2, len(stability_pathway) - 1]
                tt_line = [stability_pathway[index_stability_pathway[0]], stability_pathway[index_stability_pathway[1]]]
                plot_trunk_stability(index, tt_line, ','.join(transition_mutations_set))
                trunk_identifier.append(parent_mutations)
                mutation_identifier.append(new_number_mutations)

            elif transition_classification1 == "TB":
                # looks for the side branches that branch off a main trunk node
                if parent_mutations in trunk_identifier:
                    trunk_location = trunk_identifier.index(parent_mutations)
                    branch_mutation_position = [mutation_identifier[trunk_location], mutation_identifier[trunk_location + 1]]
                    tb_line = [stability_pathway[trunk_location], stability_pathway[trunk_location] + transition_ddG]
                    plot_branch_stability(branch_mutation_position, tb_line, ','.join(transition_mutations_set))


        transition_ddG_file.write(transition_classification1 + "\t" + transition_classification2 + "\t" + str(transition_ddG) + "\t" + ','.join(transition_mutations_set) + "\t" + parent_mutations + "\n")
        return new_number_mutations

def main():
    #make_run_file("_formatted_1")
    #os.system("./foldx3b6 -runfile stability_runfile.txt")
    stability = get_stability(pdb_name)
    stability_pathway = [0]

    ddG_mutation = pdb_name + "_total_ddG_mutations.txt"
    ddG_mutation_file = open(ddG_mutation, 'r')
    trunk_identifier = []
    mutation_identifier = [0]
    current_number_mutations = 0
    for line in ddG_mutation_file:
        if len(line.split("\t")) > 1:
            current_number_mutations = get_translation_info(line, stability_pathway, trunk_identifier, mutation_identifier, current_number_mutations)
    ddG_mutation_file.close()
    transition_ddG_file.close()
    plt.show()

pdb_name = sys.argv[1]
transition_ddG = pdb_name + "_transition_ddG_mutations.txt"
transition_ddG_file = open(transition_ddG, 'w')
main()