import sys
import os
import matplotlib.pyplot as plt
import numpy as np



# plot the main trunk pathway with annotated mutations
def plot_trunk_stability(trunk_location, tt_line, mutations):
    plt.plot(trunk_location, tt_line, marker='o', color='blue', linewidth=4)
    plt.ylabel("dG (kcal/mol)")
    plt.xlabel("number of mutations")
    epitope = epitope_classify(mutations)
    if epitope == "True":
        plt.annotate(mutations, (trunk_location[1], 15), color='green', size='xx-small', rotation='vertical', verticalalignment='center')
    elif epitope == "Receptor":
        plt.annotate(mutations, (trunk_location[1], 15), color='black', size='xx-small', rotation='vertical', verticalalignment='center')
    else:
        plt.annotate(mutations, (trunk_location[1], 15), color='purple', size='xx-small', rotation='vertical', verticalalignment='center')

# plot the side branches off the main trunk pathway with annotated mutations
def plot_branch_stability(branch_mutation_position, tb_line, mutations):
    plt.plot(branch_mutation_position, tb_line, color='red', marker='v')
    #plt.annotate(mutations, (branch_mutation_position[1], tb_line[1]), color='green', size='xx-small')

def epitope_classify(mutation):
    epitope = "False"
    if mutation != "":
        for mut in mutation.split(","):
            site = int(mut[1:len(mut) - 1])
            classification = epitope_mask[site - 1]  # subtract one because of python 0 index
            if classification == "1":
                epitope = "True"
            elif site in receptor_binding_sites:
                epitope = "Receptor"
    return epitope

# tracks the main trunk stability values, used to plot main trunk and sidebranches
def update_stability_pathway(transition_ddG, stability_pathway):
    most_recent_dG = float(stability_pathway[len(stability_pathway) - 1])
    new_stability = most_recent_dG + transition_ddG
    stability_pathway.append(new_stability)

# goes through each line of parent and child mutations and ddG to just give the corresponding transition ddGs and mutations
# organizes for plotting
def get_translation_info(line1, line2, stability_pathway, trunk_identifier, mutation_identifier, current_number_mutations):
    split_line1 = line1.strip("\n").split("\t")
    split_line2 = line2.strip("\n").split("\t")
    if len(split_line1) > 1:
        transition_classification = split_line1[0]
        ddG1 = float(split_line1[2])
        ddG2 = float(split_line2[2])
        average_ddG = np.average([ddG1, ddG2])


        transition_mutations = split_line1[3]
        parent_mutations = split_line1[4]
        number_transition_mutations = len(transition_mutations.split("\t"))
        new_number_mutations = current_number_mutations + number_transition_mutations
        index = [current_number_mutations, new_number_mutations]
        if number_transition_mutations != 0:
            if transition_classification == "TT":
                update_stability_pathway(average_ddG, stability_pathway)
                index_stability_pathway = [len(stability_pathway) - 2, len(stability_pathway) - 1]
                tt_line = [stability_pathway[index_stability_pathway[0]], stability_pathway[index_stability_pathway[1]]]
                plot_trunk_stability(index, tt_line, transition_mutations)
                trunk_identifier.append(parent_mutations)
                mutation_identifier.append(new_number_mutations)

            elif transition_classification == "TB":
                # looks for the side branches that branch off a main trunk node
                if parent_mutations in trunk_identifier:
                    trunk_location = trunk_identifier.index(parent_mutations)
                    branch_mutation_position = [mutation_identifier[trunk_location], mutation_identifier[trunk_location + 1]]
                    tb_line = [stability_pathway[trunk_location], stability_pathway[trunk_location] + average_ddG]
                    plot_branch_stability(branch_mutation_position, tb_line, transition_mutations)


        return new_number_mutations

def main():
    ddG_mutation1 = pdb_name1 + "_transition_ddG_mutations.txt"
    ddG_mutation_file1 = open(ddG_mutation1, 'r')
    ddG_mutation2 = pdb_name2 + "_transition_ddG_mutations.txt"
    ddG_mutation_file2 = open(ddG_mutation2, 'r')
    trunk_identifier = []
    mutation_identifier = [0]
    current_number_mutations = 0
    stability_pathway = [0]
    for line1, line2 in zip(ddG_mutation_file1, ddG_mutation_file2):
        if len(line1.split("\t")) > 1:
            current_number_mutations = get_translation_info(line1, line2, stability_pathway, trunk_identifier, mutation_identifier, current_number_mutations)
    plt.show()

epitope_mask = "0000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
receptor_binding_sites = [145, 155, 156, 158, 159, 189, 193]
pdb_name1 = sys.argv[1]
pdb_name2 = sys.argv[2]
main()