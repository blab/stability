import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from numpy import random

# read ddG and classification from both files, should match in number of inputs, only read in if not equal to 0.0
# put in ddG list for each file, make one classification list, average the ddG lists to make new list
# then loop through and put values in respective classifications
# give to plt.hist(list, bin='', normed='True')

def epitope_classify(mutation):
    epitope = "NA"
    if mutation != "":
        if len(mutation.split(",")) == 1:
            for mut in mutation.split(","):
                site = int(mut[1:len(mut) - 1])
                classification = epitope_mask[site - 1]  # subtract one because of python 0 index
                if classification == "1":
                    epitope = "True"
            if epitope != "True":
                epitope = "False"
    return epitope

# add each ddG value to the corresponding list of ddG values
def add_ddG_classify_list(classification, ddG, tt, tb, bc):
    if classification == "TT":
        tt.append(ddG)
    elif classification == "BB":
        tb.append(ddG)
    elif classification == "BC":
        bc.append(ddG)

def classify(file, pdb_name):
    #print("ddG averages for the structure, " + pdb_name)
    tt_list_total = []
    tb_list_total = []
    bc_list_total = []

    tt_list_e = []
    tb_list_e = []
    bc_list_e = []

    tt_list_ne = []
    tb_list_ne = []
    bc_list_ne = []
    for line in file:
        if len(line.split("\t")) > 1:
            split_line = line.split("\t")
            classification2 = split_line[1]
            ddG = float(split_line[2])
            transition_mutation = split_line[3]
            epitope = epitope_classify(transition_mutation)
            if ddG != 0.0:
                add_ddG_classify_list(classification2, ddG, tt_list_total, tb_list_total, bc_list_total)
                if epitope == "True":
                    add_ddG_classify_list(classification2, ddG, tt_list_e, tb_list_e, bc_list_e)
                elif epitope == "False":
                    add_ddG_classify_list(classification2, ddG, tt_list_ne, tb_list_ne, bc_list_ne)
    matrix = [tt_list_total, tb_list_total, bc_list_total, tt_list_e, tb_list_e, bc_list_e, tt_list_ne, tb_list_ne, bc_list_ne]
    return matrix

def average_pdb(pdb1_matrix, pdb2_matrix):
    print("ddG averages from both structures")
    combined_matrix = []
    for column in range(len(pdb1_matrix)):
        column1 = pdb1_matrix[column]
        column2 = pdb2_matrix[column]
        combined_column = []
        for row in range(len(column1)):
            average_value = (column1[row] + column2[row]) / -2
            combined_column.append(average_value)
        combined_matrix.append(combined_column)
    return combined_matrix

def average_lists(list1, list2):
    average_list = []
    for index in range(len(list1)):
        average = np.average([list1[index], list2[index]])
        average_list.append(average)
    return average_list

# Blue:   R 48 G 167 B 173
# 0.19, 0.65, 0.68
# Green:  R 118 G 188 B 49
# 0.46, 0.74, 0.19
# bins=(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2.0),

def plot_histogram(list, classification, desired_color):
    plt.ylim(0, 0.7)
    plt.hist(list, bins=(np.arange(-7, 7, 0.5)), color=desired_color, normed='True')
    plt.xlabel("Stability for " + classification + " (-ddG)")
    plt.ylabel("Probability")
    plt.show()

def main():
    transition_ddG1 = pdb_name1 + "_transition_ddG_mutations.txt"
    transition_ddG_file1 = open(transition_ddG1, 'r')
    pdb1_matrix = classify(transition_ddG_file1, pdb_name1)
    transition_ddG_file1.close()


    transition_ddG2 = pdb_name2 + "_transition_ddG_mutations.txt"
    transition_ddG_file2 = open(transition_ddG2, 'r')
    pdb2_matrix = classify(transition_ddG_file2, pdb_name2)
    transition_ddG_file2.close()

    combined_matrix = average_pdb(pdb1_matrix, pdb2_matrix)

    blue_color = [(0.19, 0.65, 0.68)]
    green_color = [(0.46, 0.74, 0.19)]
    grey_color = [(0.79608, 0.79608, 0.79608)]

    plot_histogram(combined_matrix[0], "Trunk at all sites", blue_color)
    plot_histogram(combined_matrix[1], "Side Branches at all sites", grey_color)
    plot_histogram(combined_matrix[2], "Tips at all sites", green_color)
    plot_histogram(combined_matrix[3], "Trunk at epitope sites", blue_color)
    plot_histogram(combined_matrix[4], "Side Branches at epitope sites", grey_color)
    plot_histogram(combined_matrix[5], "Tips at epitope sites", green_color)
    plot_histogram(combined_matrix[6], "Trunk at non-epitope sites", blue_color)
    plot_histogram(combined_matrix[7], "Side Branches at non-epitope sites", grey_color)
    plot_histogram(combined_matrix[8], "Tips at non-epitope sites", green_color)

epitope_mask = "0000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"


pdb_name1 = sys.argv[1]
pdb_name2 = sys.argv[2]
main()