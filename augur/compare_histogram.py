import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from numpy import random

# read ddG and classification from both files, should match in number of inputs, only read in if not equal to 0.0
# put in ddG list for each file, make one classification list, average the ddG lists to make new list
# then loop through and put values in respective classifications
# give to plt.hist(list, bin='', normed='True')

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

def plot_histogram(list, classification_list, desired_class, desired_color):
    if desired_class == "TT":
        classification = "Trunk"
    elif desired_class == "BB":
        classification = "Side Branch"
    else:
        classification = "Tips"
    desired_list = []
    for index in range(len(classification_list)):
        if classification_list[index] == desired_class:
            desired_list.append(list[index])
    plt.hist(desired_list, bins=(np.arange(-7, 7, 0.5)), color=desired_color, normed='True')
    plt.xlabel("Average ΔΔG for " + classification)
    plt.ylabel("Probability")

def main():
    classification_list = []

    transition_ddG1_list = []
    transition_ddG1 = pdb_name1 + "_transition_ddG_mutations.txt"
    transition_ddG_file1 = open(transition_ddG1, 'r')
    for line in transition_ddG_file1:
        line_split = line.split("\t")
        if len(line_split) > 1:
            classification = line_split[1]
            transition_ddG1 = float(line_split[2])
            if transition_ddG1 != 0.0:
                classification_list.append(classification)
                transition_ddG1_list.append(transition_ddG1)

    transition_ddG2_list = []
    transition_ddG2 = pdb_name2 + "_transition_ddG_mutations.txt"
    transition_ddG_file2 = open(transition_ddG2, 'r')
    for line in transition_ddG_file2:
        line_split = line.split("\t")
        if len(line_split) > 1:
            transition_ddG2 = float(line_split[2])
            if transition_ddG2 != 0.0:
                transition_ddG2_list.append(transition_ddG2)
    blue_color = [(0.19, 0.65, 0.68)]
    green_color = [(0.46, 0.74, 0.19)]
    grey_color = [(0.79608, 0.79608, 0.79608)]
    average_list = average_lists(transition_ddG1_list, transition_ddG2_list)
    #plot_histogram(average_list, classification_list, "TT", blue_color)
    #plot_histogram(average_list, classification_list, "BB", grey_color)
    plot_histogram(average_list, classification_list, "BC", green_color)
    plt.show()

pdb_name1 = sys.argv[1]
pdb_name2 = sys.argv[2]
main()