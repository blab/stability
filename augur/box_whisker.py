import sys
import os
import matplotlib.pyplot as plt
import numpy as np


def ddG_to_list(tt_ddG, tb_ddG, bb_ddG, classification, transition_ddG):
    if classification == "TT":
        tt_ddG.append(transition_ddG)
    elif classification == "TB":
        tb_ddG.append(transition_ddG)
    else:
        bb_ddG.append(transition_ddG)

def plot_box(total_ddG):
    plt.boxplot(total_ddG)
    #plt.show()

def plot_scatter(total_ddG):
    plt.scatter([1]*len(total_ddG[0]), total_ddG[0])
    plt.scatter([2]*len(total_ddG[1]), total_ddG[1])
    plt.scatter([3]*len(total_ddG[2]), total_ddG[2])
    plt.show()

def main():
    transition_ddG = pdb_name + "_transition_ddG_mutations.txt"
    transition_ddG_file = open(transition_ddG, 'r')
    tt_ddG = []
    tb_ddG = []
    bb_ddG = []
    for line in transition_ddG_file:
        line_split = line.split("\t")
        if len(line_split) > 1:
            classification = line_split[0]
            transition_ddG = float(line_split[1])
            ddG_to_list(tt_ddG, tb_ddG, bb_ddG, classification, transition_ddG)
    total_ddG = [tt_ddG, tb_ddG, bb_ddG]
    plot_box(total_ddG)
    plot_scatter(total_ddG)

pdb_name = sys.argv[1]
main()