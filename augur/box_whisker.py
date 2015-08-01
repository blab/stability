import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from numpy import random


def ddG_to_list(tt_ddG, tb_ddG, bb_ddG, classification, transition_ddG):
    if classification == "TT":
        tt_ddG.append(transition_ddG)
    elif classification == "TB":
        tb_ddG.append(transition_ddG)
    else:
        bb_ddG.append(transition_ddG)

def get_random_list(length, offset):
    random_list = []
    for number in range(length):
        random_number = random.uniform(-0.1, 0.1) + offset
        random_list.append(random_number)
    print(random_list)
    return random_list
def plot_box(total_ddG):
    plt.boxplot(total_ddG)
    #plt.show()

def plot_scatter(total_ddG):
    plt.ylim((-5, 10))
    plt.scatter(get_random_list(len(total_ddG[0]), 1), total_ddG[0], facecolors='none', color='blue')
    plt.scatter(get_random_list(len(total_ddG[1]), 2), total_ddG[1], facecolors='none', color='blue')
    plt.scatter(get_random_list(len(total_ddG[2]), 3), total_ddG[2], facecolors='none', color='blue')
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
            transition_ddG = float(line_split[2])
            ddG_to_list(tt_ddG, tb_ddG, bb_ddG, classification, transition_ddG)
    total_ddG = [tt_ddG, tb_ddG, bb_ddG]
    plot_box(total_ddG)
    plot_scatter(total_ddG)

pdb_name = sys.argv[1]
main()