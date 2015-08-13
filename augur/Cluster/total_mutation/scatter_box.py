import sys
import numpy as np
from tabulate import tabulate
from numpy import random
import matplotlib.pyplot as plt



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
    table = [["Trunk -> Trunk", str(np.average(combined_matrix[0])) + " (n=" + str(len(combined_matrix[0])) + ")", str(np.average(combined_matrix[3])) + " (n=" + str(len(combined_matrix[3])) + ")", str(np.average(combined_matrix[6])) + " (n=" + str(len(combined_matrix[6])) + ")"], ["Trunk -> Side Branch", str(np.average(combined_matrix[1])) + " (n=" + str(len(combined_matrix[1])) + ")", str(np.average(combined_matrix[4])) + " (n=" + str(len(combined_matrix[4])) + ")", str(np.average(combined_matrix[7])) + " (n=" + str(len(combined_matrix[7])) + ")"], ["Side Branch -> Tip", str(np.average(combined_matrix[2])) + " (n=" + str(len(combined_matrix[2])) + ")", str(np.average(combined_matrix[5])) + " (n=" + str(len(combined_matrix[5])) + ")", str(np.average(combined_matrix[8])) + " (n=" + str(len(combined_matrix[8])) + ")"]]
    print(tabulate(table, headers=["Classification", "All Sites", "Epitope Sites", "Non-Epitope Sites"]))
    return combined_matrix

def get_random_list(length, offset):
    random_list = []
    for number in range(length):
        random_number = random.uniform(-0.15, 0.15) + offset
        random_list.append(random_number)
    return random_list

def make_plot(matrix):
    blue_color = (0.19, 0.65, 0.68)
    green_color = (0.46, 0.74, 0.19)
    grey_color = (0.79608, 0.79608, 0.79608)
    red_color = (0.58039, 0.06667, 0.00000)
    plt.figure(figsize=(12, 10))
    #plt.ylabel("Stability Change (-ddG)")
    plt.ylim(-5, 5)
    #plt.xlim(-1, 11)
    bp = plt.boxplot(matrix, sym='', whis=0, patch_artist=True)
    for column in range(len(matrix)):
        plt.text(column + 1.3, np.average(matrix[column]) - 0.05, str(np.around(np.average(matrix[column]), 4)))
        plt.plot([column + 0.75, column + 1.25], [np.average(matrix[column]), np.average(matrix[column])], color='red')
        if column % 3 == 0:
            plt.scatter(get_random_list(len(matrix[column]), column + 1), matrix[column], color=blue_color, alpha=0.5)
        elif column %3 == 1:
            plt.scatter(get_random_list(len(matrix[column]), column + 1), matrix[column], color=grey_color, alpha=0.5)
        else:
            plt.scatter(get_random_list(len(matrix[column]), column + 1), matrix[column], color=green_color, alpha=0.5)
    for box in bp['boxes']:
        box.set(facecolor=red_color, alpha=0.25)
    for cap in bp['caps']:
        cap.set(color='grey')
    for whisker in bp['whiskers']:
        whisker.set(color='grey')
    for median in bp['medians']:
        median.set(color='purple', linewidth=2)
    #plt.xticks([1,2,3,4,5,6,7,8,9], ["Trunk at\nall sites", "Side Branch at\nall sites", "Tips at\nall sites", "Trunk at\nepitope sites", "Side Branch at\nepitope sites", "Tips at\nepitope sites", "Trunk at non-\nepitope sites", "Side Branch at non-\nepitope sites", "Tips at non-\nepitope sites"])
    plt.xticks([1,2,3,4,5,6,7,8,9], ["Trunk", "Side Branch", "Tip", "Trunk", "Side Branch", "Tip", "Trunk", "Side Branch", "Tip"])
    plt.show()


def main():
    transition_ddG1 = pdb_name1 + "_transition_ddG_mutations.txt"
    transition_ddG_file1 = open(transition_ddG1, 'r')
    pdb1_matrix = classify(transition_ddG_file1, pdb_name1)
    transition_ddG_file1.close()

    combined_matrix = average_pdb(pdb1_matrix, pdb1_matrix)
    make_plot(combined_matrix)


epitope_mask = "0000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"

pdb_name1 = sys.argv[1]


main()