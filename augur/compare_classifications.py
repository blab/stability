import sys
import numpy as np
from tabulate import tabulate


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
def add_ddG_classify1_list(classification, ddG, tt, tb, bb):
    if classification == "TT":
        tt.append(ddG)
    elif classification == "TB":
        tb.append(ddG)
    elif classification == "BB":
        bb.append(ddG)

def classify1(file, pdb_name):
    print("ddG averages for the structure, " + pdb_name)
    tt_list_total = []
    tb_list_total = []
    bb_list_total = []

    tt_list_e = []
    tb_list_e = []
    bb_list_e = []

    tt_list_ne = []
    tb_list_ne = []
    bb_list_ne = []
    for line in file:
        if len(line.split("\t")) > 1:
            split_line = line.split("\t")
            classification1 = split_line[0]
            ddG = float(split_line[2])
            transition_mutation = split_line[3]
            epitope = epitope_classify(transition_mutation)
            add_ddG_classify1_list(classification1, ddG, tt_list_total, tb_list_total, bb_list_total)
            if epitope == "True":
                add_ddG_classify1_list(classification1, ddG, tt_list_e, tb_list_e, bb_list_e)
            elif epitope == "False":
                add_ddG_classify1_list(classification1, ddG, tt_list_ne, tb_list_ne, bb_list_ne)
    matrix = np.matrix([[np.average(tt_list_total), np.average(tt_list_e), np.average(tt_list_ne)], [np.average(tb_list_total), np.average(tb_list_e), np.average(tb_list_ne)], [np.average(bb_list_total), np.average(bb_list_e), np.average(bb_list_ne)]])
    table = [["Trunk -> Trunk", str(np.average(tt_list_total)) + " (n=" + str(len(tt_list_total)) + ")", str(np.average(tt_list_e)) + " (n=" + str(len(tt_list_e)) + ")", str(np.average(tt_list_ne)) + " (n=" + str(len(tt_list_ne)) + ")"], ["Trunk -> Side Branch", str(np.average(tb_list_total)) + " (n=" + str(len(tb_list_total)) + ")", str(np.average(tb_list_e)) + " (n=" + str(len(tb_list_e)) + ")", str(np.average(tb_list_ne)) + " (n=" + str(len(tb_list_ne)) + ")"], ["Side Branch -> Side Branch", str(np.average(bb_list_total)) + " (n=" + str(len(bb_list_total)) + ")", str(np.average(bb_list_e)) + " (n=" + str(len(bb_list_e)) + ")", str(np.average(bb_list_ne)) + " (n=" + str(len(bb_list_ne)) + ")"]]
    print(pdb_name)
    print(tabulate(table, headers=["Classification", "All Sites", "Epitope Sites", "Non-Epitope Sites"]))
    return matrix

# add each ddG value to the corresponding list of ddG values
def add_ddG_classify2_list(classification, ddG, tt, tb, bc):
    if classification == "TT":
        tt.append(ddG)
    elif classification == "BB":
        tb.append(ddG)
    elif classification == "BC":
        bc.append(ddG)

def classify2(file, pdb_name):
    print("ddG averages for the structure, " + pdb_name)
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
            add_ddG_classify2_list(classification2, ddG, tt_list_total, tb_list_total, bc_list_total)
            if epitope == "True":
                add_ddG_classify2_list(classification2, ddG, tt_list_e, tb_list_e, bc_list_e)
            elif epitope == "False":
                add_ddG_classify2_list(classification2, ddG, tt_list_ne, tb_list_ne, bc_list_ne)
    matrix = np.matrix([[np.average(tt_list_total), np.average(tt_list_e), np.average(tt_list_ne)], [np.average(tb_list_total), np.average(tb_list_e), np.average(tb_list_ne)], [np.average(bc_list_total), np.average(bc_list_e), np.average(bc_list_ne)]])
    table = [["Trunk -> Trunk", str(np.average(tt_list_total)) + " (n=" + str(len(tt_list_total)) + ")", str(np.average(tt_list_e)) + " (n=" + str(len(tt_list_e)) + ")", str(np.average(tt_list_ne)) + " (n=" + str(len(tt_list_ne)) + ")"], ["Trunk -> Side Branch", str(np.average(tb_list_total)) + " (n=" + str(len(tb_list_total)) + ")", str(np.average(tb_list_e)) + " (n=" + str(len(tb_list_e)) + ")", str(np.average(tb_list_ne)) + " (n=" + str(len(tb_list_ne)) + ")"], ["Side Branch -> Tip", str(np.average(bc_list_total)) + " (n=" + str(len(bc_list_total)) + ")", str(np.average(bc_list_e)) + " (n=" + str(len(bc_list_e)) + ")", str(np.average(bc_list_ne)) + " (n=" + str(len(bc_list_ne)) + ")"]]
    print(pdb_name)
    print(tabulate(table, headers=["Classification", "All Sites", "Epitope Sites", "Non-Epitope Sites"]))
    return matrix

def average_pdb(pdb1_matrix, pdb2_matrix, classification_system):
    print("ddG averages from both structures")
    average_matrix = (pdb1_matrix + pdb2_matrix)
    average_matrix /= 2
    if classification_system == 1:
        table = [["Trunk -> Trunk", str(average_matrix[0, 0]) + " (n=" + str(64) + ")", str(average_matrix[0, 1]) + " (n=" + str(13) + ")", str(average_matrix[0, 2]) + " (n=" + str(15) + ")"], ["Trunk -> Side Branch", str(average_matrix[1, 0]) + " (n=" + str(248) + ")", str(average_matrix[1, 1]) + " (n=" + str(34) + ")", str(average_matrix[1, 2]) + " (n=" + str(53) + ")"], ["Side Branch -> Side Branch", str(average_matrix[2, 0]) + " (n=" + str(1039) + ")", str(average_matrix[2, 1]) + " (n=" + str(130) + ")", str(average_matrix[2, 2]) + " (n=" + str(212) + ")"]]
    else:
        table = [["Trunk -> Trunk", str(average_matrix[0, 0]) + " (n=" + str(64) + ")", str(average_matrix[0, 1]) + " (n=" + str(13) + ")", str(average_matrix[0, 2]) + " (n=" + str(15) + ")"], ["Trunk -> Side Branch", str(average_matrix[1, 0]) + " (n=" + str(498) + ")", str(average_matrix[1, 1]) + " (n=" + str(71) + ")", str(average_matrix[1, 2]) + " (n=" + str(115) + ")"], ["Side Branch -> Tip", str(average_matrix[2, 0]) + " (n=" + str(789) + ")", str(average_matrix[2, 1]) + " (n=" + str(93) + ")", str(average_matrix[2, 2]) + " (n=" + str(150) + ")"]]
    print(tabulate(table, headers=["Classification", "All Sites", "Epitope Sites", "Non-Epitope Sites"]))

def main():
    transition_ddG1 = pdb_name1 + "_transition_ddG_mutations.txt"
    transition_ddG_file1 = open(transition_ddG1, 'r')
    pdb1_summary1 = classify1(transition_ddG_file1, pdb_name1)
    transition_ddG_file1.close()
    transition_ddG_file1 = open(transition_ddG1, 'r')
    pdb1_summary2 = classify2(transition_ddG_file1, pdb_name1)
    transition_ddG_file1.close()

    print("\n")
    print("\n")
    transition_ddG2 = pdb_name2 + "_transition_ddG_mutations.txt"
    transition_ddG_file2 = open(transition_ddG2, 'r')
    pdb2_summary1 = classify1(transition_ddG_file2, pdb_name2)
    transition_ddG_file2.close()
    transition_ddG_file2 = open(transition_ddG2, 'r')
    pdb2_summary2 = classify2(transition_ddG_file2, pdb_name2)
    transition_ddG_file2.close()

    print("\n")
    print("\n")
    average_pdb(pdb1_summary1, pdb2_summary1, 1)
    average_pdb(pdb1_summary2, pdb2_summary2, 2)



epitope_mask = "0000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
example_sequ = "QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFINEDFNWTGVAQDGGSYACKRGSVNSFFSRLNWLHKSEYKYPALNVTMPNNGKFDKLYIWGVHHPSTDRDQTSLYVRASGRVTVSTKRSQQTVTPNIGSRPWVRGLSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRLIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQKGNIRCNICI*"

pdb_name1 = sys.argv[1]
pdb_name2 = sys.argv[2]


main()
