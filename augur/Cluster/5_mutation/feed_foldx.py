'''
script to feed PDB and mutation information to FoldX in order to determine delta delta G
for trunk vs non-trunk mutations on the tree generated from next flu
'''

import os
import sys


# makes a mutation run file for the specified pdb file.
def make_run_file(pdb_name, name_extension):
    pdb_file_name = pdb_name + name_extension + ".pdb"
    runfileName = "mutate_runfile.txt"
    runfile = open(runfileName, 'w')
    runfile.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<BuildModel>mutant1,%s;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;' % (pdb_file_name, "individual_list.txt"))

# Checks list of mutations so that only one occurs at each site.
# So if in list had S9A, A9K. Would only include S9K in list.
def check_multiple_mutations(mutations):
    mutation_dictionary = {}
    complete_mutations = ""
    mutation_list = mutations.split(",")
    if len(mutation_list) > 1:
        for mut in mutation_list:
            if len(mut) > 0:  # to get rid of that last , ""
                site = int(mut[1:len(mut) - 1])
                if site not in mutation_dictionary:
                    mutation_dictionary[site] = mut[0] + mut[len(mut) - 1]
                else:
                    mutation_dictionary[site] = mutation_dictionary.get(site)[0] + mut[len(mut) - 1]
    for site, mutation in mutation_dictionary.items():
        complete_mutations += mutation[0] + str(site) + mutation[1] + ","
    return complete_mutations[:len(complete_mutations) - 1]  # to get rid of the comma
    #print(mutations)

# Need to limit mutation sites to range of protein structure, so for 4WE4 to sites 9-501.
def check_siterange(mut):
    lowerRange = 9
    upperRange = 502
    missing_lower = 328
    missing_upper = 333
    site = int(mut[1:len(mut) - 1])
    if missing_lower <= site <= missing_upper:
        return ""
    elif lowerRange <= site <= upperRange:
        return mut + ","
    else:
        return ""


# creates a new file listing all mutations needed for current FoldX run, overwrites the current file
def overwrite_mutation_file(mutations):
    if mutations != "":
        mutations_list = mutations.split(",")
        current_run_mutations = ""
        for mut in mutations_list:
                current_run_mutations += check_siterange(mut)  # will have extra comma at end
        current_run_mutations = check_multiple_mutations(current_run_mutations)  # does not have extra comma
        mutationFileName = "individual_list.txt"
        mutationFile = open(mutationFileName, 'w')  # overwrites the current file for FoldX
        mutationFile.write(current_run_mutations + ";")
        mutationFile.close()

# opens the output of the mutation command in foldX and gets the ddG value for the mutation that was just performed
def get_ddG(pdb_name):
    ddGFileName = "Average_mutant1"
    ddGFile = open(ddGFileName, 'r')
    for line in ddGFile:
        if line.startswith(pdb_name):
            ddGline = line.split()
            ddG = ddGline[2]
    ddGFile.close()
    return float(ddG)
    # probably want to delete this file so that if it's not found can then write some error

# writes the resulting trunk, ddG and mutation information to a file called "ddG_mutations.txt"
def write_final_doc(parent_trunk, parent_ddG, parent_mutations, child_trunk, child_ddG, child_mutations, child_tip):
    finalFile = open(finalFileName, 'a')
    finalFile.write(parent_trunk + "\t" + str(parent_ddG) + "\t" + parent_mutations + "\t" + child_trunk + "\t" + str(child_ddG) + "\t" + child_mutations + "\t" + child_tip + "\n")
    finalFile.close()

# sometimes the mutations are not ordered the same in the file, so in order to compare them in the dictionary
# we order them alphabetically in a list so that they are easily compared
def make_ordered_list(mutations):
    mutations = mutations[:len(mutations) -1]
    list_mutations = mutations.split(",")
    list_mutations = sorted(list_mutations)
    ordered_mutations = ','.join(list_mutations)
    return ordered_mutations

def run_mutations(pdb_name, mutations_run_dictionary, mutations):
    ddG = 0.0
    if mutations != "":
        if mutations in mutations_run_dictionary:
            ddG = mutations_run_dictionary[mutations]
        else:
            make_run_file(pdb_name, "_formatted_1")
            overwrite_mutation_file(mutations)
            os.system("./foldx3b6 -runfile mutate_runfile.txt")
            ddG = get_ddG(pdb_name)
            mutations_run_dictionary[mutations] = ddG
    return ddG

# aligns the structure and root sequence so that the protein structure can be used for all mutations in the tree
def align(root_sequence, structure_sequence):
    print("Aligning the protein structure sequence to the root sequence")
    list_mutations = []
    for index in range(len(root_sequence)):
        site = index + 9
        if root_sequence[index] != structure_sequence[index]:
            mutation = structure_sequence[index] + str(site) + root_sequence[index]
            list_mutations.append(mutation)
    report_mutations = ','.join(list_mutations)
    print(report_mutations)
    return report_mutations

def main(pdb_name, structure_sequence):
    foldx_round = 0
    mutations_run_dictionary = {}
    for line in mutation_trunk_file:
        if line != "":
            split_line = line.split("\t")
            # this makes the first mutations needed to make the root == structure so that all subsequent mutations are found
            if split_line[0] == "Root_seq":
                print("*** Making mutations needed to make the root equal to the protein structure we are using.")
                child_root_sequence = split_line[1][8:502]
                mutations = align(child_root_sequence, structure_sequence)
                overwrite_mutation_file(mutations)
                make_run_file(pdb_name, "_formatted")
                os.system("./foldx3b6 -runfile mutate_runfile.txt")
            else:
                foldx_round += 1
                print("*** This is round " + str(foldx_round) + " for foldX mutations")
                parent_trunk = split_line[0]
                child_trunk = split_line[2]
                parent_mutations = split_line[1].strip("\n")
                child_mutations = split_line[3].strip("\n")
                print("Running parent mutations: " + parent_trunk + " " + parent_mutations)
                print("Running child mutations: " + child_trunk + " " + child_mutations)
                parent_ddG = run_mutations(pdb_name, mutations_run_dictionary, parent_mutations)
                child_ddG = run_mutations(pdb_name, mutations_run_dictionary, child_mutations)
                print("*** Parent_DDG: " + str(parent_ddG))
                print("*** Child_DDG: " + str(parent_ddG))
                child_tip = split_line[4]
                write_final_doc(parent_trunk, parent_ddG, parent_mutations, child_trunk, child_ddG, child_mutations, child_tip)



pdb_name = sys.argv[1]
part_of_file = sys.argv[2]

# assumes that the PDB file has already been repaired
print("Using FoldX to estimate ddG values for trunk vs non-trunk mutations from the tree generated by nextflu")

# list of mutation and trunk information from the tree
mutation_trunk_fileName = part_of_file + "_mutation_trunk.txt"
mutation_trunk_file = open(mutation_trunk_fileName, 'r')

if pdb_name.startswith("2YP7_trimer"):
    print("Using the 2YP7_trimer structure to calculate ddG for the given tree")
    # Runs the following protein structures with the given tree file
    structure_sequence2 = "STATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRKSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTQGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQI"
    # final output file where mutation, trunk and ddG information will be stored
    finalFileName = pdb_name + "_" + part_of_file + "_ddG_mutations.txt"
    finalFile = open(finalFileName, 'w')
    main(pdb_name, structure_sequence2)
    finalFile.close()
elif pdb_name.startswith("1HA0_trimer"):
    print("Using the 1HA0_trimer structure to calculate ddG for the given tree")
    structure_sequence1 = "STATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"
    finalFileName = pdb_name + "_" + part_of_file + "_ddG_mutations.txt"
    finalFile = open(finalFileName, 'w')
    main(pdb_name, structure_sequence1)
    finalFile.close()





