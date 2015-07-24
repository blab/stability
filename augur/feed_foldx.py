'''
script to feed PDB and mutation information to FoldX in order to determine delta delta G
for trunk vs non-trunk mutations on the tree generated from next flu
'''

import os
import sys

# assumes that the PDB file has already been repaired
print("Using FoldX to estimate ddG values for trunk vs non-trunk mutations from the tree generated by nextflu")

# list of mutation and trunk information from the tree
mutation_trunk_fileName = "mutation_trunk.txt"
mutation_trunk_file = open(mutation_trunk_fileName, 'r')
# File where mutation information for foldx will be stored

# final output file where mutation, trunk and ddG information will be stored
finalFileName = "ddG_mutations.txt"
finalFile = open(finalFileName, 'w')
finalFile.close()  # in order to create a new blank document
# runfile
runfileName = "mutate_runfile.txt"
runfile = open(runfileName, 'w')

# makes a mutation run file for the specified pdb file.
def make_run_file(name_extension):
    pdb_file_name = pdb_name + name_extension + ".pdb"
    runfileName = "mutate_runfile.txt"
    runfile = open(runfileName, 'w')
    runfile.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<BuildModel>mutant1,%s;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;' % (pdb_file_name, "individual_list.txt"))

# creates a new file listing all mutations needed for current FoldX run, overwrites the current file
def overwrite_mutation_file(mutations):
    mutationFileName = "individual_list.txt"
    mutationFile = open(mutationFileName, 'w')  # overwrites the current file for FoldX
    mutationFile.write(mutations + ";")
    mutationFile.close()

# opens the output of the mutation command in foldX and gets the ddG value for the mutation that was just performed
def get_ddG():
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
def write_final_doc(parent_trunk, parent_ddG, parent_mutations, child_trunk, child_ddG, child_mutations):
    finalFile = open(finalFileName, 'a')
    finalFile.write(parent_trunk + "\t" + str(parent_ddG) + "\t" + parent_mutations + "\t" + child_trunk + "\t" + str(child_ddG) + "\t" + child_mutations + "\n")
    finalFile.close()

# sometimes the mutations are not ordered the same in the file, so in order to compare them in the dictionary
# we order them alphabetically in a list so that they are easily compared
def make_ordered_list(mutations):
    mutations = mutations[:len(mutations) -1]
    list_mutations = mutations.split(",")
    list_mutations = sorted(list_mutations)
    ordered_mutations = ','.join(list_mutations)
    return ordered_mutations

def run_mutations(mutations_run_dictionary, mutations):
    ddG = 0.0
    if mutations != "":
        if mutations in mutations_run_dictionary:
            ddG = mutations_run_dictionary[mutations]
        else:
            make_run_file("_formatted_1")
            overwrite_mutation_file(mutations)
            os.system("./foldx3b6 -runfile mutate_runfile.txt")
            ddG = get_ddG()
            mutations_run_dictionary[mutations] = ddG
    return ddG

def main():
    foldx_round = 0
    mutations_run_dictionary = {}
    for line in mutation_trunk_file:
        if line != "":
            split_line = line.split("\t")
            # this makes the first mutations needed to make the root == structure so that all subsequent mutations are found
            if split_line[0] == "RootMut":
                print("*** Making mutations needed to make the root equal to the protein structure we are using.")
                print(line)
                mutations = split_line[1].strip("\n")
                overwrite_mutation_file(mutations)
                make_run_file("_formatted")
                os.system("./foldx3b6 -runfile mutate_runfile.txt")
                ddG = get_ddG()
                write_final_doc("Root", ddG, mutations, "Root", ddG, mutations)
            else:
                foldx_round += 1
                print("*** This is round " + str(foldx_round) + " for foldX mutations")
                print(line)
                parent_trunk = split_line[0]
                child_trunk = split_line[2]
                parent_mutations = split_line[1].strip("\n")
                child_mutations = split_line[3].strip("\n")
                print("Running parent mutations: " + parent_trunk + " " + parent_mutations)
                print("Running child mutations: " + child_trunk + " " + child_mutations)
                parent_ddG = run_mutations(mutations_run_dictionary, parent_mutations)
                child_ddG = run_mutations(mutations_run_dictionary, child_mutations)
                print("*** Parent_DDG: " + str(parent_ddG))
                print("*** Child_DDG: " + str(parent_ddG))
                write_final_doc(parent_trunk, parent_ddG, parent_mutations, child_trunk, child_ddG, child_mutations)

pdb_name = sys.argv[1]
main()