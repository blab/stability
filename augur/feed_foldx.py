'''
script to feed PDB and mutation information to FoldX in order to determine delta delta G
for trunk vs non-trunk mutations on the tree generated from next flu
'''

import os

# assumes that the PDB file has already been repaired


mutation_trunk_fileName = "mutation_trunk.txt"
mutation_trunk_file = open(mutation_trunk_fileName, 'r')

count = 1
for line in mutation_trunk_file:
    trunk = line
    mutation = mutation_trunk_file.next()

    mutationRun_fileNameRead = "mutate_runfile" + str(count - 1)+ ".txt"
    mutationRun_fileNameWrite = "mutate_runfile" + str(count)+ ".txt"
    individualMutationFile = "individual_list.txt"
    repairedPDB_fileName = "1MBN_repaired.pdb"

    mutationRun_fileRead = open(mutationRun_fileNameRead, 'r')
    mutationRun_fileWrite = open(mutationRun_fileNameWrite, 'w')

    for line in mutationRun_fileRead:
        # this changes the pdb file
        if line[:5] == "<PDBS>":
            line = "<PDBS>1MBN_repaired.pdb;"
        # this changes the mutation you want to make
        elif line[:6] == "<Build":
            line = "<BuildModel>mutant1," + individualMutationFile + ";\n"
        mutationRun_fileWrite.write(line)

    mutationRun_fileWrite.close()
    os.system("./foldx3b6 -runfile mutate_runfile2.txt")
    count += 1
    print(count)