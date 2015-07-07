'''
script to feed PDB and mutation information to FoldX in order to determine delta delta G
for trunk vs non-trunk mutations on the tree generated from next flu
'''

import os

# assumes that the PDB file has already been repaired
print("Using FoldX to estimate ddG values for trunk vs non-trunk mutations from the tree")

# list of mutation and trunk information from the tree
mutation_trunk_fileName = "mutation_trunk.txt"
mutation_trunk_file = open(mutation_trunk_fileName, 'r')

# File where mutation information for foldx will be stored
mutationFileName = "individual_list.txt"

# final output file where mutation, trunk and ddG information will be stored
finalFileName = "ddG_mutations.txt"


for line in mutation_trunk_file:
	print(line)
	if line[:4] == "True" or line[:4] == "False":
		trunk = line
		finalFile = open(finalFileName, 'a')
		finalFile.write(trunk)
		finalFile.close()
	else:	
		mutation = line
		mutationFile = open(mutationFileName, 'w')  # overwrites the current file for FoldX
		mutationFile.write(mutation)
		mutationFile.close()
		# runs foldx to get ddG value which it outputs to a seperate file.
		os.system("./foldx3b6 -runfile mutate_runfile.txt")
		# need to obtain the resulting ddG value from the output file
		ddGFileName = "Average_mutant1"
		ddGFile = open(ddGFileName, 'r')
		for line in ddGFile:
			if line[:5] == "1MBN_":
				print(line)
				ddGline = line.split()
				print(ddGline)
				ddG = ddGline[2]
				print(ddG)
				# write the final information to a seperate file.
				finalFile = open(finalFileName, 'a')
				finalFile.write(mutation)
				finalFile.write(ddG)
				finalFile.close()