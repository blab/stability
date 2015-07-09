pdb_fileName = "4WE4_repaired.pdb"
pdb_file = open(pdb_fileName, 'r')

new_pdb_fileName = "4WE4_repaired_nochain.pdb"
new_pdb_file = open(new_pdb_fileName, 'w')


# this removes the chain information A or B in the 1968 4WE4 structure that was 
# giving FoldX problems in finding the mutation
for line in pdb_file:
	newLine = line
	if len(line) > 21:
		if line[21] == "A" or line[21] == "B":
			newLine = line[:21] + " " + line[22:]
	new_pdb_file.write(newLine)
new_pdb_file.close()

pdb_fileName = "4WE4_repaired_nochain.pdb"
pdb_file = open(pdb_fileName, 'r')

new_pdb_fileName = "4WE4_repaired_formatted.pdb"
new_pdb_file = open(new_pdb_fileName, 'w')

# this should renumber the second half of the pdb file since it used to be separated into 
# 2 chains. 
foundTER = False
for line in pdb_file:
	newLine = line
	if len(line) > 21 and foundTER == True:
		number = int(line[23:26])
		number += 329 
		newLine = line[:23] + str(number) + line[26:]
	if line[0:3] == "TER":
		foundTER = True
	new_pdb_file.write(newLine)
new_pdb_file.close()