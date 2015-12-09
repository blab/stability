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
    print(mutation_list)
    if len(mutation_list) > 1:
        for mut in mutation_list:
            print(mut)
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

# Need to limit mutation sites to range of protein structure, so for 4WE4 to sites 9-501. A392S
def check_siterange(mut):
    # compiled from 2YP7 and 1HA0 to run the same residues
    lowerRange = 9
    upperRange = 502
    missing_lower = 328
    missing_upper = 333
    '''
    if pdb_name.startswith("2YP7"):
    	#2YP7
    	lowerRange = 8
    	upperRange = 503
    	missing_lower = 328
    	missing_upper = 333
    elif pdb_name.startswith("1HA0"):
    	#1HA0
    	lowerRange = 9
    	upperRange = 502
	'''
    '''
    if pdb_name.startswith("2HMG"):
        lowerRange = 1
        upperRange = 503
    elif pdb_name.startswith("4WE4"):
        lowerRange = 9
        upperRange = 501
    '''
    #print("Mut: " + mut)
    site = int(mut[1:len(mut) - 1])
    if missing_lower <= site <= missing_upper:
        return ""
    elif lowerRange <= site <= upperRange:
        #adjust1 = str(site + upperRange)
        #adjust2 = str(site + (2 * upperRange))
        return mut + ","  #+ mut[0] + adjust1 + mut[len(mut) - 1] + "," + mut[0] + adjust2 + mut[len(mut) - 1] + ","
    else:
        return ""

def include_chain_info(mutations):
    mutations_list = mutations.split(",")
    foldx_mutations = ""
    if pdb_name.startswith("1HA0"):
        chain1 = "A"
        chain2 = "M"
        chain3 = "Y"
    elif pdb_name.startswith("2YP7"):
        chain1 = "A"
        chain2 = "P"
        chain3 = "E"
    '''
    if pdb_name.startswith("4WE4"):
        start1 = 9
        end1 = 329
        start2 = 330
        end2 = 501
        chain1 = "A"
        chain2 = "U"
        chain3 = "C"
        chain4 = "W"
        chain5 = "E"
        chain6 = "Y"
    elif pdb_name.startswith("2HMG"):
        start1 = 1
        end1 = 328
        start2 = 330
        end2 = 503
        chain1 = "A"
        chain2 = "B"
        chain3 = "C"
        chain4 = "D"
        chain5 = "E"
        chain6 = "F"
    
    for mut in mutations_list:
        site = int(mut[1:len(mut) - 1])
        if start1 <= site <= end1:
            site = str(site)
            foldx_mutations += mut[0] + chain1 + site + mut[len(mut) - 1] + "," + mut[0] + chain3 + site + mut[len(mut) - 1] + "," + mut[0] + chain5 + site + mut[len(mut) - 1] + ","
        elif start2 <= site <= end2:
            if pdb_name.startswith("2HMG"):
                site -= 1  # since it was missing a residue in the middle that we subsequently added? Do we only this when aligning?
            site = str(site - end1)
            foldx_mutations += mut[0] + chain2 + site + mut[len(mut) - 1] + "," + mut[0] + chain4 + site + mut[len(mut) - 1] + "," + mut[0] + chain6 + site + mut[len(mut) - 1] + ","
    return foldx_mutations[:len(foldx_mutations) - 1]
	'''
    for mut in mutations_list:
        site = mut[1:len(mut) - 1]
        foldx_mutations += mut[0] + chain1 + site + mut[len(mut) - 1] + "," + mut[0] + chain2 + site + mut[len(mut) - 1] + "," + mut[0] + chain3 + site + mut[len(mut) - 1] + ","
    return foldx_mutations[:len(foldx_mutations) - 1]
	
# creates a new file listing all mutations needed for current FoldX run, overwrites the current file
def overwrite_mutation_file(mutations):
    if mutations != "":
        mutations_list = mutations.split(",")
        current_run_mutations = ""
        print(mutations_list)
        for mut in mutations_list:
                current_run_mutations += check_siterange(mut)  # will have extra comma at end
        current_run_mutations = check_multiple_mutations(current_run_mutations)  # does not have extra comma
        current_run_mutations = include_chain_info((current_run_mutations))  # does not have extra comma
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
def write_final_doc(parent_trunk, parent_ddG, parent_mutations, parent_hash, child_trunk, child_ddG, child_mutations, child_hash, child_tip):
    finalFile = open(finalFileName, 'a')
    finalFile.write(parent_trunk + "\t" + str(parent_ddG) + "\t" + parent_mutations + "\t" + parent_hash + "\t"+ child_trunk + "\t" + str(child_ddG) + "\t" + child_mutations + "\t" + child_hash + "\t" + child_tip + "\n")
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
            make_run_file(pdb_name, "_repaired_1")
            overwrite_mutation_file(mutations)
            os.system("./foldx4 -runfile mutate_runfile.txt")
            ddG = get_ddG(pdb_name)
            mutations_run_dictionary[mutations] = ddG
    return ddG

# aligns the structure and root sequence so that the protein structure can be used for all mutations in the tree
def align(root_sequence, structure_sequence):
    print("Aligning the protein structure sequence to the root sequence")
    list_mutations = []
    for root_index in range(len(root_sequence)):
        structure_index = root_index
        site = root_index + 9  # for both 1HA0 and 2YP7
        '''
        if pdb_name.startswith("4WE4"):
            site = root_index + 9
        elif pdb_name.startswith("2HMG"):
            site = root_index + 1
        '''
        '''
            if root_index > 328:
                structure
        '''
        #print(root_sequence[index] + " : " + structure_sequence[index])
        if root_sequence[root_index] != structure_sequence[structure_index]:
            mutation = structure_sequence[structure_index] + str(site) + root_sequence[root_index]
            list_mutations.append(mutation)
    report_mutations = ','.join(list_mutations)
    return report_mutations

def main(pdb_name, structure_sequence):
    foldx_round = 0
    mutations_run_dictionary = {}
    for line in mutation_trunk_file:
        #print(line)
        if line != "":
            split_line = line.split("\t")
            # this makes the first mutations needed to make the root == structure so that all subsequent mutations are found
            if split_line[0] == "Root_seq":
                print("*** Making mutations needed to make the root equal to the protein structure we are using.")
                '''
                if pdb_name.startswith("4WE4"):
                    child_root_sequence = split_line[1][8:501]
                elif pdb_name.startswith("2HMG"):
                    child_root_sequence = split_line[1][0:503]
                '''
                child_root_sequence = split_line[1][8:502]  #for 2YP7 and 1HA0
                child_root_hash = split_line[3]
                print(structure_sequence)
                print()
                print(child_root_sequence)
                mutations = align(child_root_sequence, structure_sequence)
                overwrite_mutation_file(mutations)
                make_run_file(pdb_name, "_repaired")
                os.system("./foldx4 -runfile mutate_runfile.txt")
                finalFile = open(finalFileName, 'a')
                finalFile.write("Root" + "\t" + child_root_hash + "\n")
                finalFile.close()
            else:
                foldx_round += 1
                print("*** This is round " + str(foldx_round) + " for foldX mutations")
                parent_trunk = split_line[0]
                child_trunk = split_line[3]
                parent_mutations = split_line[1].strip("\n")
                child_mutations = split_line[4].strip("\n")
                print("Running parent mutations: " + parent_trunk + " " + parent_mutations)
                print("Running child mutations: " + child_trunk + " " + child_mutations)
                parent_ddG = run_mutations(pdb_name, mutations_run_dictionary, parent_mutations)
                child_ddG = run_mutations(pdb_name, mutations_run_dictionary, child_mutations)
                print("*** Parent_DDG: " + str(parent_ddG))
                print("*** Child_DDG: " + str(parent_ddG))
                parent_hash = split_line[2]
                child_hash = split_line[5]
                child_tip = split_line[6]
                write_final_doc(parent_trunk, parent_ddG, parent_mutations, parent_hash, child_trunk, child_ddG, child_mutations, child_hash, child_tip)

# assumes that the PDB file has already been repaired
print("Using FoldX to estimate ddG values for trunk vs non-trunk mutations from the tree generated by nextflu")
pdb_name = sys.argv[1]
print(pdb_name)
part_of_file = sys.argv[2]
print(part_of_file)
run_number = sys.argv[3]

# list of mutation and trunk information from the tree
os.chdir("multiple_runs_" + run_number + "/" + part_of_file + "_foldx_split")
mutation_trunk_fileName = part_of_file + "_mutation_trunk.txt"
print("File Name: " + mutation_trunk_fileName)
print("CWD: " + os.getcwd())
mutation_trunk_file = open(mutation_trunk_fileName, 'r')


if pdb_name == "2YP7_trimer":
    print("Using the 2YP7_trimer structure to calculate ddG for the given tree")
    # Runs the following protein structures with the given tree file
    structure_sequence = "STATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRKSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTQGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQI"
elif pdb_name == "1HA0_trimer":
    print("Using the 1HA0_trimer structure to calculate ddG for the given tree")
    structure_sequence = "STATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"
elif pdb_name.startswith("4WE9_trimer"):
    print("Using the 4WE9_trimer structure to calculate ddG for the given tree")
    structure_sequence = "STATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSNNSFFSRLNWLTHLNFKYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQI"
elif pdb_name.startswith("2HMG_trimer"):
    print("Using the 2HMG_trimer structure to calculate ddG for the given tree")
    structure_sequence = "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSDFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKG"
elif pdb_name.startswith("4WE4_trimer"):
    print("Using the 4WE4_trimer structure to calculate ddG for the given tree")
    structure_sequence = "STATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQ"
# final output file where mutation, trunk and ddG information will be stored   

finalFileName = part_of_file + "_" + pdb_name + "_" + "ddG_mutations.txt"
finalFile = open(finalFileName, 'w')
main(pdb_name, structure_sequence)
finalFile.close()

                          #"QKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFINEDFNWTGVAQDGGSYACKRGSVNSFFSRLNWLHKSEYKYPALNVTMPNNGKFDKLYIWGVHHPSTDRDQTSLYVRASGRVTVSTKRSQQTVTPNIGSRPWVRGLSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRLIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKG"