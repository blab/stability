from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mutation_stability import mutation_stability
import os



class virus_stability(object):


    def __init__(self, hash, strain, trunk, tip, date, seq, foldx_directory, outgroup_directory):
        self.hash_code = hash
        self.strain = strain
        self.trunk = trunk
        self.tip = tip
        self.date = date
        self.seq = seq
        self.foldx_directory = foldx_directory
        self.find_outgroup(outgroup_directory)


        self.mutations_from_outgroup = set()  #does not include chain info
        self.sorted_mutation_string = ""

        self.mutations_from_parent = set()
        self.sorted_parent_mutation_string = ""

        self.pdb_structures = ["1HA0", "2YP7"]
        self.formatted_mut = {}

        self.ddg_outgroup = {}
        self.ddg_parent= {}
        self.parent_strain = ""



    def __str__(self):
        # ["hash code", "strain", "trunk (T/F)", "tip (T/F)", date, "ddG to outgroup (1HA0)", "ddG to outgroup (2YP7)", "ddG to parent (1HA0)", "ddG to parent (2YP7)", "mutation from parent" "parent strain", "aa_sequence"]
        return [self.hash_code, self.strain, self.trunk, self.tip, self.date, self.ddg_outgroup["1HA0"], self.ddg_outgroup["2YP7"], self.ddg_parent["1HA0"], self.ddg_parent["2YP7"], self.mutations_from_parent, self.parent_strain, self.seq]
        #attributes = [str(self.hash_code), str(self.strain), (str(self.trunk)), (str(self.tip)), str(self.date), str(self.seq)]
        #return "\t".join(attributes)

    def find_outgroup(self, directory):
        # get outgroup amino_acid sequence
        try:
            #print("Using H3N2_outgroup.gb for the outgroup, assumed to be Beijing 1992 strain")
            temp_outgroup = SeqIO.read(directory + 'H3N2_outgroup.gb', 'genbank')
        except:
            print("could not find H3N2_outgroup.gb which contained the outgroup file")
            raise
        dna_seq = str(temp_outgroup.seq).upper()
        coding_dna = Seq(dna_seq, generic_dna)
        protein = coding_dna.translate()
        self.outgroup_seq = str(protein)

    def align_to_outgroup(self):
        '''
        :mutates: self.mutations_from_outgroup
        aligns outgroup and self.seq, returns mutations from outgroup -> self.seq. Modifies self.mutations_from_outgroup
        :return: comma-seperated string of mutations from outgroup
        '''
        #
        mutations_set = set()
        outgroup_align_seq = self.outgroup_seq[24:]
        virus_align_seq = self.seq[24:]
        for index in range(len(virus_align_seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")
            if outgroup_align_seq[index] != virus_align_seq[index]:
                mutation = outgroup_align_seq[index] + str(site) + virus_align_seq[index]
                mutations_set.add(mutation)
        self.mutations_from_outgroup = mutations_set
        self.sorted_mutation_string = " ".join(sorted(list(self.mutations_from_outgroup)))
        return ','.join(mutations_set)

    def align_outgroup_to_sequence(self):
        '''
        :mutates: self.mutations_from_outgroup
        aligns self.seq and outgroup, returns mutations from sequence -> outgroup
        :return: comma-seperated string of mutations from structure
        '''
        #
        mutations_set = set()
        outgroup_align_seq = self.outgroup_seq[24:]
        for index in range(len(self.seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")
            if outgroup_align_seq[index] != self.seq[index]:
                mutation = self.seq[index] + str(site) + outgroup_align_seq[index]
                mutations_set.add(mutation)
        self.mutations_from_outgroup = mutations_set
        return ','.join(mutations_set)

    def find_mutations(self):
        '''
        Finds and stores mutations that are valid for each of the structures specified
        '''
        for structure in self.pdb_structures:
            list_of_mutations = list(self.mutations_from_outgroup)
            mut_stability = mutation_stability(list_of_mutations, structure)
            self.formatted_mut[structure] = mut_stability.get_formatted_mutations()

    def get_parent_mutations(self, parent):

        self.mutations_from_parent = self.mutations_from_outgroup - parent.mutations_from_outgroup
        self.sorted_parent_mutation_string = " ".join(sorted(list(self.mutations_from_parent)))


    def overwrite_mutation_file(self, structure):
        '''
        creates or overwrites "mutation_runfile_list.txt" which includes foldx formated, comma-seperated list of mutations for current virus
        :param structure: specify structure to be used in order to use right list of mutations
        '''
        try:
            mutations = self.formatted_mut[structure]
        except:
            print("mutations were not formatted for this structure")
            raise

        mutationFileName = "individual_list.txt"
        mutationFile = open(mutationFileName, 'w')  # overwrites the current file for FoldX
        mutationFile.write(mutations + ";")
        mutationFile.close()

    def make_run_file(self, structure, extension):
        '''
        makes the run file needed for foldx to run.
        :param structure: specify the structure to be used by foldx
        '''
        pdb_file_name = structure + extension
        runfileName = "mutate_runfile.txt"
        runfile = open(runfileName, 'w')
        runfile.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<BuildModel>mutant1,%s;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;' % (pdb_file_name, "individual_list.txt"))


    def read_ddG_output(self, structure):
        '''
        opens the output of the mutation command in foldX and gets the ddG value for the mutation that was just performed
        :param structure: specify the structure that was used by foldx
        '''
        ddGFileName = "Average_mutant1"
        ddGFile = open(ddGFileName, 'r')
        try:
            for line in ddGFile:
                if line.startswith(structure):
                    ddGline = line.split()
                    ddG = ddGline[2]
        except:
            print("couldn't find ddG output")
            raise
        ddGFile.close()
        self.ddg_outgroup[structure] = ddG

    def calculate_ddg_outgroup(self, calculated_stabilities):
        '''
        calls appropriate functions to calculate  ddG using foldx for the specified structure and ddG gets assigned to self.ddG_outgroup
        '''
        stabilities = calculated_stabilities[self.seq]
        self.ddg_outgroup['1HA0'] = stabilities[0]
        self.ddg_outgroup['2YP7'] = stabilities[1]
        '''
        for structure in self.pdb_structures:

            self.align_to_outgroup()
            self.find_mutations()
            os.chdir(self.foldx_directory)
            self.make_run_file(structure, "_trimer_repaired_1.pdb")
            self.overwrite_mutation_file(structure)
            if os.path.exists('foldx3b6'):
                os.system("./foldx3b6 -runfile mutate_runfile.txt")
            else:
                print("could not call foldx")
                raise FileNotFoundError
            self.read_ddG_output(structure)
            os.chdir("../")
        '''

    def check_valid_structure(self, structure):
        '''
        :raises Exception if the structure given is not in the allowed structures for this pipeline
        '''
        if structure not in self.pdb_structures:
            raise Exception("This pipeline does not work for that structure, only works for 1HA0 or 2YP7")

    def calculate_ddg_parent(self, parent):
        '''
        calculate the change in stability from the parent to the current virus
        '''
        for pdb in self.pdb_structures:
            try:
                self.ddg_parent[pdb] = self.ddg_outgroup[pdb] - parent.ddg_outgroup[pdb]
            except:
                print("could not calculate ddg from parent for this pair")
                print(self.strain + " | " + parent.strain)
                raise Exception("Could not calculate ddG from parent")
        self.get_parent_mutations(parent)
        self.parent_strain = parent.strain

    def mutate_pdb_to_outgroup(self, structure):
        print("- Mutating the " + structure + " structure to the outgroup sequence")
        print("--Checking valid structure given")
        self.check_valid_structure(structure)
        print("--Aligning the structure to the outgroup")
        self.align_outgroup_to_sequence()
        print("--Checking that mutations are valid for foldx and the structure")
        self.find_mutations(structure)
        print("--Making runfile and mutation file for foldx run")
        self.make_run_file(structure, "_trimer_repaired.pdb")
        self.overwrite_mutation_file(structure)
        print("--Running foldx!")
        os.system("./foldx3b6 -runfile mutate_runfile.txt")



'''
    def check_length_sequence(self):
        # check that the length sequence is always the same. Since by the time it comes out of tree_mutations,
        # should always be the same
'''
