from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mutation_stability import mutation_stability
import os




class virus_stability:


    def __init__(self, hash_code, accension, trunk, tip, seq):

        self.accession = accension
        self.hash_code = hash_code
        self.trunk = trunk
        self.tip = tip
        self.seq = seq
        self.mutations_from_outgroup = set()  #does not include chain info

        self.structures = ["1HA0", "2YP7"]
        self.mut1 = ""
        self.mut2 = ""

        self.ddg_outgroup1 = None
        self.ddg_outgroup2 = None
        self.ddg_parent1 = None
        self.ddg_parent2 = None

        # get outgroup amino_acid sequence
        try:
            temp_outgroup = SeqIO.read('source-data/H3N2_outgroup.gb', 'genbank')
        except:
            print("could not find source-data/H3N2_outgroup.gb which contained the outgroup file")
        dna_seq = str(temp_outgroup.seq).upper()
        coding_dna = Seq(dna_seq, generic_dna)
        protein = coding_dna.translate()
        self.outgroup_seq = str(protein)

    def __str__(self):
        return ("|hash: %s |accession: %s |trunk: %s |tip: %s\n\t|sequence: %s") %(self.hash_code, self.accession, self.trunk, self.tip, self.seq)

    def align_to_outgroup(self):
        '''
        :mutates: self.mutations_from_outgroup
        aligns outgroup and self.seq, returns mutations from outgroup -> self.seq. Modifies self.mutations_from_outgroup
        :return: comma-seperated string of mutations from outgroup
        '''
        #
        mutations_set = set()
        outgroup_align_seq = self.outgroup_seq[24:]
        for outgroup_index in range(len(self.outgroup_seq)):
            seq_index = outgroup_index
            site = outgroup_index + 9  # for both 1HA0 and 2YP7
            if self.outgroup_seq[outgroup_index] != self.seq[seq_index]:
                mutation = self.outgroup_seq[outgroup_index] + str(site) + self.seq[seq_index]
                mutations_set.add(mutation)
        self.mutations_from_outgroup = mutations_set
        return ','.join(mutations_set)

    def find_mutations(self):
        '''
        intializes self.mut1 and self.mut2 for mutations that are valid for each of the structures 1HA0 and 2YP7
        '''
        list_of_mutations = self.align_to_outgroup()
        mutations1 = mutation_stability(list_of_mutations, "1HA0")
        mutations2 = mutation_stability(list_of_mutations, "2YP7")
        self.mut1 = mutations1.get_formatted_mutations()
        self.mut2 = mutations2.get_formatted_mutations()

    def overwrite_mutation_file(self, structure):
        '''
        creates or overwrites "mutation_runfile_list.txt" which includes foldx formated, comma-seperated list of mutations for current virus
        :param structure: specify structure to be used in order to use right list of mutations
        '''

        if structure == "1HA0":
            mutations = self.mut1
        elif structure == "2YP7":
            mutations = self.mut2
        mutationFileName = "mutation_runfile_list.txt"
        mutationFile = open(mutationFileName, 'w')  # overwrites the current file for FoldX
        mutationFile.write(mutations + ";")
        mutationFile.close()

    def make_run_file(self, structure):
        '''
        makes the run file needed for foldx to run.
        :param structure: specify the structure to be used by foldx
        '''
        pdb_file_name = structure + "_trimer_repaired_1.pdb"
        runfileName = "mutate_runfile.txt"
        runfile = open(runfileName, 'w')
        runfile.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<BuildModel>mutant1,%s;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;' % (pdb_file_name, "mutation_runfile_list.txt"))


    def read_ddG_output(self, structure):
        '''
        opens the output of the mutation command in foldX and gets the ddG value for the mutation that was just performed
        :param structure: specify the structure that was used by foldx
        '''
        ddGFileName = "Average_mutant1"
        ddGFile = open(ddGFileName, 'r')
        for line in ddGFile:
            if line.startswith(structure):
                ddGline = line.split()
                ddG = ddGline[2]
        ddGFile.close()
        if structure == "1HA0":
            self.ddg_outgroup1 = ddG
        elif structure == "2YP7":
            self.ddg_outgroup2 = ddG


    def calculate_ddg_outgroup(self, structure):
        '''
        calls appropriate functions to calculate the ddG using foldx for the specified structure
        '''
        self.check_valid_structure(structure)
        self.make_run_file(structure)
        self.overwrite_mutation_file(structure)
        os.system("./foldx3b6 -runfile mutate_runfile.txt")
        self.read_ddG_output(structure)

    def check_valid_structure(self, structure):
        '''
        :raises Exception if the structure given is not in the allowed structures for this pipeline
        '''
        if structure not in self.structures:
            raise Exception("This pipeline does not work for that structure, only works for 1HA0 or 2YP7")
'''
    def check_length_sequence(self):
        # check that the length sequence is always the same. Since by the time it comes out of tree_mutations,
        # should always be the same

    def calculate_ddg_parent(self, parent):
        self.ddg_parent= self.ddg_outgroup - parent.ddg_outgroup

    def check_same_virus(self, other):
        # just compare the sequences
'''
