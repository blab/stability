import os, re, time
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability

class tree_stability(object):

    output_file_name = "0_ddG_mutations.txt"
    try:
        output_file = open(output_file_name, 'w')
    except:
        print("can't create output file in current directory")
        raise



    def __init__(self, **kwargs):
        self.pdb_structures = ["1HA0", "2YP7"] # can add functionality to parser later
        fout_headings_begin = ["hash_code", "accession", "trunk_(T/F)", "tip_(T/F)"]
        fout_headings_end = ["parent accession", "aa_sequence"]
        fout_ddg_headings = []
        for pdb in self.pdb_structures:
            fout_ddg_headings.append("ddg_to_outgroup_" + pdb)
            fout_ddg_headings.append("ddg_to_parent_" + pdb)
        self.fout_headings = fout_headings_begin + fout_ddg_headings + fout_headings_end

        self.local_storage_name = "ddg_output_database"
        self.ddg_calculated = {}


    def calculate_stability(self):
        '''
        Calculate ddG between virus and outgroup and virus and parent
        '''
        self.read_current_database()
        self.viruses_outgroup_ddG()
        self.viruses_parent_ddG()

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.hash_to_virus.values():
            virus.find_mutations()
            virus.calculate_ddg_outgroup()
            self.store_calculation(virus.sorted_mutation_string, virus.ddg_outgroup['1HA0'], virus.ddg_outgroup['2YP7'])

    def viruses_parent_ddG(self):
        '''
        go through each virus object and calculate ddG from the parent to the current virus
        '''
        for pair in self.virus_and_parent:
            virus = pair[0]
            parent = pair[1]
            virus.calculate_ddg_parent(parent)
            virus.get_parent_mutations(parent)

    def read_current_database(self):
        read_local_storage_file = open(self.local_storage_name, 'r')
        for line in read_local_storage_file:
            mutations, ddg_1HA0, ddg_2YP7 = line.split("\t")
            self.ddg_calculated[mutations] = [ddg_1HA0, ddg_2YP7]
        read_local_storage_file.close()


    def store_calculation(self, sorted_mutation, ddg_1HA0, ddg_2YP7):
        local_storage_name = "ddg_output_database"
        local_storage_file = open(local_storage_name, 'a')
        if sorted_mutation not in self.ddg_calculated.keys():
            local_storage_file.write("\t".join([sorted_mutation, ddg_1HA0, ddg_2YP7]))