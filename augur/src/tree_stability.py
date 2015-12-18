import os, re, time
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability

class tree_stability(object):

    output_file_name = split_number + "_" + "ddG_mutations.txt"
    try:
        output_file = open(output_file_name, 'w')
    except:
        print("can't create output file in current directory")
        raise



    def __init__(self, **kwargs):
        #self.pdb_structures = ["1HA0", "2YP7"] # can add functionality to parser later
        fout_headings_begin = ["hash_code", "accession", "trunk_(T/F)", "tip_(T/F)"]
        fout_headings_end = ["parent accession", "aa_sequence"]
        fout_ddg_headings = []
        for pdb in self.pdb_structures:
            fout_ddg_headings.append("ddg_to_outgroup_" + pdb)
            fout_ddg_headings.append("ddg_to_parent_" + pdb)
        self.fout_headings = fout_headings_begin + fout_ddg_headings + fout_headings_end

    def calculate_stability(self):
        '''
        Calculate ddG between virus and outgroup and virus and parent
        '''
        self.viruses_outgroup_ddG()
        self.viruses_parent_ddG()

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.hash_to_virus.values():
            virus.find_mutations()
            for pdb in self.pdb_structures:
                virus.calculate_ddG_outgroup(pdb)

    def viruses_parent_ddG(self):
        '''
        go through each virus object and calculate ddG from the parent to the current virus
        '''
        for pair in self.virus_and_parent:
            virus = pair[0]
            parent = pair[1]
            virus.calculate_ddG_parent(parent)

