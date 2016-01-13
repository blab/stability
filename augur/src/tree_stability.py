import os, re, time, shutil, boto3
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability

class tree_stability(object):
    '''
    Goes back through all virus objects, looks up their calculated ddg from outgroup. Virus_stability determines all
    their other meta data. Prints all this information to /stability-data/ddg_output.txt. Assigns ddg to the 'ep'
    attribute of all nodes.
    '''
    def __init__(self, **kwargs):
        self.pdb_structures = ["1HA0", "2YP7"] # can add functionality to parser later

        self.stability_output = "stability-data/"
        self.output_file_name = "ddg_output.txt"
        try:
            self.output_file = open(self.stability_output + self.output_file_name, 'w')
        except:
            print("can't create output file in current directory")
            raise

        dynamodb = boto3.resource('dynamodb')
        self.table=dynamodb.Table('stability')

    def calculate_stability(self):
        print("Reading in new calculated stabilities for sequences")
        self.viruses_outgroup_ddG()
        self.viruses_parent_ddG()
        self.print_viruses()
        self.assign_node_ddG()

    def print_viruses(self):
        for virus in self.hash_to_virus.values():
            print(virus)
            self.output_file.write(virus)
        self.output_file.close()

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.hash_to_virus.values():
            ddg_list = self.get_dynamodb_stability(virus.seq)
            virus.calculate_ddg_outgroup(ddg_list)

    def get_dynamodb_stability(self, sequence):
        '''
        checks the stability table to see if the sequence already has had stability calculated for it
        :return returns a list containing the stability output for that sequence, if it can't find the stability, raises an exception
        '''
        response = self.table.get_item(
            ProjectionExpression='ddg',
            Key={'sequence':sequence}
        )
        try:
            ddg_list = response['Item']['ddg']
            return ddg_list
        except:
            print("couldn't get ddg from the table for this sequence")
            print(sequence)
            raise

    def viruses_parent_ddG(self):
        '''
        go through each virus object and calculate ddG from the parent to the current virus, also get mutations from parent to new virus
        '''
        for pair in self.virus_and_parent:
            virus = pair[0]
            parent = pair[1]
            virus.calculate_ddg_parent(parent)
            virus.get_parent_mutations(parent)

    def assign_node_ddG(self):
        print("assigning ddg attribute to nodes")
        for node in self.tree.postorder_node_iter():
            hash = str(node)
            virus = self.hash_to_virus[hash]
            average_ddg = (virus.ddg_outgroup['1HA0'] + virus.ddg_outgroup['2YP7']) / 2
            try:
                setattr(node, 'ep', average_ddg)
            except:
                print("couldn't assign ddg attribute to current node")