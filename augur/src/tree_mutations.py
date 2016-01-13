import os, re, time
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability
import boto3


class tree_mutations(object):
	'''
	Run through and add each node to dictionary from hash code to virus object, and list of current virus and parent virus pairs.
    Also determine which sequences have not yet had stability calculated for them relative to Beijing outgroup. 
    Print those new sequences to /stability-data/new_seq_file.txt
    '''
	
    def __init__(self, **kwargs):

        self.universal_attributes = ['trunk', 'aa_seq']
        self.sample_attributes = ['date', 'strain']

        self.hash_to_virus = {}  # dictionary from accession to virus object
        self.virus_and_parent = []  # set of lists; containing each virus and it's parent

        self.local_storage_name = "ddg_output_database.txt"
        self.sequences_calculated = set()
        self.new_sequences = []
        self.stability_output = "stability-data/"
        new_seq_fname = self.stability_output + "new_seq_file.txt"
        self.new_seq_file = open(new_seq_fname, 'w')

        dynamodb = boto3.resource('dynamodb')
        self.table=dynamodb.Table('stability')

    def tip_attribute(self, node):
        '''
        checks if the current node has is a tip. True if it has no children, false otherwise
        '''
        children = 0
        for child in node.child_nodes():
            children += 1
        if children == 0:
            return True
        else:
            return False

    def get_node_info(self, node):
        '''
        extracts attribute information from the node
        :returns a virus_stability object with attribute information
        '''
        attr_node = {}
        for attr in self.universal_attributes:
            try:
                attr_node[attr] = (getattr(node, attr))
            except:
                print("Node was missing universal attribute hash code or trunk")
                raise
        attr_node['tip'] = self.tip_attribute(node)
        for attr in self.sample_attributes:
            if hasattr(node, attr):
                attr_node[attr] = (getattr(node, attr))
            else:
                attr_node[attr] = None
        try:  # since aa_seq is split between SigPep, HA1 and HA2, need to combine them
            v_seperated_seq = attr_node['aa_seq']
            attr_node['aa_seq'] = v_seperated_seq['SigPep'] + v_seperated_seq['HA1'] + v_seperated_seq['HA2']
        except:
            print("Couldn't combine a samples amino acid sequence chains")
            raise
        return virus_stability(str(node), attr_node['strain'], attr_node['trunk'], attr_node['tip'], attr_node['date'], attr_node['aa_seq'], self.stability_output, "source-data/")

    def read_current_database(self):
        '''
        read the dynamodb 'stability' database
        read the current local database to see which sequences already have had their stability calculated
        '''




        if os.path.isfile(self.stability_output + self.local_storage_name):
            read_local_storage_file = open(self.stability_output + self.local_storage_name, 'r')
            for line in read_local_storage_file:
                ddg_1HA0, ddg_2YP7, sequence = line.split("\t")
                if sequence not in self.sequences_calculated:
                    self.sequences_calculated.add(sequence)
        else:
            print("No local ddG storage, creating new file")
            read_local_storage_file = open(self.stability_output + self.local_storage_name, 'w')
        read_local_storage_file.close()

    def check_dynamodb(self, sequence):
        '''
        checks the stability table to see if the sequence already has had stability calculated for it
        :return returns true if in database, false if not in database
        '''
        response = self.table.get_item(
            ProjectionExpression='ddg',
            Key={'sequence':sequence}
        )
        return 'Item' in response.keys()


    def determine_new_sequences(self):
        '''
        Determine which unique sequences that have not yet had stability calculated for them.
        Then print all those sequences to stability-data/new_seq_file.txt
        '''
        print("Determining which sequences need to have stability calculated before continuing")
        for virus in self.hash_to_virus.values():
            if virus.seq not in self.new_sequences and not self.check_dynamodb(virus.seq):
                self.new_sequences.append(virus.seq)
        for seq in self.new_sequences:
            self.new_seq_file.write(seq + "\n")
        print(str(len(self.hash_to_virus.values())-len(self.new_sequences)) + " sequences were found in the database")
        print("There were " + str(len(self.new_sequences)) + " new sequences, please calculate their stabilities on the cluster before continuing")


    def catalog_mutations(self):
        '''
        Run through and add each node to dictionary from hash code to virus object, and list of current virus and parent virus pairs.
        Also determine which sequences have not yet had stability calculated for them relative to Beijing outgroup. 
        Print those new sequences to /stability-data/new_seq_file.txt

        '''
        for parent in self.tree.postorder_node_iter():
            for node in parent.child_nodes():
                if str(node) in self.hash_to_virus.keys():  # already created virus_stability object
                    node_virus = self.hash_to_virus[str(node)]
                else:  # create new virus_stability object
                    node_virus = self.get_node_info(node)
                    self.hash_to_virus[str(node)] = node_virus
                if str(parent) in self.hash_to_virus.keys():
                    parent_virus = self.hash_to_virus[str(parent)]
                else:
                    parent_virus = self.get_node_info(parent)
                self.virus_and_parent.append([node_virus, parent_virus])
                #self.mutations_file.write(str(node_virus) + " | " + str(parent_virus) +"\n")  # won't need this once dump implemented
        self.read_current_database()
        self.determine_new_sequences()