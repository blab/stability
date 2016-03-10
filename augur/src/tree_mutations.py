import os, re, time
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability
import rethinkdb as r
import hashlib


'''
    Run through and add each node to dictionary from hash code to virus object, and list of current virus and parent virus pairs.
    Also determine which sequences have not yet had stability calculated for them relative to Beijing outgroup.
    Print those new sequences to /stability-data/new_seq_file.txt
'''
class tree_mutations(object):
    def __init__(self, **kwargs):
        self.universal_attributes = ['trunk', 'aa_seq']
        self.sample_attributes = ['date', 'strain']

        self.hash_to_virus = {}  # dictionary from accession to virus object
        self.virus_and_parent = []  # set of lists; containing each virus and it's parent

        self.structures = ['2YP2', '2YP7']

        self.sequences_calculated = set()
        self.new_sequences = {}
        self.stability_output = "stability-data/"
        new_seq_fname = self.stability_output + "new_seq_file.txt"
        self.new_seq_file = open(new_seq_fname, 'w')
        self.outgroup_strain = self.outgroup['strain']

        if 'RETHINK_AUTH_KEY' in os.environ:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")
        self.database='test'
        self.table='stability'
        self.connect_rethink()

    def connect_rethink(self):
        '''
        Connect to rethink database,
        Check for existing table, otherwise create it
        '''
        try:
            r.connect(host="ec2-52-90-204-136.compute-1.amazonaws.com", port=28015, db=self.database, auth_key=self.auth_key).repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

        existing_tables = r.db(self.database).table_list().run()
        if self.table not in existing_tables:
            #create virus table with primary key 'strain'
            existing_tables = r.db(self.database).table_list().run()
            if self.table not in existing_tables:
                r.db(self.database).table_create(self.table, primary_key='md5_sequence').run()

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

    def check_rethinkdb(self, seq):
        '''
        checks the stability table to see if the sequence already has had stability calculated for it
        :return returns list of structures that have already been calculated for the sequence
        '''
        hash_function = hashlib.md5()
        hash_function.update(seq)
        hash_sequence = hash_function.hexdigest()
        document = r.table(self.table).get(hash_sequence).run()
        if document is None:
                return []
        else:
            if self.outgroup_strain in document:
                return document[self.outgroup_strain].keys()
            else:
                return[]


    def determine_new_sequences(self):
        '''
        Determine which unique sequences and structures that have not yet had stability calculated for them.
        Then print all those sequences and structuresto stability-data/new_seq_file.txt
        '''
        print("Determining which sequences need to have stability calculated before continuing")
        for virus in self.hash_to_virus.values():
            if virus.seq not in self.new_sequences:
                calculated_structures = self.check_rethinkdb(virus.seq)
                new_structures =[]
                for structure in self.structures:
                    if structure not in calculated_structures:
                        new_structures.append(structure)
                if len(new_structures) > 0:
                    self.new_sequences[virus.seq] = new_structures
                    print(self.new_sequences[virus.seq])

        for seq in self.new_sequences:
            self.new_seq_file.write(",".join(self.new_sequences[seq]) + "\t" + seq + "\n")
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
                    self.hash_to_virus[str(parent)] = parent_virus
                self.virus_and_parent.append([node_virus, parent_virus])
                #self.mutations_file.write(str(node_virus) + " | " + str(parent_virus) +"\n")  # won't need this once dump implemented
        self.determine_new_sequences()