import os, re, time
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability


class tree_mutations(object):
    def __init__(self, **kwargs):

        self.mutations_fname = "0_mutation_file.txt"
        self.mutations_file = open(self.mutations_fname, 'w')

        self.universal_attributes = ['trunk', 'aa_seq']
        self.sample_attributes = ['date', 'strain']

        self.hash_to_virus = {}  # dictionary from accession to virus object
        self.virus_and_parent = []  # set of lists; containing each virus and it's parent

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
        return virus_stability(str(node), attr_node['strain'], attr_node['trunk'], attr_node['tip'], attr_node['date'], attr_node['aa_seq'])



    def catalog_mutations(self):
        '''
        run through and print the mutations needed to get to each node in the tree
        print to the document mutation_trunk.txt

        '''
        for parent in self.tree.postorder_node_iter():
            for node in parent.child_nodes():
                if str(node) in self.hash_to_virus.keys():  # already created virus_stability object
                    node_virus = self.hash_to_virus[str(node)]
                else:  # create new virus_stability object
                    node_virus = self.get_node_info(node)
                if str(parent) in self.hash_to_virus.keys():
                    parent_virus = self.hash_to_virus[str(parent)]
                else:
                    parent_virus = self.get_node_info(parent)
                self.virus_and_parent.append([node_virus, parent_virus])
                self.mutations_file.write(str(node_virus) + " | " + str(parent_virus) +"\n")  # won't need this once dump implemented