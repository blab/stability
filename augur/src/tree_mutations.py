import os, re, time
import dendropy
from seq_util import *
from date_util import *

class tree_mutations(object):
    def __init__(self, **kwargs):

        self.mutations_fname = "0_mutation_file.txt"
        self.mutations_file = open(self.mutations_fname, 'w')

        self.universal_attributes = ['trunk', 'aa_seq']
        self.sample_attributes = ['date', 'strain']

    # check if the current node has any children/is a tip
    def tip_attribute(self, node):
        children = 0
        for child in node.child_nodes():
            children += 1
        if children == 0:
            return True
        else:
            return False

    def get_node_info(self, node):
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
        return node_info(str(node), attr_node['strain'], attr_node['trunk'], attr_node['tip'], attr_node['date'], attr_node['aa_seq'])



    def catalog_mutations(self):
        '''
        run through and print the mutations needed to get to each node in the tree
        print to the document mutation_trunk.txt

        '''
        for parent in self.tree.postorder_node_iter():
            for node in parent.child_nodes():
                node_information = self.get_node_info(node)
                parent_node_information = self.get_node_info(parent)
                self.mutations_file.write(str(node_information) + " | " + str(parent_node_information) +"\n")


class node_info(object):

    def __init__(self, v_hash, v_strain, v_trunk, v_tip, v_date, v_seperated_seq):
        #v_hash, v_strain, v_trunk, v_tip, v_seq = virus_info.split("\t")
        self.v_hash = v_hash
        self.v_strain = v_strain
        self.v_trunk = v_trunk
        self.v_tip = v_tip
        self.v_date = v_date
        self.v_seperated_seq = v_seperated_seq
        self.v_seq = self.v_seperated_seq['SigPep'] + self.v_seperated_seq['HA1'] + self.v_seperated_seq['HA2']
        self.attributes = [str(self.v_hash), str(self.v_strain), (str(self.v_trunk)), (str(self.v_tip)), str(self.v_date), str(self.v_seq)]

    def __str__(self):
        return "\t".join(self.attributes)