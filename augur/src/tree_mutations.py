import os, re, time
import dendropy
from seq_util import *
from date_util import *

class tree_mutations(object):
    def __init__(self, **kwargs):
		'''
        '''

    def catalog_mutations(self):
        '''
        run through and print the mutations needed to get to each node in the tree
        print to the document mutation_trunk.txt

        '''

        # file that mutation and trunk information will be printed to
        mutation_trunk_fileName = "0_mutation_trunk.txt"
        mutation_trunk_file = open(mutation_trunk_fileName, 'w')

        # Add the current nodes mutations to the total list of mutations that were needed to get to that
        # node from the root
        def update_mutations(current_node_mutations, current_total_mutations):
            # for formatting purposes only want to add the new mutation to the list if it is not blank
            if str(current_node_mutations) != "":
                currentMutation = str(current_node_mutations).split(",")
                for mut in currentMutation:
                    current_total_mutations += mut + ","
            return current_total_mutations

        # check if the current node has any children/is a tip
        def determine_tip(node):
            tip = False
            children = 0
            for child in node.child_nodes():
                children += 1
            if children == 0:
                tip = True
            return tip

        # Go through all the nodes and print the mutations that were needed to get to that node from the root.
        # Print to "mutation_trunk.txt" trunk and mutation information.
        def node_foldx(self, node, current_total_mutations):
            if node.parent_node is None:  # root of the tree
                print("Ultimate Root : " + node.aa_seq)
                print("Ultimate Root Hash: " + str(node))
                print(len(node.child_nodes()))
                for child in node.child_nodes():
                    print("Mutations!!: " + child.aa_muts)
                    print("Hash: " + str(child))
                    print(len(child.child_nodes()))
                    if len(child.child_nodes()) != 0:
                    #if child.aa_muts != "":  # only want to print out root mutations if not the outgroup
                        print("Root Child: " + child.aa_seq)
                        print("Root Child: " + child.aa_muts)
                        print("Root Child Hash: " + str(child))
                        root_child_seq = child.aa_seq
                        root_child_hash = str(child)
                        current_total_mutations = ""
                        #current_total_mutations = update_mutations(child.aa_muts, current_total_mutations)
                        # these are the mutations needed to make the structure equal to the first node
                        mutation_trunk_file.write("Root_seq" + "\t" + root_child_seq + "\t" + "Root_hash" + "\t" + root_child_hash + "\n")
                    node_foldx(self, child, "")
            else:  # internal or leaf node
                local_parent_mutation = str(current_total_mutations)
                local_parent_trunk = str(node.trunk)
                local_parent_hash = str(node)
                for child in node.child_nodes():
                    trunk = str(child.trunk)
                    tip = determine_tip(child)
                    hash = str(child)
                    current_total_mutations = update_mutations(child.aa_muts, local_parent_mutation)
                    mutation_trunk_file.write(local_parent_trunk + "\t" + local_parent_mutation[:len(local_parent_mutation) - 1] + "\t" + local_parent_hash + "\t" + trunk + "\t" + current_total_mutations[:len(current_total_mutations) - 1] + "\t" + hash + "\t" + str(tip) +"\n")
                    node_foldx(self, child, current_total_mutations)


        node_foldx(self, self.tree.seed_node, "")