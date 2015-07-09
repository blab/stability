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
        print to the document mutation_trunk.txt with trunk information on the first line
        followed by mutation information on the second line

        '''
        mutation_trunk_fileName = "mutation_trunk.txt"
        mutation_trunk_file = open(mutation_trunk_fileName, 'w')

        # change depending on protein structure/outgroup you are using
        lowerRange = 9
        upperRange = 501


        def check_siterange(mut):
            site = int(mut[1:len(mut) - 1])
            #print(mut + "->" + str(site))
            if lowerRange <= site <= upperRange:
                return mut + ","
            else:
                return ""

        def update_mutations(current_node_mutations, current_total_mutations):
            # for formatting purposes only want to add the new mutation to the list if it is not blank
            if str(current_node_mutations) != "":
                currentMutation = str(current_node_mutations).split(",")
                for mut in currentMutation:
                    current_total_mutations += check_siterange(mut)
            return current_total_mutations

        def node_foldx(self, node, mutation):
            if node.parent_node is None:  # root of the tree
                for child in node.child_nodes():
                    #print("hi")
                    #print(child.aa_seq)
                    #print(child.aa_muts)
                    mutation = update_mutations(child.aa_muts, mutation)
                    # these are the mutations needed to make the structure equal to the first node
                    mutation_trunk_file.write("RootMut" + "\t" + mutation[:len(mutation) - 1] + ";\n")
                    print("RootMut" + "\t" + mutation[:len(mutation) - 1] + ";\n")
                    node_foldx(self, child, "")
            else:  # internal or leaf node
                for child in node.child_nodes():
                    trunk = str(child.trunk)
                    # for formatting purposes only want to add the new mutation to the list if it is not blank
                    mutation = update_mutations(child.aa_muts, mutation)
                    print(trunk + "\t" + mutation[:len(mutation) - 1] +";\n")
                    mutation_trunk_file.write(trunk + "\t" + mutation[:len(mutation) - 1] +";\n")
                    node_foldx(self, child, mutation)


        node_foldx(self, self.tree.seed_node, "")