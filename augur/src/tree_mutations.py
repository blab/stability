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

        '''
        mutation_trunk_fileName = "mutation_trunk.txt"
        mutation_trunk_file = open(mutation_trunk_fileName, 'w')

        def node_foldx(self, node, mutation):
            if node.parent_node is None:  # root of the tree
                    for child in node.child_nodes():
                        node_foldx(self, child, "")
            else:  # internal or leaf node
                for child in node.child_nodes():
                        trunk = str(child.trunk)
                        # for formatting purposes only want to add the new mutation to the list if it is not blank
                        if str(child.aa_muts) != "":
                            currentMutation = str(child.aa_muts).split(",")
                            #if len(currentMutation) > 0:
                            for mut in currentMutation:
                                mutation += mut + ";"
                        print(trunk)
                        print(mutation)
                        mutation_trunk_file.write(trunk + "\n" + mutation + "\n")
                        #mutation_trunk_file.write(mutation)
                        node_foldx(self, child, mutation)


        node_foldx(self, self.tree.seed_node, "")