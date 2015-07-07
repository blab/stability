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
		
        def node_foldx(self, node, mutation):
            if node.parent_node is None:  # root of the tree
                for child in node.child_nodes():
                    #print("hi")
                    print(child.aa_seq)
                    print(child.aa_muts)
                    '''
                    print(child.aa_muts) 
                    this will print out all the mutations needed to get from the 
                    root/structure to the first child nodes so that foldX should find all 
                    mutations needed. 
                    '''
                    node_foldx(self, child, "")
            else:  # internal or leaf node
                for child in node.child_nodes():
                        trunk = str(child.trunk)
                        # for formatting purposes only want to add the new mutation to the list if it is not blank
                        if str(child.aa_muts) != "":
                            currentMutation = str(child.aa_muts).split(",")
                            for mut in currentMutation:
                                # need to limit to range of protein structure, 9-501
                                site = mut[1:len(mut)]
                                if site >= lowerRange and site <= upperRange:
                                	mutation += mut + ","
                        print(mutation)
                        mutation_trunk_file.write(trunk + "\n" + mutation + ";\n")
                        node_foldx(self, child, mutation)


        node_foldx(self, self.tree.seed_node, "")