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
        mutation_trunk_fileName = "mutation_trunk.txt"
        mutation_trunk_file = open(mutation_trunk_fileName, 'w')

        # Change depending on protein structure/outgroup you are using
        '''
        # 4WE4 1968 Hong Kong protein structure
        lowerRange = 9
        upperRange = 501
        '''
        '''
        # 1HA0 1968 Aichi protein structure
        lowerRange = 9
        upperRange = 502
        '''
        '''
        # 4WE5
        lowerRange = 4
        upperRange = 496
        '''
        '''
        # 4WE6
        lowerRange = 38
        upperRange = 580
        '''
        '''
        # 4WE9
        lowerRange = 8
        upperRange = 502
        '''
        
        # 2YP7
        lowerRange = 8
        upperRange = 503


        # Need to limit mutation sites to range of protein structure, so for 4WE4 to sites 9-501.
        def check_siterange(mut):
            site = int(mut[1:len(mut) - 1])
            #print(mut + "->" + str(site))
            if lowerRange <= site <= upperRange:
                return mut + ","
            else:
                return ""

        # Checks list of mutations so that only one occurs at each site.
        # So if in list had S9A, A9K. Would only include S9K in list.
        def check_multiple_mutations(mutations):
            mutation_dictionary = {}
            complete_mutations = ""
            mutation_list = mutations.split(",")
            if len(mutation_list) > 1:
                for mut in mutation_list:
                    if len(mut) > 0:
                        site = int(mut[1:len(mut) - 1])
                        if site not in mutation_dictionary:
                            mutation_dictionary[site] = mut[0] + mut[len(mut) - 1]
                        else:
                            mutation_dictionary[site] = mutation_dictionary.get(site)[0] + mut[len(mut) - 1]
            for site, mutation in mutation_dictionary.items():
                complete_mutations += mutation[0] + str(site) + mutation[1] + ","
            return complete_mutations
            #print(mutations)

        # Add the current nodes mutations to the total list of mutations that were needed to get to that
        # node from the root
        def update_mutations(current_node_mutations, current_total_mutations):
            # for formatting purposes only want to add the new mutation to the list if it is not blank
            if str(current_node_mutations) != "":
                currentMutation = str(current_node_mutations).split(",")
                for mut in currentMutation:
                    current_total_mutations += check_siterange(mut)
            current_total_mutations = check_multiple_mutations(current_total_mutations)
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
                print("Root : " + node.aa_seq)
                for child in node.child_nodes():
                    print("Child: " + child.aa_seq)
                    print("Child: " + child.aa_muts)
                    if child.aa_muts != "":  # only want to print out root mutations if not the outgroup
                        current_total_mutations = ""
                        current_total_mutations = update_mutations(child.aa_muts, current_total_mutations)
                        # these are the mutations needed to make the structure equal to the first node
                        mutation_trunk_file.write("RootMut" + "\t" + current_total_mutations[:len(current_total_mutations) - 1] + "\n")
                    node_foldx(self, child, "")
            else:  # internal or leaf node
                local_parent_mutation = str(current_total_mutations)
                local_parent_trunk = str(node.trunk)
                for child in node.child_nodes():
                    trunk = str(child.trunk)
                    tip = determine_tip(child)
                    current_total_mutations = update_mutations(child.aa_muts, local_parent_mutation)
                    mutation_trunk_file.write(local_parent_trunk + "\t" + local_parent_mutation[:len(local_parent_mutation) - 1] + "\t" + trunk + "\t" + current_total_mutations[:len(current_total_mutations) - 1] + "\t" + str(tip) +"\n")
                    node_foldx(self, child, current_total_mutations)


        node_foldx(self, self.tree.seed_node, "")