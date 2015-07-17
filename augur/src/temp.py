def check_siterange_samesite(current_mutations, list_mutations, lowerRange, upperRange):
            mutation_dictionary = {}
            current_list = current_mutations.split(",")
            for mut in current_list:
                if len(mut) > 0:
                    site = int(mut[1:len(mut) - 1])
                    mutation_dictionary[site] = mut[0] + mut[len(mut) - 1]
            for mut in list_mutations:
                site = int(mut[1:len(mut) - 1])
                #if site >= lowerRange and site <= upperRange:
                if lowerRange <= site <= upperRange:
                    if site not in mutation_dictionary:
                        mutation_dictionary[site] = mut[0] + mut[len(mut) - 1]
                    else:
                        mutation_dictionary[site] = mutation_dictionary.get(site)[0] + mut[len(mut) - 1]
                        print("hi")
            '''
            print("hi")
            for site, mutation in mutation_dictionary.items():
                print(str(site) + " -> " + mutation)
            '''
            for site, mutation in mutation_dictionary.items():
                current_mutations += mutation[0] + str(site) + mutation[1] + ","
            print(current_mutations)
            return current_mutations


def node_foldx(self, node, current_mutations):
            if node.parent_node is None:  # root of the tree
                print("Root: " + node.aa_seq)
                for child in node.child_nodes():
                    print("Child Seq: " + child.aa_seq)
                    print("Child Muts: " + child.aa_muts)
                    if str(child.aa_muts) != "":
                        list_mutations = str(child.aa_muts).split(",")
                        current_mutations = check_siterange_samesite(current_mutations, list_mutations, lowerRange, upperRange)
                        # these are the mutations needed to make the structure equal to the first node
                        mutation_trunk_file.write("RootMut" + "\t" + current_mutations[:len(current_mutations) - 1] + ";\n")
                    node_foldx(self, child, "")
            else:  # internal or leaf node
                for child in node.child_nodes():
                        trunk = str(child.trunk)
                        # only want to use if mutation exists
                        if str(child.aa_muts) != "":
                            list_mutations = str(child.aa_muts).split(",")
                            print(list_mutations)
                            current_mutations = check_siterange_samesite(current_mutations, list_mutations, lowerRange, upperRange)
                        mutation_trunk_file.write(trunk + "\t" + current_mutations[:len(current_mutations) - 1] + ";\n")
                        node_foldx(self, child, current_mutations)


        node_foldx(self, self.tree.seed_node, "")






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

        # need to limit mutation sites to range of protein structure, 9-501. Also checks list of mutations so that only
        # one occurs at each site. So if in list had S9A, A9K. Would only include S9K in list.
        def check_siterange_samesite(current_mutations, list_mutations, lowerRange, upperRange):
            mutation_dictionary = {}
            #print(current_mutations)
            current_list = current_mutations.split(",")  # mutations that were passed to this node
            #print(current_list)
            if len(current_list) > 1:
                for mut in current_list:
                    if len(mut) > 0:
                        site = int(mut[1:len(mut) - 1])
                        mutation_dictionary[site] = mut[0] + mut[len(mut) - 1]
            #print(mutation_dictionary)
            #print(list_mutations + lowerRange + upperRange)
            '''
            for mut in list_mutations:
                site = int(mut[1:len(mut) - 1])
                #if site >= lowerRange and site <= upperRange:
                if lowerRange <= site <= upperRange:
                    if site not in mutation_dictionary:
                        mutation_dictionary[site] = mut[0] + mut[len(mut) - 1]
                    else:
                        mutation_dictionary[site] = mutation_dictionary.get(site)[0] + mut[len(mut) - 1]
            '''
            '''
            print("hi")
            for site, mutation in mutation_dictionary.items():
                print(str(site) + " -> " + mutation)
            '''
            for site, mutation in mutation_dictionary.items():
                current_mutations += mutation[0] + str(site) + mutation[1] + ","
            #print(mutation_dictionary)
            #print(current_mutations)
            return current_mutations


        def node_foldx(self, node, current_mutations):
            if node.parent_node is None:  # root of the tree
                print("Root: " + node.aa_seq)
                for child in node.child_nodes():
                    print("Child Seq: " + child.aa_seq)
                    print("Child Muts: " + child.aa_muts)
                    if str(child.aa_muts) != "":
                        list_mutations = str(child.aa_muts).split(",")
                        current_mutations = check_siterange_samesite(current_mutations, list_mutations, lowerRange, upperRange)
                        # these are the mutations needed to make the structure equal to the first node
                        mutation_trunk_file.write("RootMut" + "\t" + current_mutations[:len(current_mutations) - 1] + ";\n")
                    node_foldx(self, child, "")
            else:  # internal or leaf node
                for child in node.child_nodes():
                        trunk = str(child.trunk)
                        # only want to use if mutation exists
                        if str(child.aa_muts) != "":
                            list_mutations = str(child.aa_muts).split(",")
                            print(list_mutations)
                            current_mutations = check_siterange_samesite(current_mutations, list_mutations, lowerRange, upperRange)
                        mutation_trunk_file.write(trunk + "\t" + current_mutations[:len(current_mutations) - 1] + ";\n")
                        node_foldx(self, child, current_mutations)


        node_foldx(self, self.tree.seed_node, "")


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
            print(mut)
            site = mut[1:len(mut)]
            if lowerRange <= site <= upperRange:
                return mut + ","
            else:
                return ""


        def node_foldx(self, node, mutation):
            if node.parent_node is None:  # root of the tree
                for child in node.child_nodes():
                    #print(child.aa_seq)
                    #print(child.aa_muts)
                    '''
                    print(child.aa_muts)
                    this will print out all the mutations needed to get from the
                    root/structure to the first child nodes so that foldX should find all
                    mutations needed.
                    '''
                    mutation = "Mutation: "
                    if str(child.aa_muts) != "":
                        current_mutations = str(child.aa_muts).split(",")
                        #print(current_mutations)
                        for mut in current_mutations:
                            # need to limit to range of protein structure, 9-501
                            mutation += check_siterange(mut)
                            print(check_siterange(mut))
                            print(mutation)
                        print("Hi" + mutation)
                    mutation_trunk_file.write("RootMut" + "\t" + mutation[:len(mutation) - 1] + ";\n")
                    node_foldx(self, child, "")
            else:  # internal or leaf node
                for child in node.child_nodes():
                        trunk = str(child.trunk)
                        mutation = "Mutation: "
                        # for formatting purposes only want to add the new mutation to the list if it is not blank
                        if str(child.aa_muts) != "":
                            currentMutation = str(child.aa_muts).split(",")
                            for mut in currentMutation:
                                # need to limit to range of protein structure, 9-501
                                mutation += check_siterange(mut)
                        mutation_trunk_file.write(trunk + "\n" + mutation + ";\n")
                        node_foldx(self, child, mutation)


        node_foldx(self, self.tree.seed_node, "")