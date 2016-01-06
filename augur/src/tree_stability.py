import os, re, time, shutil
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability

class tree_stability(object):





    def __init__(self, **kwargs):
        self.pdb_structures = ["1HA0", "2YP7"] # can add functionality to parser later
        self.local_storage_name = "ddg_output_database.txt"
        self.ddg_calculated = {}
        self.directory = "foldx-output/"
        self.stability_output = "stability-data/"

        self.output_file_name = "ddg_output_database.txt"
        try:
            self.output_file = open(self.stability_output + self.output_file_name, 'w')
        except:
            print("can't create output file in current directory")
            raise
        self.num_splits = 20

    def calculate_stability(self):
        print("Reading in new calculated stabilities for sequences")
        #self.read_current_database()
        #self.read_new_calc_ddg()
        self.viruses_outgroup_ddG()
        self.viruses_parent_ddG()
        self.print_viruses()
        self.assign_node_ddG()
        #self.write_current_database()

    def print_viruses(self):
        for virus in self.hash_to_virus.values():
            print(virus)
            self.output_file.write(virus)
        self.output_file.close()

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.hash_to_virus.values():
            ddg_list = self.get_dynamodb_stability(virus.seq)
            virus.calculate_ddg_outgroup(ddg_list)

    def get_dynamodb_stability(self, sequence):
        '''
        checks the stability table to see if the sequence already has had stability calculated for it
        :return returns a list containing the stability output for that sequence, if it can't find the stability, raises an exception
        '''
        response = self.table.get_item(
            ProjectionExpression='ddg',
            Key={'sequence':sequence}
        )
        try:
            ddg_list = response['Item']['ddg']
            return ddg_list
        except:
            print("couldn't get ddg from the table for this sequence")
            raise

    def viruses_parent_ddG(self):
        '''
        go through each virus object and calculate ddG from the parent to the current virus, also get mutations from parent to new virus
        '''
        for pair in self.virus_and_parent:
            virus = pair[0]
            parent = pair[1]
            virus.calculate_ddg_parent(parent)
            virus.get_parent_mutations(parent)

    def assign_node_ddG(self):
        print("assigning ddg attribute to nodes")
        for node in self.tree.postorder_node_iter():
            hash = str(node)
            virus = self.hash_to_virus[hash]
            average_ddg = (virus.ddg_outgroup['1HA0'] + virus.ddg_outgroup['2YP7']) / 2
            try:
                setattr(node, 'ep', average_ddg)
            except:
                print("couldn't assign ddg attribute to current node")

    def read_current_database(self):
        '''
        read the updated current database and add to self.ddg_calculated if it's not there
        :return:
        '''
        if os.path.isfile(self.stability_output):
            read_local_storage_file = open(self.stability_output + self.local_storage_name, 'r')
            for line in read_local_storage_file:
                try:
                    ddg_1HA0, ddg_2YP7, sequence = line.split("\t")
                    if sequence not in self.ddg_calculated.keys():
                        self.ddg_calculated[sequence] = [ddg_1HA0, ddg_2YP7]
                except:
                    print("couldn't read this line")
                    print(line)
                    pass
        else:
            print("No local ddG storage, creating new file")
            read_local_storage_file = open(self.directory + self.local_storage_name, 'w')
        read_local_storage_file.close()

    def read_new_calc_ddg(self):
        '''
        Go to each split folder and read in new calculated stability changes and sequences
        :return:
        '''
        for i in range(self.num_splits):
            i += 1
            folder_name = str(i) + "_foldx_split"
            file_name = str(i) + "_sequences_ddg.txt"
            if os.path.exists(folder_name + "/" + file_name):
                new_file = open(folder_name + "/" + file_name, 'r')
            else:
                print("Couldn't find the split folders calculated ddg and sequences")
                raise FileNotFoundError
            for line in new_file:
                ddg_1HA0, ddg_2YP7, sequence = line.split("\t")
                if sequence not in self.ddg_calculated.keys():
                    self.ddg_calculated[sequence] = [ddg_1HA0, ddg_2YP7]


    def write_current_database(self):
        '''
        write to local storage file sequence and ddg
        :return:
        '''
        try:
            local_storage_file = open(self.stability_output + self.local_storage_name, 'a')
        except:
            local_storage_file = open(self.stability_output + self.local_storage_name, 'w')
            print("made a new database storage file")
        for seq in self.ddg_calculated.keys():
            self.store_calculation(local_storage_file, self.ddg_calculated[seq], seq)
        local_storage_file.close()

    def store_calculation(self, file, ddg_list, sequence):
        file.write("\t".join([ddg_list[0], ddg_list[1], sequence]))
