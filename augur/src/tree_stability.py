import os, re, time, shutil, boto3
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability

class tree_stability(object):
    '''
    Goes back through all virus objects, looks up their calculated ddg from outgroup. Virus_stability determines all
    their other meta data. Prints all this information to /stability-data/ddg_output.txt. Assigns ddg to the 'ep'
    attribute of all nodes.
    '''
    def __init__(self, **kwargs):
        self.pdb_structures = ["1HA0", "2YP7"] # can add functionality to parser later

        self.stability_output = "stability-data/"
        self.output_file_name = "ddg_output.txt"
        self.ddg_database_file_name = "ddg_database.txt"

        self.local_database_exists = False

        try:
            self.output_file = open(self.stability_output + self.output_file_name, 'w')
            if not self.local_database_exists:
                self.ddg_database_file = open(self.stability_output + self.ddg_database_file_name, 'w')
        except:
            print("can't create output file in current directory")
            raise
        self.sequence_to_stability = {}

        dynamodb = boto3.resource('dynamodb')
        self.table=dynamodb.Table('stability_1968')

    def calculate_stability(self):
        print("Reading in new calculated stabilities for sequences")
        if not self.local_database_exists:
            self.load_local_database()
        else:
            self.read_local_database()
            self.viruses_outgroup_ddG()
            self.viruses_parent_ddG()
            self.print_viruses()
            self.assign_node_ddG()

    def print_viruses(self):
        print("Printing Viruses")
        for virus in self.hash_to_virus.values():
            print(virus.__str__())
            self.output_file.write(virus.__str__())
        self.output_file.close()

    def viruses_ddG(self):
        num = 0
        for pair in self.virus_and_parent:
            num += 1
            print("virus " + str(num) + " of " + str(len(self.virus_and_parent)) + " viruses")
            virus = pair[0]
            parent = pair[1]
            print("Virus: ")
            virus.calculate_ddg_outgroup(self.get_stability(virus.seq))
            print("Parent: ")
            parent.calculate_ddg_outgroup(self.get_stability(parent.seq))
            virus.calculate_ddg_parent(parent)

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        num = 0
        for virus in self.hash_to_virus.values():
            num += 1
            print("virus " + str(num) + " of " + str(len(self.hash_to_virus.keys())) + " viruses")
            ddg_list = self.get_stability(virus.seq)
            print(ddg_list)
            virus.align_to_outgroup()
            virus.calculate_ddg_outgroup(ddg_list)
            print(virus.ddg_outgroup.keys())
            virus.determine_relative_time()

    def viruses_parent_ddG(self):
        '''
        go through each virus object and calculate ddG from the parent to the current virus, also get mutations from parent to new virus
        '''
        print("Calculating parent to virus ddg")
        print(len(self.virus_and_parent))
        for pair in self.virus_and_parent:
            virus = pair[0]
            parent = pair[1]
            try:
                virus_ddg = self.sequence_to_stability[virus.seq]
            except:
                virus_ddg = ['0.0', '0.0']
                print("Couldn't find in dictionary")
                print(virus.seq)
            try:
                parent_ddg = self.sequence_to_stability[parent.seq]
            except:
                parent_ddg = ['0.0', '0.0']
                print("Couldn't find in dictionary")
                print(parent.seq)
            virus.parent_strain = parent.strain
            virus.calculate_ddg_parent(virus_ddg, parent_ddg)
            virus.get_parent_mutations(virus.mutations_from_outgroup, parent.mutations_from_outgroup)

    def get_stability(self, sequence):
        '''
        checks the stability table to see if the sequence already has had stability calculated for it
        :return returns a list containing the stability output for that sequence, if it can't find the stability, raises an exception
        '''
        if sequence in self.sequence_to_stability.keys():
            ddg_list = self.sequence_to_stability[sequence]
            return ddg_list
        else:
            ddg_list = self.get_dynamodb_stability(sequence)
            print("couldn't get ddg from the table for this sequence")
            print(sequence)
            return [0, 0]


    def get_dynamodb_stability(self, sequence):
        response = self.table.get_item(
            ProjectionExpression='ddg',
            Key={'sequence':sequence}
        )
        try:
            ddg_list = response['Item']['ddg']
            self.sequence_to_stability[sequence] = ddg_list
            print(ddg_list)
            return ddg_list
        except:
            print("couldn't get ddg from the table for this sequence")
            print(sequence)
            return [0, 0]
            pass

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

    def read_local_database(self):
        if os.path.exists(self.stability_output+self.ddg_database_file_name):
            print("found local database")
            self.ddg_database_file = open(self.stability_output + self.ddg_database_file_name, 'r')
            print(self.stability_output + self.ddg_database_file_name)
            for line in self.ddg_database_file:
                split_line = line.split("\t")
                if len(split_line) == 3:
                    try:
                        sequence = split_line[2].strip()
                        self.sequence_to_stability[sequence] = [split_line[0].strip(), split_line[1].strip()]
                    except:
                        raise
            self.sequence_to_stability['MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGRICDSPHRILDGKNCTLIDALLGDPHCDGFQNKEWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFINEDFNWTGVAQDGGSYACKRGSVNSFFSRLNWLHKSEYKYPALNVTMPNNGKFDKLYIWGVHHPSTDRDQTSLYVRASGRVTVSTKRSQQTVTPNIGSRPWVRGQSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRLIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQKGNIRCNICI'] = ['0.0', '0.0']
            print("Number of sequences in dictionary")
            print(len(self.sequence_to_stability.keys()))
        else:
            print("couldn't find local database")
            print(self.stability_output+self.ddg_database_file_name)
            raise Exception

    def load_local_database(self):
        num = 0
        for virus in self.hash_to_virus.values():
            num += 1
            print("virus " + str(num) + " of " + str(len(self.hash_to_virus.keys())) + " viruses")
            self.get_dynamodb_stability(virus.seq)
        self.print_database()

    def print_database(self):
        print("Printing out the database")
        for sequence in self.sequence_to_stability:
            ddg_list = self.sequence_to_stability[sequence]
            self.ddg_database_file.write(ddg_list[0] + "\t" + ddg_list[1] + "\t" + sequence + "\n")