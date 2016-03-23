import os, re, time, shutil, boto3
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability
import rethinkdb as r
import hashlib

class tree_stability(object):
    '''
    Goes back through all virus objects, looks up their calculated ddg from outgroup. Virus_stability determines all
    their other meta data. Prints all this information to /stability-data/ddg_output.txt. Assigns ddg to the 'ep'
    attribute of all nodes.
    '''
    def __init__(self, **kwargs):
        self.stability_output = "stability-data/"
        self.output_file_name = "ddg_output.txt"

        self.sequence_to_stability = {}

        if 'RETHINK_AUTH_KEY' in os.environ:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")
        self.database='test'
        self.table='stability'
        self.connect_rethink()

    def connect_rethink(self):
        '''
        Connect to rethink database,
        Check for existing table, otherwise create it
        '''
        try:
            r.connect(host="ec2-52-90-204-136.compute-1.amazonaws.com", port=28015, db=self.database, auth_key=self.auth_key).repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

    def calculate_stability(self):
        print("Reading in new calculated stabilities for sequences")
        self.viruses_ddG()
        self.print_viruses()
        self.assign_node_ddG()

    def print_viruses(self):
        print("Printing Viruses")
        for virus in self.hash_to_virus.values():
            print(virus.__str__())
            self.output_file.write(virus.__str__())
        self.output_file.close()

    def get_stability(self, sequence):
        '''
        checks the stability table to see if the sequence already has had stability calculated for it
        :return returns a list containing the stability output for that sequence, if it can't find the stability, raises an exception
        '''
        print("Getting Stability")
        hash_function = hashlib.md5()
        hash_function.update(sequence)
        hash_sequence = hash_function.hexdigest()
        document = r.table(self.table).get(hash_sequence).run()
        print(sequence)
        print(hash_sequence)
        print(document)
        '''
        if sequence in self.sequence_to_stability.keys():
            ddg_list = self.sequence_to_stability[sequence]
            return ddg_list
        else:
            ddg_list = self.get_dynamodb_stability(sequence)
            print("couldn't get ddg from the table for this sequence")
            print(sequence)
            return [0, 0]
        '''
    '''
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
    '''
    def viruses_ddG(self):
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
            print("-------")

            '''
            virus.align_to_outgroup()
            virus.calculate_ddg_outgroup(ddg_list)
            print(virus.ddg_outgroup.keys())
            virus.determine_relative_time()
            '''

    def sum_ddg(self, virus):
        '''
        sum up individual ddg mutation effects compared to each structure
        '''



    '''
    def viruses_parent_ddG(self):

        go through each virus object and calculate ddG from the parent to the current virus, also get mutations from parent to new virus

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
    '''


    '''
    def get_dynamodb_stability(self, sequence):
        response = self.table.get_item(
            ProjectionExpression='ddg_1968',
            Key={'sequence':sequence}
        )
        try:
            ddg_list = response['Item']['ddg_1968']
            self.sequence_to_stability[sequence] = ddg_list
            print(ddg_list)
            return ddg_list
        except:
            print("couldn't get ddg from the table for this sequence")
            print(sequence)
            return [0, 0]
            pass
    '''
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
    '''
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
    '''