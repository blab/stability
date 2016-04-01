import os, re, time, shutil, json
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
        self.output_file = open(self.stability_output + "ddg_output.txt", 'w')

        self.sequence_to_stability = {}

        if 'RETHINK_AUTH_KEY' in os.environ:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")
        self.database = 'test'
        self.table = 'stability'
        self.connect_rethink()

        self.structure_seqs = {}
        self.structure_seqs['1HA0'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"
        self.structure_seqs['2YP7'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRKSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"
        self.structure_seqs['2YP2'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"

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
        self.sum_ddg()
        #self.epistasis_ddG()
        self.print_viruses()
        #self.assign_node_ddG()

    def print_viruses(self):
        print("Printing Viruses")
        json.dump(self.hash_to_virus, self.output_file, indent=1)
        self.output_file.close()

    def epistasis_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.hash_to_virus.values():
            ddg_list = self.get_stability(virus)

    def get_stability(self, virus):
        '''
        checks the stability table to see if the sequence already has had stability calculated for it
        :return returns a list containing the stability output for that sequence, if it can't find the stability, raises an exception
        '''
        sequence = virus['seq']
        hash_function = hashlib.md5()
        hash_function.update(sequence)
        hash_sequence = hash_function.hexdigest()
        document = r.table(self.table).get(hash_sequence).run()
        if document is not None:
            for structure in self.structures:
                virus[structure]['epistasis_ddg'] = document[structure]
        else:
            raise Exception("Couldn't find this sequence in rethinkdb")

    def sum_ddg(self):
        '''
        sum up individual ddg mutation effects compared to each structure
        '''
        print("Determining stability change by summing individual mutation effects on each structure")
        self.open_mutator()
        for virus in self.hash_to_virus.values():
            for structure in self.structures:
                virus[structure] = self.align_to_structure(virus, structure)
                ddg = 0
                for mut in virus[structure]['mutations']:
                    ddg += self.mutator_ddg[structure][mut]
                virus[structure]['sum_ddg'] = ddg

    def open_mutator(self):
        self.mutator_ddg = {}
        for structure in self.structures:
            file = open("source-data/" + structure + "_mutator_ddg.txt", 'r')
            mut_ddg = {}
            for line in file:
                info = line.split()
                mut_ddg[info[1]] = float(info[0])
            self.mutator_ddg[structure] = mut_ddg

    def align_to_structure(self, virus, structure):
        '''
        aligns to structure sequence to virus sequence
        :return: mutations from structure to self.seq
        '''
        mutations = []
        structure_align_seq = self.structure_seqs[structure][24:]
        virus_align_seq = virus['seq'][24:]
        if (len(structure_align_seq)>virus_align_seq):
            print("Outgroup Sequence longer than the virus sequence")
            raise Exception
        for index in range(len(structure_align_seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")
            if structure_align_seq[index] != virus_align_seq[index]:
                mutation = structure_align_seq[index] + str(site) + virus_align_seq[index]
                mutations.append(mutation)
        mutations = filter(lambda mut: self.site_range_valid(mut), mutations)
        return {'mutations': mutations}
        #self.structure_muts[structure] = list(mutations_set)

    def site_range_valid(self, mutation):
       '''
       protein structures (1HA0, 2YP7) are missing certain amino acid sites, method checks that mutation is in structure
       :param mutation: mutation in standard format
       :return: true if site is in structure, false if site range is not in structure
       '''
       lowerRange = 9
       upperRange = 502
       missing_lower = 328
       missing_upper = 333
       site = int(mutation[1:len(mutation) - 1])
       if missing_lower <= site <= missing_upper:  # in missing middle section
            return False
       elif lowerRange <= site <= upperRange: # in range of protein structure besides middle section
            return True
       else:
            return False

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