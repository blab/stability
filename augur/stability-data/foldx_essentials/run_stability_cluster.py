import os, sys
import rethinkdb as r
import hashlib
from mutation_stability import mutation_stability

class run_stability():

    def __init__(self, split_number):
        try:
            self.sequence_file = open(split_number + "_sequences_file.txt", 'r')
        except:
            print("That sequence file does not exist in the current directory")
            raise
        self.sequence_to_ddG = {}  # dictionary from accession to virus object
        self.viruses = {}
        self.virus_list = []
        self.split_number = split_number

    def read_sequence_file(self):
        '''
        reads through the sequence file and stores in self.sequence.list
        '''
        for line in self.sequence_file:
            split_line = line.split("\t")
            structures = split_line[0].split(",")
            sequence = split_line[1]
            self.virus_list.append(virus_stability_cluster(sequence, self.split_number, structures))

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.virus_list:
            virus.calculate_ddg()
            virus.upload_to_database()
            virus.output_file.write(virus.__str__())

    def calculate_ddg_for_sequences(self):
        self.read_sequence_file()
        self.viruses_outgroup_ddG()

class virus_stability_cluster():
    def __init__(self, sequence, split_number, structures):
        self.output_file_name = split_number + "_sequences_ddg.txt"
        try:
            self.output_file = open(self.output_file_name, 'a')
        except:
            print("can't open output file in current directory")
            raise
        self.pdb_structures = structures
        self.seq = sequence
        self.ddg_outgroup = {}
        self.formatted_mut = {}
        self.structure_muts = {}

        self.database = 'test'
        self.table = 'stability'
        if 'RETHINK_AUTH_KEY' in os.environ:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")

        self.structure_seqs = {}
        self.structure_seqs['1HA0'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"
        self.structure_seqs['2YP7'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRKSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"
        self.structure_seqs['2YP2'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"

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
            raise Exception("Failed to connect to the database" + self.database)

    def __str__(self):
        return "\t".join([str(self.ddg_outgroup), self.seq]) + "\n"

    def calculate_ddg(self):
        '''
        calls appropriate functions to calculate  ddG using foldx for the specified structure and ddG gets assigned to self.ddG_outgroup
        '''
        for structure in self.pdb_structures:
            self.find_mutations(structure)
            if len(self.formatted_mut[structure]) > 0:
                self.overwrite_mutation_file(structure)
                if os.path.exists('foldx'):
                    #os.system("./foldx --command=BuildModel --pdb=" + structure + "_trimer_repaired_" + self.outgroup + ".pdb --mutant-file=individual_list.txt")
                    os.system("./foldx --command=BuildModel --pdb=" + structure + "_trimer_repaired.pdb --mutant-file=individual_list.txt")
                else:
                    print("could not call foldx")
                    raise FileNotFoundError
                self.read_ddG_output(structure)

    def find_mutations(self, structure):
        '''
        Finds and stores mutations that are valid for each of the structures specified
        '''

        self.align_to_structure(structure)
        #list_of_mutations = list(self.mutations_from_outgroup)
        list_of_mutations = self.structure_muts[structure]
        mut_stability = mutation_stability(list_of_mutations, structure)
        self.formatted_mut[structure] = mut_stability.get_formatted_mutations()

    def align_to_structure(self, structure):
        '''
        aligns to structure sequence to virus sequence
        :return: mutations from structure to self.seq
        '''

        mutations_set = set()
        structure_align_seq = self.structure_seqs[structure][24:]
        virus_align_seq = self.seq[24:]
        if (len(structure_align_seq)>virus_align_seq):
            print("Outgroup Sequence longer than the virus sequence")
            raise Exception
        for index in range(len(structure_align_seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")
            if structure_align_seq[index] != virus_align_seq[index]:
                mutation = structure_align_seq[index] + str(site) + virus_align_seq[index]
                mutations_set.add(mutation)
        self.structure_muts[structure] = list(mutations_set)

    def overwrite_mutation_file(self, structure):
        '''
        creates or overwrites "mutation_runfile_list.txt" which includes foldx formated, comma-seperated list of mutations for current virus
        :param structure: specify structure to be used in order to use right list of mutations
        '''
        try:
            mutations = self.formatted_mut[structure]
        except:
            print("mutations were not formatted for this structure")
            raise

        mutationFileName = "individual_list.txt"
        mutationFile = open(mutationFileName, 'w')  # overwrites the current file for FoldX
        mutationFile.write(mutations + ";")
        mutationFile.close()

    def read_ddG_output(self, structure):
        '''
        opens the output of the mutation command in foldX and gets the ddG value for the mutation that was just performed
        :param structure: specify the structure that was used by foldx
        '''
        ddGFileName = "Average_" + structure + "_trimer_repaired.fxout"
        ddGFile = open(ddGFileName, 'r')
        try:
            for line in ddGFile:
                if line.startswith(structure):
                    ddGline = line.split()
                    ddG = ddGline[2]
        except:
            print("couldn't find ddG output")
            raise
        ddGFile.close()
        self.ddg_outgroup[structure] = ddG
        os.remove(ddGFileName)

    def upload_to_database(self):
        '''
        upload to dynamodb 'stability' table, also write to local file
        :return:
        '''
        print("uploading calculation to rethinkdb...")
        hash_function = hashlib.md5()
        hash_function.update(self.seq)
        hash_sequence = hash_function.hexdigest()
        document = r.table(self.table).get(hash_sequence).run()
        # Sequence doesn't exist in table yet so add it
        if document is None:
            print("Inserting new sequence document")
            self.ddg_outgroup['md5'] = hash_sequence
            r.db(self.database).table(self.table).insert(self.ddg_outgroup).run()
        else:
            print("Updating existing sequence document")
            for structure in self.pdb_structures:
                r.db(self.database).table(self.table).get({'md5': hash_sequence}).update({structure: self.ddg_outgroup[structure]}).run()
        print("upload successful")

def main(index):
    os.chdir(index + "_foldx_split")
    new_stability_run = run_stability(index)
    new_stability_run.calculate_ddg_for_sequences()
    os.chdir("../")


if __name__ == "__main__":
    index = sys.argv[1]
    main(index)