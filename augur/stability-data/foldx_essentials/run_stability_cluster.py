import os, sys
from virus_stability import virus_stability
import rethinkdb as r
import hashlib

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

class virus_stability_cluster(virus_stability):
    def __init__(self, sequence, split_number, structures):
        virus_stability.__init__(self, None, None, None, None, None, sequence, "")
        self.output_file_name = split_number + "_sequences_ddg.txt"
        try:
            self.output_file = open(self.output_file_name, 'a')
        except:
            print("can't open output file in current directory")
            raise
        self.pdb_structures = structures
        self.ddg_outgroup = {}

        self.database = 'test'
        self.table = 'stability'
        if 'RETHINK_AUTH_KEY' in os.environ:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")

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

        hash_function = hashlib.md5()
        hash_function.update(self.seq)
        hash_sequence = hash_function.hexdigest()
        self.ddg_outgroup['md5_sequence'] = hash_sequence

        for structure in self.pdb_structures:
            if len(self.structure_muts[structure]) > 0:
                print("uploading calculation to rethinkdb...")
                document = r.db(self.database).table(self.table).get(hash_sequence).run()
                # Sequence doesn't exist in table yet so add it
                if document is None:
                    print("Inserting new sequence document")
                    r.db(self.database).table(self.table).insert(self.ddg_outgroup).run()
                # Sequence exists in table so just add stability information
                else:
                    print("Updating existing sequence document")
                    r.db(self.database).table(self.table).get(hash_sequence).update({structure: self.ddg_outgroup[structure]}).run()
                print("upload successful")
            else:
                print("skipping this sequence, no valid mutations")

def main(index):
    os.chdir(index + "_foldx_split")
    new_stability_run = run_stability(index)
    new_stability_run.calculate_ddg_for_sequences()
    os.chdir("../")


if __name__ == "__main__":
    index = sys.argv[1]
    main(index)