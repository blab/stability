import os, sys
from virus_stability import virus_stability
import boto3

class run_stability():

    def __init__(self, split_number):
        try:
            self.sequence_file = open(split_number + "_sequences_file.txt", 'r')
        except:
            print("That sequence file does not exist in the current directory")
            raise
        self.sequence_to_ddG = {}  # dictionary from accession to virus object
        self.virus_list = []
        self.split_number = split_number

    def __str__(self):
        return("Current mutation file is " + str(self.sequence_file))


    def read_sequence_file(self):
        '''
        reads through the sequence file and stores in self.sequence.list
        '''

        for line in self.sequence_file:
            self.virus_list.append(virus_stability_cluster(line.strip(), self.split_number))

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.virus_list:
            virus.calculate_ddg_outgroup()
            virus.upload_to_database()
            virus.output_file.write(virus.__str__())

    def calculate_ddg_for_sequences(self):
        self.read_sequence_file()
        self.viruses_outgroup_ddG()

class virus_stability_cluster(virus_stability):
    def __init__(self, sequence, split_number):
        virus_stability.__init__(self, None, None, None, None, None, sequence, "", "")
        self.output_file_name = split_number + "_sequences_ddg.txt"
        try:
            self.output_file = open(self.output_file_name, 'a')
        except:
            print("can't open output file in current directory")
            raise
        try:
            dynamodb = boto3.resource('dynamodb')
            self.table=dynamodb.Table('stability_1968')
        except:
            print("Couldn't connect to dynamodb or the stability table")
            raise

    def __str__(self):
        return "\t".join([self.ddg_outgroup["1HA0"], self.ddg_outgroup["2YP7"], self.seq]) + "\n"

    def calculate_ddg_outgroup(self):
        '''
        calls appropriate functions to calculate  ddG using foldx for the specified structure and ddG gets assigned to self.ddG_outgroup
        '''
        for structure in self.pdb_structures:
            self.align_to_outgroup()
            self.find_mutations()
            if len(self.mutations_from_outgroup) > 0:
                self.overwrite_mutation_file(structure)
                if os.path.exists('foldx'):
                    os.system("./foldx --command=BuildModel --pdb=" + structure + "_trimer_repaired_1968.pdb --mutant-file=individual_list.txt")
                else:
                    print("could not call foldx")
                    raise FileNotFoundError
                self.read_ddG_output(structure)


    def upload_to_database(self):
        '''
        upload to dynamodb 'stability' table, also write to local file
        :return:
        '''

        if len(self.mutations_from_outgroup) > 0:
            print("uploading calculation to dynamodb...")
            self.table.put_item(
               Item={
                   'sequence': self.seq,
                   'ddg_1968': [self.ddg_outgroup["1HA0"], self.ddg_outgroup["2YP7"]],
                }
            )
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