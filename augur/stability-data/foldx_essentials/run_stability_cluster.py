import os, sys
from virus_stability import virus_stability

class run_stability():

    def __init__(self, split_number):
        try:
            self.sequence_file = open(split_number + "_sequences_file.txt", 'r')
        except:
            print("That sequence file does not exist in the current directory")
            raise
        self.sequence_to_ddG = {}  # dictionary from accession to virus object
        self.virus_list = []

        self.output_file_name = split_number + "_sequences_ddg.txt"
        self.output_file_format_headings = ["sequence", "ddG to outgroup (1HA0)", "ddG to outgroup (2YP7)"]
        try:
            self.output_file = open(self.output_file_name, 'w')
        except:
            print("can't create output file in current directory")
            raise

    def __str__(self):
        return("Current mutation file is " + str(self.sequence_file))


    def read_sequence_file(self):
        '''
        reads through the sequence file and stores in self.sequence.list
        '''

        for line in self.sequence_file:
            self.virus_list.append(virus_stability_cluster(line.strip()))

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.virus_list:
            virus.calculate_ddg_outgroup()

    def write_to_txt(self):
        '''
        writes information for each virus to _ddG_mutations.txt with header on 1st line and tab de-limited
        '''
        for virus in self.virus_list:
            self.output_file.write(virus + "\n")

    def calculate_ddg_for_sequences(self):
        self.read_sequence_file()
        self.viruses_outgroup_ddG()
        self.write_to_txt()

class virus_stability_cluster(virus_stability):
    def __init__(self, sequence):
        virus_stability.__init__(self, None, None, None, None, None, sequence, "", "")

    def __str__(self):
        return "\t".join([self.ddg_outgroup["1HA0"], self.ddg_outgroup["2YP7"], self.seq])

    def calculate_ddg_outgroup(self):
        '''
        calls appropriate functions to calculate  ddG using foldx for the specified structure and ddG gets assigned to self.ddG_outgroup
        '''
        for structure in self.pdb_structures:
            self.align_to_outgroup()
            self.find_mutations()
            self.make_run_file(structure, "_trimer_repaired_1.pdb")
            self.overwrite_mutation_file(structure)
            if os.path.exists('foldx3b6'):
                os.system("./foldx3b6 -runfile mutate_runfile.txt")
            else:
                print("could not call foldx")
                raise FileNotFoundError
            self.read_ddG_output(structure)

def main(index):
    os.chdir(index + "_foldx_split")
    new_stability_run = run_stability(index)
    new_stability_run.calculate_ddg_for_sequences()
    os.chdir("../")


if __name__ == "__main__":
    index = sys.argv[1]
    main(index)