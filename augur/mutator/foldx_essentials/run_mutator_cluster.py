import os, sys
from virus_stability import virus_stability

class run_mutator():

    def __init__(self, split_number, structure):
        try:
            self.sequence_file = open(split_number + "_sequences_file.txt", 'r')
        except:
            print("That sequence file does not exist in the current directory")
            raise
        self.sequence_to_ddG = {}  # dictionary from accession to virus object
        self.virus_list = []
        self.split_number = split_number
        self.structure = structure



    def __str__(self):
        return("Current mutation file is " + str(self.sequence_file))


    def read_sequence_file(self):
        '''
        reads through the sequence file and stores in self.sequence.list
        '''

        for line in self.sequence_file:
            self.virus_list.append(virus_stability_cluster("".join(line.strip()), self.split_number, self.structure))

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.virus_list:
            virus.calculate_ddg()

    def calculate_ddg_for_sequences(self):
        self.read_sequence_file()
        self.viruses_outgroup_ddG()

class virus_stability_cluster(virus_stability):
    def __init__(self, sequence, split_number, structure):
        virus_stability.__init__(self, None, None, None, None, None, sequence, "")
        self.split_number = split_number
        self.structure = structure

    def __str__(self):
        return "\t".join([self.ddg_outgroup[self.structure], "".join(self.mutations_from_outgroup)]) + "\n"

    def calculate_ddg(self):
        '''
        calls appropriate functions to calculate  ddG using foldx for the specified structure and ddG gets assigned to self.ddG_outgroup
        '''
        self.find_mutations(self.structure)
        if len(self.formatted_mut[self.structure]) > 0:
            self.overwrite_mutation_file(self.structure)
            if os.path.exists('foldx'):
                #os.system("./foldx --command=BuildModel --pdb=" + structure + "_trimer_repaired_" + self.outgroup + ".pdb --mutant-file=individual_list.txt")
                os.system("./foldx --command=BuildModel --pdb=" + self.structure + "_trimer_repaired.pdb --mutant-file=individual_list.txt")
            else:
                print("could not call foldx")
                raise FileNotFoundError
            self.read_ddG_output(self.structure)
            self.write_output()
        else:
            print("skipping this sequence, no valid mutations")
    def write_output(self):
        self.output_file_name = self.split_number + "_" + self.structure + "_sequences_ddg.txt"
        try:
            self.output_file = open(self.output_file_name, 'a')
        except:
            print("can't create output file in current directory")
            raise
        self.output_file.write(self.__str__())
        self.output_file.close()

def main(index, structure):
    os.chdir(index + "_foldx_split")
    os.system("chmod 755 foldx")
    new_stability_run = run_mutator(index, structure)
    new_stability_run.calculate_ddg_for_sequences()
    os.chdir("../")


if __name__ == "__main__":
    index = sys.argv[1]
    structure = sys.argv[2]
    main(index, structure)