import os, re, time, shutil
import dendropy
from seq_util import *
from date_util import *
from virus_stability import virus_stability

class tree_stability(object):





    def __init__(self, **kwargs):
        self.pdb_structures = ["1HA0", "2YP7"] # can add functionality to parser later
        fout_headings_begin = ["hash_code", "accession", "trunk_(T/F)", "tip_(T/F)"]
        fout_headings_end = ["parent accession", "aa_sequence"]
        fout_ddg_headings = []
        for pdb in self.pdb_structures:
            fout_ddg_headings.append("ddg_to_outgroup_" + pdb)
            fout_ddg_headings.append("ddg_to_parent_" + pdb)
        self.fout_headings = fout_headings_begin + fout_ddg_headings + fout_headings_end

        self.local_storage_name = "ddg_output_database.txt"
        self.ddg_calculated = {}
        self.directory = "foldx-output/"
        self.stability_output = "stability-data/"

        self.output_file_name = "0_ddG_mutations.txt"
        try:
            self.output_file = open(self.directory + self.output_file_name, 'w')
        except:
            print("can't create output file in current directory")
            raise

        self.new_sequences = []
        self.splits = 20
        self.run_on_cluster = True


    def calculate_stability(self):
        '''
        Calculate ddG between virus and outgroup and virus and parent
        '''
        self.viruses_outgroup_ddG()
        self.viruses_parent_ddG()

    def calculate_stability_cluster(self):
        print("Reading in current database of calculated stabilities for sequences")
        self.read_current_database()
        os.chdir(self.stability_output)
        print("Determining what new sequences will need stabilities calculated")
        list_sequences = self.split_list(self.get_new_sequences(), self.splits)
        print(len(list_sequences))
        if len(list_sequences) > 0:
            bash_file_name = "bash_file.sh"
            bash_file = open(bash_file_name, 'w')
            bash_file.write("#!/bin/sh" + "\n")

            for i in range(len(list_sequences)):
                self.make_sequence_file(i+1, list_sequences[i])
                self.add_bash_file(bash_file, i+1)
            bash_file.write("wait")
            bash_file.close()
        os.system("chmod 755 " + bash_file_name)
        os.system("./" + bash_file_name)
        os.chdir("../")
        print("Reading in new calculated stabilities for sequences")
        self.read_new_calc_ddg()
        self.viruses_outgroup_ddG()
        self.viruses_parent_ddG()
        self.print_viruses()
        self.write_current_database()

        #self.viruses_outgroup_ddG()
        #self.viruses_parent_ddG()

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
            virus.calculate_ddg_outgroup()
            #self.store_calculation(virus.seq, virus.ddg_outgroup['1HA0'], virus.ddg_outgroup['2YP7'])

    def viruses_parent_ddG(self):
        '''
        go through each virus object and calculate ddG from the parent to the current virus, also get mutations from parent to new virus
        '''
        for pair in self.virus_and_parent:
            virus = pair[0]
            parent = pair[1]
            virus.calculate_ddg_parent(parent)
            virus.get_parent_mutations(parent)

    def get_new_sequences(self):
        '''
        check if sequence of all viruses have been calculated yet, otherwise add them to list of new sequences to
        calculate ddg for.
        :return:
        '''
        new_sequences = []
        for virus in self.hash_to_virus.values():
            if virus.seq not in self.ddg_calculated.keys():
                new_sequences.append(virus.seq)
        self.new_sequences = new_sequences
        return new_sequences

    def split_list(self, initial, n):
        '''
        Split initial list into n lists
        :return list of split lists
        '''
        out = []
        new_index = int(1.0 * len(initial) / n + 0.5)
        for i in range(n-1):
            out.append(initial[i*new_index:i*new_index+new_index])
        out.append(initial[n*new_index-new_index:])
        return out

    def make_sequence_file(self, index, list_sequences):
        folder_name = str(index) + "_foldx_split"
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
            self.move_foldx_files(folder_name)

        sequence_file_name = str(index) + "_sequences_file.txt"
        sequence_file = open(folder_name + "/" +sequence_file_name, 'w')
        for sequence in list_sequences:
            sequence_file.write(sequence + "\n")

    def move_foldx_files(self, folder):
        '''
        :param destination_path: the path of where you want to move the files
        :return:moves the files from /src to the destination path
        '''
        essential_files_directory = os.getcwd() + "/foldx_essentials"
        essential_files = os.listdir(essential_files_directory)
        for file in essential_files:
            file_name = os.path.join(essential_files_directory, file)
            if os.path.isfile(file_name):
                shutil.copy(file_name, folder)

    def add_bash_file(self, bash_file, index):
        if self.run_on_cluster:
            bash_file.write("srun -n 1 -c 1 -t 24:00:00 -o output_" + str(index) + ".txt python " + str(index) + "_foldx_split/run_stability_cluster.py " + str(index) + " & \n")
        else:
            bash_file.write("python " + str(index) + "_foldx_split/run_stability_cluster.py " + str(index) + " & \n")

    def read_current_database(self):
        if os.path.isfile(self.directory + self.local_storage_name):
            read_local_storage_file = open(self.directory + self.local_storage_name, 'r')
            for line in read_local_storage_file:
                ddg_1HA0, ddg_2YP7, sequence = line.split("\t")
                self.ddg_calculated[sequence] = [ddg_1HA0, ddg_2YP7]
        else:
            print("No local ddG storage, creating new file")
            read_local_storage_file = open(self.directory + self.local_storage_name, 'w')
        read_local_storage_file.close()

    def read_new_calc_ddg(self):
        '''
        Go to each split folder and read in new calculated stability changes and sequences
        :return:
        '''
        for i in range(self.splits):
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
                self.ddg_calculated[sequence] = [ddg_1HA0, ddg_2YP7]


    def write_current_database(self):
        '''
        write to local storage file sequence and ddg
        :return:
        '''
        local_storage_file = open(self.directory + self.local_storage_name, 'w')
        for seq in self.ddg_calculated.keys():
            self.store_calculation(local_storage_file, self.ddg_calculated[seq], seq)
        local_storage_file.close()

    def store_calculation(self, file, ddg_list, sequence):
        file.write("\t".join([ddg_list[0], ddg_list[1], sequence]))
