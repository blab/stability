import os, shutil

class calculate_stability(object):

    def __init__(self):
        self.new_sequences_list = []
        self.split_new_seq_list = []
        self.num_splits = 4
        self.bash_file_name = "bash_file.sh"

    def run(self):
        self.read_sequence_file()
        self.prep_split_folders()
        os.system("chmod 755 " + self.bash_file_name)
        os.system("./" + self.bash_file_name)

    def read_sequence_file(self):
        '''
        read the sequence file containing new sequences from tree_mutations.py. Then add those new sequences to self.new_sequences_list
        :return:
        '''
        try:
            sequence_file = open("new_seq_file.txt", 'r')
        except:
            print("couldn't find file that contains the new sequences")
            print(os.getcwd())
            raise
        for seq in sequence_file:
            if len(seq.strip()) > 0:
                self.new_sequences_list.append(seq.strip())
        self.split_new_seq_list = self.split_list(self.new_sequences_list, self.num_splits)

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

    def prep_split_folders(self):
        '''
        Create split folders and sequence files to run on the cluster and a bash file to call
        :return:
        '''
        if len(self.split_new_seq_list) > 0:
            bash_file = open(self.bash_file_name, 'w')
            bash_file.write("#!/bin/sh" + "\n")

            num = 1
            for list in self.split_new_seq_list:
                self.make_split_sequence_files(num, list)
                self.add_bash_run(bash_file, num)
                num+=1
            bash_file.write("wait")
            bash_file.close()
        else:
            print("no new sequences")
            raise Exception

    def add_bash_run(self, bash_file, index):
        '''
        Add a line to the bash file to calculate stabilities in the designated split
        :param bash_file: the bash file to write to
        :param index: the split to work with
        :return:
        '''
        bash_file.write("srun -n 1 -c 1 -t 120:00:00 -o output_" + str(index) + ".txt python " + str(index) + "_foldx_split/run_stability_cluster.py " + str(index) + " & \n")
        #bash_file.write("python " + str(index) + "_foldx_split/run_stability_cluster.py " + str(index) + " & \n")

    def move_foldx_files(self, folder):
        '''
        :param folder: the path of where you want to move the files
        :return:moves the files from /src to the folder
        '''
        essential_files_directory = os.getcwd() + "/foldx_essentials"
        essential_files = os.listdir(essential_files_directory)
        for file in essential_files:
            file_name = os.path.join(essential_files_directory, file)
            if os.path.isfile(file_name):
                shutil.copy(file_name, folder)

    def make_split_sequence_files(self, index, list_sequences):
        '''
        makes the split folder and moves all essential foldx files into it, then makes the sequence file
        :param index: identifying which split folder we are currently working with
        :param list_sequences: list of sequences to put in this split folder
        :return:
        '''
        folder_name = str(index) + "_foldx_split"
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
            self.move_foldx_files(folder_name)

        sequence_file_name = str(index) + "_sequences_file.txt"
        sequence_file = open(folder_name + "/" +sequence_file_name, 'w')
        for sequence in list_sequences:
            sequence_file.write(sequence + "\n")

def main():
    new_sequences = calculate_stability()
    new_sequences.run()

if __name__ == "__main__":
    main()