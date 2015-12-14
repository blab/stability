from virus_stability import virus_stability

class run_stability:

    def __init__(self, split_number):
        try:
            self.mutation_file = open("foldx_pipeline/" + split_number + "_mutation_trunk.txt", 'r')
        except:
            print("That mutation file does not exist in the current directory")
            raise
        self.accension_to_virus = {}  # dictionary from accension to virus object
        self.virus_and_parent = []  # set of lists; containing each virus and it's parent

    def __str__(self):
        return("Current mutation file is " + str(self.mutation_file))


    def read_mutation_file(self):
        '''
        reads through the mutation file and stores pair of virus and it's parent
        :mutates: appends new virus pairs to self.virus_and_parent
        '''

        for line in self.mutation_file:
            print(line.split("|"))
            v_info, pv_info = line.split("|")
            virus = self.check_virus_exists(v_info.strip())
            parent_virus = self.check_virus_exists(pv_info.strip())
            self.virus_and_parent.append([virus, parent_virus])
            print(virus)
            print(parent_virus)

    def check_virus_exists(self,virus_info):
        '''
        :param virus_info: string of tab sepearated information about the current virus
        :return: checks whether the current virus was already seen, returns that virus from self.accension_to_virus
        otherwise it returns a new object for the virus given
        '''

        try:
            print(virus_info.split("\t"))
            v_hash, v_accension, v_trunk, v_tip, v_seq = virus_info.split("\t")
        except:
            print("The mutation file line did not contain a hash code, accension number, trunk, tip and sequence.")
            raise

        if v_accension in self.accension_to_virus.keys():
            return self.accension_to_virus(v_accension)
        else:
            new_virus = virus_stability(v_hash, v_accension, v_trunk, v_tip, v_seq)
            self.accension_to_virus[v_accension] = new_virus
            return new_virus


def main():
    new_stability_run = run_stability("0")
    print(new_stability_run)
    new_stability_run.read_mutation_file()


if __name__=="__main__":
    main()