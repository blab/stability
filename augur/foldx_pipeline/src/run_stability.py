from virus_stability import virus_stability

class run_stability:

    def __init__(self, split_number):
        try:
            self.mutation_file = open("foldx_pipeline/" + split_number + "_mutation_file.txt", 'r')
        except:
            print("That mutation file does not exist in the current directory")
            raise
        self.accession_to_virus = {}  # dictionary from accession to virus object
        self.virus_and_parent = []  # set of lists; containing each virus and it's parent

        self.output_file_name = split_number + "_" + "ddG_mutations.txt"
        self.output_file_format_headings = ["hash code", "accession", "trunk (T/F)", "tip (T/F)", "ddG to outgroup (1HA0)", "ddG to outgroup (2YP7)", "ddG to parent (1HA0)", "ddG to parent (2YP7)", "parent accession", "aa_sequence"]
        try:
            self.output_file = open(self.output_file_name, 'w')
        except:
            print("can't create output file in current directory")
            raise

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
        :param virus_info: string of tab separated information about the current virus
        :return: checks whether the current virus was already seen, returns that virus from self.accession_to_virus
        otherwise it returns a new object for the virus given
        '''

        try:
            print(virus_info.split("\t"))
            v_hash, v_accession, v_trunk, v_tip, v_seq = virus_info.split("\t")
        except:
            print("The mutation file line did not contain a hash code, accession number, trunk, tip and sequence.")
            raise

        if v_accession in self.accession_to_virus.keys():
            return self.accession_to_virus[v_accession]
        else:
            new_virus = virus_stability(v_hash, v_accession, v_trunk, v_tip, v_seq)
            self.accession_to_virus[v_accession] = new_virus
            return new_virus

    def viruses_outgroup_ddG(self):
        '''
        go through each virus object and determine the list of foldx formatted mutations for each structure. Also calculate
        the ddG from the outgroup to the current virus for each structure
        '''
        for virus in self.accession_to_virus.values():
            virus.find_mutations()
            virus.calculate_ddG_outgroup("1HA0")
            virus.calculate_ddG_outgroup("2YP7")

    def viruses_parent_ddG(self):
        '''
        go through each virus object and calculate ddG from the parent to the current virus
        '''
        for pair in self.virus_and_parent:
            virus = pair[0]
            parent = pair[1]
            virus.calculate_ddG_parent(parent)

    def write_to_txt(self):
        '''
        writes information for each virus to _ddG_mutations.txt with header on 1st line and tab de-limited
        '''
        self.output_file.write("\t".join(self.output_file_format_headings))
        for virus in self.accession_to_virus:
            self.output_file.write("\t".join(virus.output_file_format()) + "\n")
'''
    def write_to_database(self):
        # store in database
'''
def main():
    new_stability_run = run_stability("0")
    new_stability_run.read_mutation_file()
    new_stability_run.viruses_outgroup_ddG()
    new_stability_run.viruses_parent_ddG()
    new_stability_run.write_to_txt()


if __name__=="__main__":
    main()