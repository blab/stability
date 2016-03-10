from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mutation_stability import mutation_stability
import os



class virus_stability(object):


    def __init__(self, hash, strain, trunk, tip, date, seq, foldx_directory):
        self.hash_code = hash
        self.strain = strain
        self.trunk = trunk
        self.tip = tip
        self.date = date
        self.seq = seq
        self.foldx_directory = foldx_directory

        self.mutations_from_outgroup = set()  #does not include chain info
        self.mutations_from_1968 = set()
        self.mutations_time = 0

        self.mutations_from_parent = set()

        self.pdb_structures = []
        self.formatted_mut = {}

        self.ddg_outgroup = {}
        self.ddg_parent= {}
        self.parent_strain = ""

        self.structure_seqs = {}
        self.structure_seqs['1HA0'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"
        self.structure_seqs['2YP7'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRKSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"
        self.structure_seqs['2YP2'] = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"

        self.structure_muts = {}



    def __str__(self):
        # ["hash code", "strain", "trunk (T/F)", "tip (T/F)", date, "ddG to outgroup (1HA0)", "ddG to outgroup (2YP7)", "ddG to parent (1HA0)", "ddG to parent (2YP7)", "mutation from parent" "parent strain", "aa_sequence"]
        if len(self.ddg_outgroup.keys()) == 0:
            print("Readjusting parent ddg!!!")
            self.ddg_outgroup['1HA0'] = None
            self.ddg_outgroup['2YP7'] = None
        if len(self.ddg_parent.keys()) == 0:
            print("Readjusting parent ddg!!!")
            self.ddg_parent['1HA0'] = None
            self.ddg_parent['2YP7'] = None
        return "\t".join([self.hash_code, str(self.strain), str(self.trunk), str(self.tip), str(self.date), str(self.ddg_outgroup['1HA0']), str(self.ddg_outgroup['2YP7']), str(self.ddg_parent['1HA0']), str(self.ddg_parent['2YP7']), str(" ".join(list(self.mutations_from_parent))), str(" ".join(list(self.mutations_from_1968))), str(self.mutations_time), str(self.parent_strain), self.seq, "\n"])

    '''
    def find_outgroup(self, directory):
        # get outgroup amino_acid sequence

        try:
            #print("Using H3N2_outgroup.gb for the outgroup, assumed to be Beijing 1992 strain")
            temp_outgroup = SeqIO.read(directory + 'H3N2_outgroup.gb', 'genbank')
        except:
            print("could not find H3N2_outgroup.gb which contained the outgroup file")
            raise
        dna_seq = str(temp_outgroup.seq).upper()
        coding_dna = Seq(dna_seq, generic_dna)
        protein = coding_dna.translate()
        self.outgroup_seq = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"

    '''
    '''
    def align_to_outgroup(self):
        :mutates: self.mutations_from_outgroup
        aligns outgroup and self.seq, returns mutations from outgroup -> self.seq. Modifies self.mutations_from_outgroup
        :return: comma-seperated string of mutations from outgroup

        mutations_set = set()
        outgroup_align_seq = self.outgroup_seq[24:]
        virus_align_seq = self.seq[24:]
        if (len(outgroup_align_seq)>virus_align_seq):
            print("Outgroup Sequence longer than the virus sequence")
            raise Exception
        for index in range(len(outgroup_align_seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")
            if outgroup_align_seq[index] != virus_align_seq[index]:
                mutation = outgroup_align_seq[index] + str(site) + virus_align_seq[index]
                mutations_set.add(mutation)
        self.mutations_from_outgroup = mutations_set
        return ','.join(mutations_set)
    '''
    '''
    def align_outgroup_to_sequence(self):

        :mutates: self.mutations_from_outgroup
        aligns self.seq and outgroup, returns mutations from sequence -> outgroup
        :return: comma-seperated string of mutations from structure
        mutations_set = set()
        outgroup_align_seq = self.outgroup_seq[24:]
        align_seq = self.seq[24:]
        if (len(outgroup_align_seq)>self.seq):
            print("Outgroup Sequence longer than the virus sequence")
            raise Exception
        for index in range(len(outgroup_align_seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")

            if outgroup_align_seq[index] != align_seq[index]:
                mutation = align_seq[index] + str(site) + outgroup_align_seq[index]
                print(mutation)
                mutations_set.add(mutation)
        self.mutations_from_outgroup = mutations_set
        return ','.join(mutations_set)
    '''
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
        if (len(structure_align_seq)>len(virus_align_seq)):
            raise Exception("structure sequence longer than the virus sequence")
        for index in range(len(structure_align_seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")
            if structure_align_seq[index] != virus_align_seq[index]:
                print(structure_align_seq[index] + str(site) + virus_align_seq[index])
                mutation = structure_align_seq[index] + str(site) + virus_align_seq[index]
                mutations_set.add(mutation)
        self.structure_muts[structure] = list(mutations_set)

    def get_parent_mutations(self, outgroup_mutations, parent_mutations):
        '''
        determine what mutations were needed to get from the parent to the current virus
        :param parent:
        :return:
        '''
        self.mutations_from_parent = outgroup_mutations - parent_mutations

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
        ddGFileName = "Average_" + structure + "_trimer_repaired_1968.fxout"
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

    def calculate_ddg_outgroup(self, calculated_stabilities):
        '''
        calls appropriate functions to calculate  ddG using foldx for the specified structure and ddG gets assigned to self.ddG_outgroup
        '''
        if len(calculated_stabilities) > 2:
            print("more than 2 stability calculations given")
            print(calculated_stabilities)
            raise Exception
        print("Calculated Stabilities: " + str(calculated_stabilities))
        self.ddg_outgroup['1HA0'] = float(calculated_stabilities[0])
        print("Outgroup ddg 1HA0: " + str(self.ddg_outgroup['1HA0']))
        self.ddg_outgroup['2YP7'] = float(calculated_stabilities[1])
        print("Outgroup ddg 2YP7: " + str(self.ddg_outgroup['2YP7']))

    def check_valid_structure(self, structure):
        '''
        :raises Exception if the structure given is not in the allowed structures for this pipeline
        '''
        if structure not in self.pdb_structures:
            raise Exception("This pipeline does not work for that structure, only works for 1HA0 or 2YP7")

    def calculate_ddg_parent(self, virus_ddg, parent_ddg):
        '''
        calculate the change in stability from the parent to the current virus
        '''
        self.ddg_parent['1HA0'] = float(virus_ddg[0]) - float(parent_ddg[0])
        self.ddg_parent['2YP7'] = float(virus_ddg[1]) - float(parent_ddg[1])
        '''
        for pdb in self.pdb_structures:
            try:
                self.ddg_parent[pdb] = self.ddg_outgroup[pdb] - parent.ddg_outgroup[pdb]
                print("ddg from parent for " + pdb + ": " + str(self.ddg_parent[pdb]))
            except:
                print("could not calculate ddg from parent for this pair")
                print(str(self.strain) + " | " + str(parent.strain))
                print(self.ddg_outgroup[pdb])
                print(parent.ddg_outgroup[pdb])
                raise Exception("Could not calculate ddG from parent")
        '''
    def determine_relative_time(self):
        '''
        Need to determine number of mutations from 1968 sequence as a proxy for time
        '''
        mutations_set = set()
        sequence_1968 = "MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"
        old_align_seq = sequence_1968[24:]
        virus_align_seq = self.seq[24:]
        for index in range(len(old_align_seq)):
            site = index + 9  # for both 1HA0 and 2YP7, start at site number 9 in structure ("STAT...")
            if old_align_seq[index] != virus_align_seq[index]:
                mutation = old_align_seq[index] + str(site) + virus_align_seq[index]
                mutations_set.add(mutation)
        self.mutations_from_1968 = mutations_set
        self.mutations_time = len(mutations_set)
        return ','.join(mutations_set)


'''
    def check_length_sequence(self):
        # check that the length sequence is always the same. Since by the time it comes out of tree_mutations,
        # should always be the same
'''
