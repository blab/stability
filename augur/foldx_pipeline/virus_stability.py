from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mutation_stability import mutation_stability




class virus_stability:


    def __init__(self, hash_code, accension, trunk, tip, seq):

        self.accession = accension
        self.hash_code = hash_code
        self.trunk = trunk
        self.tip = tip
        self.seq = seq
        self.mutations_from_outgroup = set()


        self.ddg_outgroup = None
        self.ddg_parent = None

        self.structure_1HA0_seq = "STATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTQGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRNEALNNRFQI"
        self.structure_2YP7_seq = "STATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRKSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTQGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQI"

        # get outgroup amino_acid sequence
        try:
            temp_outgroup = SeqIO.read('source-data/H3N2_outgroup.gb', 'genbank')
        except:
            print("could not find source-data/H3N2_outgroup.gb which contained the outgroup file")
        dna_seq = str(temp_outgroup.seq).upper()
        coding_dna = Seq(dna_seq, generic_dna)
        protein = coding_dna.translate()
        self.outgroup_seq = str(protein)

    def __str__(self):
        return ("|hash: %s |accession: %s |trunk: %s |tip: %s\n\t|sequence: %s") %(self.hash_code, self.accession, self.trunk, self.tip, self.seq)

    def align_to_outgroup(self):
        '''
        :mutates: self.mutations_from_outgroup
        aligns outgroup and self.seq, returns mutations from outgroup -> self.seq. Modifies self.mutations_from_outgroup
        :return: comma-seperated string of mutations from outgroup
        '''
        #
        mutations_set = set()
        outgroup_align_seq = self.outgroup_seq[24:]
        for outgroup_index in range(len(self.outgroup_seq)):
            seq_index = outgroup_index
            site = outgroup_index + 9  # for both 1HA0 and 2YP7
            if self.outgroup_seq[outgroup_index] != self.seq[seq_index]:
                mutation = self.outgroup_seq[outgroup_index] + str(site) + self.seq[seq_index]
                mutations_set.add(mutation)
        self.mutations_from_outgroup = mutations_set
        return ','.join(mutations_set)

'''
    def check_length_sequence(self):
        # check that the length sequence is always the same. Since by the time it comes out of tree_mutations,
        # should always be the same

    def calculate_ddg_outgroup(self):
        # read the ddG from the output file of most recent foldx run

    def calculate_ddg_parent(self, parent):
        self.ddg_parent= self.ddg_outgroup - parent.ddg_outgroup

    def check_same_virus(self, other):
        # just compare the sequences
'''
