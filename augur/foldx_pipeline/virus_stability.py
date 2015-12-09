from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna





class virus_stability:


    def __init__(self, hash_code, accension, trunk, tip, seq):

        self.accession = accension
        self.hash_code = hash_code
        self.seq = seq
        self.trunk = trunk
        self.tip = tip

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

        print("accession: " + self.accession)
        print("hash: " + self.hash_code)
        print("seq: " + self.seq)
        print("trunk: " + self.trunk)
        print("tip: " + self.tip)

    def align_to_outgroup(self):
        # align outgroup and self.seq, return mutations from outgroup -> seq
        print("Aligning the protein structure sequence to the root sequence")
        mutations_set = set()
        for root_index in range(len(root_sequence)):
            structure_index = root_index
            site = root_index + 9  # for both 1HA0 and 2YP7
            '''
            if pdb_name.startswith("4WE4"):
                site = root_index + 9
            elif pdb_name.startswith("2HMG"):
                site = root_index + 1
            '''
            '''
                if root_index > 328:
                    structure
            '''
            #print(root_sequence[index] + " : " + structure_sequence[index])
            if root_sequence[root_index] != structure_sequence[structure_index]:
                mutation = structure_sequence[structure_index] + str(site) + root_sequence[root_index]
                list_mutations.append(mutation)
        report_mutations = ','.join(list_mutations)
        return report_mutations

    def calculate_ddg_outgroup(self):
        # read the ddG from the output file of most recent foldx run

    def calculate_ddg_parent(self, parent):
        self.ddg_parent= self.ddg_outgroup - parent.ddg_outgroup

    def check_same_virus(self, other):


