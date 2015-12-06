





class virus_stability:


    def __init__(self, hash_code, accension, seq, trunk, tip):

        self.accession = accension
        self.hash_code = hash_code
        self.seq = seq
        self.trunk = trunk
        self.tip = tip

        self.ddg_outgroup = None
        self.ddg_parent = None

        #self.outgroup_seq =

    def __str__(self):

        print("accession: " + self.accession)
        print("hash: " + self.hash)
        print("seq: " + self.seq)
        print("trunk: " + self.trunk)
        print("tip: " + self.tip)

    def align_to_outgroup(self):
        # align outgroup and self.seq, return mutations from outgroup -> seq

    def calculate_ddg_outgroup(self):
        # read the ddG from the output file of most recent foldx run

    def calculate_ddg_parent(self, parent):
        self.ddg_parent= self.ddg_outgroup - parent.ddg_outgroup

    def check_same_virus(self, other):


