


class mutation_stability(object):
    '''
    check  mutation and format it so that it's compatible with foldx structure 1HA0 and 2YP7
    '''
    def __init__(self, mut, structure):
        self.mut = mut  # list of mutations
        self.mut_set = set(mut)  # set of mutations
        self.mut_chain_info_set = set()
        self.structure = structure  # either 1HA0 or 2YP7

        if self.structure not in ["1HA0", "2YP7"]:
            raise ValueError("This program only works for pdb structures 1HA0 or 2YP7")


    def __str__(self):
        return ", ".join(self.mut_chain_info_set)

    def site_range_valid(self, mutation):
       '''
       protein structures (1HA0, 2YP7) are missing certain amino acid sites, method checks that mutation is in structure
       :param mutation: mutation in standard format
       :return: true if site is in structure, false if site range is not in structure
       '''
       lowerRange = 9
       upperRange = 502
       missing_lower = 328
       missing_upper = 333
       site = int(mutation[1:len(mutation) - 1])
       if missing_lower <= site <= missing_upper:  # in missing middle section
            return False
       elif lowerRange <= site <= upperRange: # in range of protein structure besides middle section
            return True
       else:
            return False

    def include_chain_info(self, mutation):
        '''
        includes chain information for each mutation passed to function. HA is a trimer so need to specify chain for
        foldx
        :param mutation: mutation in standard format
        '''

        set_with_chain_mutations = set()
        if self.structure == "1HA0":
            chains = ["A", "M", "Y"]
        elif self.structure == "2YP7":
            chains = ["A", "P", "E"]
        site = mutation[1:len(mutation) - 1]
        aa1 = mutation[0]
        aa2 = mutation[len(mutation)-1]
        for chain in chains:
            self.mut_chain_info_set.add(aa1+chain+site+aa2)

    def check_valid_mutation(self):
        '''
        checks each mutation in mut_set that it is a valid mutation for the structures 1HA0, 2YP7. Calls
        include_chain_info, which adds each mutation with chain info to self.mut_chain_info_set.
        '''
        for mutation in self.mut_set:
            site_valid = self.site_range_valid(mutation)
            if site_valid:
                self.include_chain_info(mutation)

    def get_formatted_mutations(self):
        self.check_valid_mutation()
        return ';'.join(self.mut_chain_info_set)
