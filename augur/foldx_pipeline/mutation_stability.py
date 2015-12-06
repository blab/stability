import collections

class mutation_stability:
    '''
    check  mutation and format it so that it's compatible with foldx structure 1HA0 and 2YP7
    '''
    def __init__(self, mut):
        self.mut = mut  # list of mutations
        self.mut_set = set(mut)  # set of mutations
    def __str__(self):
        return (str(self.mut_set))
        #return ", ".join(self.mut_set)

    def site_range_valid(self, mutation):
       '''
       :task: protein structures (1HA0, 2YP7) are missing certain amino acid sites, method checks that mutation is in structure
       :param mutation: mutation in standard foldx format
       :return: true if site is in structure, false if site range is not in structure
       '''
       # compiled from 2YP7 and 1HA0 to run the same residues
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

    def check_valid_mutation(self):
        invalid_mutations = set()
        for mutation in self.mut_set:
            site_valid = self.site_range_valid(mutation)
            if not site_valid:
                invalid_mutations.add(mutation)
                #self.mut_set.remove(mutation)
        self.mut_set.difference_update(invalid_mutations)