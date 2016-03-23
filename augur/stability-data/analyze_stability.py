import matplotlib.pyplot as plt
import json
import itertools
from scipy.stats.stats import pearsonr

class analyze_stability(object):


    def __init__(self):
        input_file = open("ddg_output.txt", 'r')
        self.hash_to_virus = json.load(input_file)
        self.structures = ['2YP2', '2YP7', '1HA0']

    def scatter_ddg(self, attribute):
        ddg = {}
        for structure in self.structures:
            ddg[structure] = []
        for virus in self.hash_to_virus.values():
            for structure in self.structures:
                ddg[structure].append(virus[structure][attribute])
        combinations = list(itertools.combinations(range(len(self.structures)), 2))
        plt.figure(1)
        subplot = 311
        for pair in combinations:
            plt.subplot(subplot)
            subplot += 1
            x = ddg[self.structures[pair[0]]]
            y = ddg[self.structures[pair[1]]]
            plt.scatter(x,y)
            plt.xlabel(self.structures[pair[0]])
            plt.ylabel(self.structures[pair[1]])
        plt.show()

if __name__ == "__main__":
    run = analyze_stability()
    run.scatter_ddg('sum_ddg')