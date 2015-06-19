import os, re, time
import dendropy
from seq_util import *
from date_util import *

class tree_mutations(object):
	def __init__(self, **kwargs):
		'''
		'''
		
	def catalog_mutations(self):
		'''
		run through and print mutations

		'''                     
		
		import pprint
		pp = pprint.PrettyPrinter(indent=4)
		for node in self.tree.postorder_internal_node_iter():
			sequence = str(node.seq)
                        for child in node.child_nodes():
        					mutation = str(child.aa_muts)
        					trunk = str(child.trunk)
        					#pp.pprint(mutation + '\n' + trunk + '\n' + sequence)
        					pp.pprint(mutation)
        					pp.pprint(trunk)
        					pp.pprint(sequence)