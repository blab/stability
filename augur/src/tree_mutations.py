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
		for node in self.tree.postorder_node_iter():
			pp.pprint(vars(node))
