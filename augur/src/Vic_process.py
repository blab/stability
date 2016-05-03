import matplotlib as mpl
mpl.use('pdf')
import time, re, os
from virus_filter import flu_filter, fix_name
from virus_clean import virus_clean
from tree_refine import tree_refine
from tree_titer import HI_tree
from H3N2_process import H3N2_refine as BVic_refine
from process import process, virus_config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

sp = 15
epitope_mask = np.array(['1' if pos in [141,142,145,146,172,176,178,179,180,181,183,184,185, #Sa
										170,173,174,177,206,207,210,211,212,214,216,		 #Sb
										183,187,191,196,221,225,254,258,288,				 #Ca1
										154,157,158,159,161,163,238,239,242,243,			 #Ca2
										87, 88, 90, 91, 92, 95, 96, 98, 99, 100, 132, 139	 #Cb
									   ]
						else '0' for pos in xrange(1,1725)])

receptor_binding_sites = [159,169,170,172,173,203,207]


virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'virus':'Vic',
	'alignment_file':'data/Vic_gisaid_epiflu_sequence.fasta',
	'outgroup':'B/HongKong/02/1993',
	'force_include':'data/Vic_HI_strains.txt',
	'force_include_all':False,
	'date_spec':'year',
	'max_global':True,   # sample as evenly as possible from different geographic regions
	# define relevant clades in canonical HA1 numbering (+1)
	# numbering starting at methionine including the signal peptide
	'clade_designations': {
		'1A': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K')],
		'1B': [('HA1', 75,'K'), ('HA1', 58, 'P'), ('HA1', 165, 'K')],
		'117V': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K'), ('HA1', 129, 'D'), ('HA1', 117, 'V')]
	},
	'HI_fname':'data/Vic_HI_titers.txt',
	'html_vars': {'coloring': 'lbi, dfreq, region, date, cHI, HI_dist',
				  'gtplaceholder': 'HA1 positions...',
				  'freqdefault': '1A, 1B'},
	'js_vars': {'LBItau': 0.0005, 'LBItime_window': 0.5, 'dfreq_dn':2},
	'layout':'auspice',
	})


class BVic_filter(flu_filter):
	def __init__(self,min_length = 987, **kwargs):
		'''
		parameters
		min_length  -- minimal length for a sequence to be acceptable
		'''
		flu_filter.__init__(self, **kwargs)
		self.min_length = min_length
		self.vaccine_strains =[
			{
				'strain':    	'B/Shangdong/7/97',
				'isolate_id':	'EPI_ISL_1790',
				'date':    		'1997-07-01', #(Month and day unknown)
				'region':   	'China',
				'seq':'GATCGAATCTGCACTGGGATAACATCGTCAAACTCACCCCATGTGGTCAAAACTGCTACTCAAGGGGAGGTCAATGTGACTGGTGTGATACCACTGACAACAACACCCACCAAATCTCATTTTGCAAATCTCAAAGGAACAAAAACCAGAGGGAAACTATGCCCAAAATGCCTCAACTGTACAGATCTGGACGTGGCCTTGGGCAGACCAAAATGCACGGGGAACATACCTTCGGCAAAAGTTTCAATACTCCATGAAGTCAGACCTGTTACATCTGGGTGCTTTCCTATAATGCACGACAGAACAAAAATTAGACAGCTGCCCAATCTTCTCAGAGGATACGAACATATCAGGTTATCAATTCATAACGTTATCAATGCAGAAAAGGCACCAGGAGGACCCTACAAAATTGGAACCTCAGGGTCTTGCCCTAACGTTACCAATGGAAACGGATTCTTCGCAACAATGGCTTGGGCCGTCCCAAAAAACGACAACAACAAAACAGCAACAAATTCATTAACAATAGAAGTACCATACATTTGTACAGAAGGAGAAGACCAAATTACCGTTTGGGGGTTCCACTCTGATAACGAAAACCAAATGGCAAAACTCTATGGGGACTCAAAGCCCCAGAAGTTCACCTCATCTGCCAACGGAGTGACCACACATTACGTTTCACAGATTGGTGGCTTCCCAAATCAAACAGAAGACGGAGGACTACCACAAAGTGGTAGAATTGTTGTTGATTACATGGTGCAAAAATCTGGGAAAACAGGAACAATTACCTATCAAAGAGGTATTTTATTGCCTCAAAAAGTGTGGTGCGCAAGTGGCAGGAGCAAGGTAATAAAAGGGTCCTTGCCTTTAATTGGAGAAGCAGATTGCCTCCACGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTATTACACAGGGGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTGAAAACACCCTTGAAGCTGGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCTTAGAAGGAGGATGGGAAGGAATGATTGCAGGTTGGCACGGATACACATCCCATGGAGCACATGGAGTAGCAGTGGCAGCAGACCTTAAGAGTACTCAAGAAGCCATAAACAAGATAACAAAAAATCTCAACTCTTTGAGTGAGCTGGAAGTAAAGAATCTTCAAAGACTAAGCGGTGCCATGGATGAACTCCACAACGAAATACTAGAACTAGACGAGAAAGTGGATGATCTCAGAGCTGATACAATAAGCTCGCAAATAGAACTCGCAGTCTTGCTTTCCAAT',
			},
			{
				'strain':   'B/HongKong/330/2001',
				'isolate_id': 'EPI_ISL_2342',
				'date':    	'2001-07-01', 	#(Month and day unknown)
				'region':	'China',
				'seq':   	'GATCGAATCTGCACTGGAATAACATCGTCAAACTCACCCCATGTGGTCAAAACTGCTACTCAAGGGGAAGTCAATGTGACTGGTGTGATACCACTGACAACAACACCCACCAAATCTCATTTTGCAAATCTCAAAGGAACAAAAACCAGAGGGAAACTATGCCCAAAATGTCTCAACTGCACAGATCTGGACGTGGCCTTGGGCAGACCAAAATGCACGGGGAACATACCTTCGGCAAAAGTTTCAATACTCCATGAAGTAAGACCTGTTACATCTGGGTGCTTTCCTATAATGCACGACAGAACAAAAATTAGACAGCTGCCCAATCTTCTCAGAGGATACGAACGTATCAGGTTATCAAACCATAACGTTATCAATGCAGAAAAAGCACCAGGAGGACCCTACAAAATTGGAACCTCAGGGTCTTGCCCTAACGTTACCAATGGAAACGGATTCTTCGCAACAATGGCTTGGGCTGTCCCAAAAAACGAAAACAACAAAACAGCAACAAATTCATTAACAATAGAAGTACCATACATTTGTACAGAAGGAGAAGACCAAATTACCGTTTGGGGGTTCCACTCTGATAGCGAAACCCAAATGGCAAAACTCTATGGAGACTCAAAGCCTCAGAAGTTCACTTCATCTGCCAACGGAGTGACCACACATTACGTTTCACAGATTGGTGGCTTCCCAAATCAAACAGAAGACGGAGGACTACCACAAAGTGGTAGAATTGTTGTTGATTACATGGTGCAAAAATCTGGAAAAACAGGAACAATTACCTATCAAAGAGGTATTTTATTGCCTCAAAAAGTGTGGTGCGCAAGTGGCAGGAGCAAGGTAATAAAAGGATCCTTGCCTTTAATTGGAGAAGCAGATTGCCTCCACGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTATACAGGGGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTGAAAACACCCTTGAAGCTGGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCTTCGGAGCTTTGGGCTG',
			},
			{
				'strain': 'B/Malaysia/2506/2004',
				'isolate_id': 'EPI_ISL_21142',
				'date':'2004-07-01', # (Month and day unknown) |   |
				'region':'SouthEast Asia',
				'seq':'ATTGTACTACTCATGGTAGTAACATCCAATGCAGATCGAATCTGCACTGGGATAACATCGTCAAACTCACCACATGTTGTCAAAACTGCTACTCAAGGGGAGGTCAATGTGACTGGTGTAATACCACTGACAACAACACCCACCAAATCTCATTTTGCAAATCTCAAAGGAACAGAAACCAGAGGGAAACTATGCCCAAAATGTCTCAACTGCACAGATCTGGACGTGGCCTTGGGCAGACCAAAATGCACGGGGAACATACCCTCGGCAAGAGTTTCAATACTCCATGAAGTCAGACCTGTTACATCTGGGTGCTTTCCTATAATGCACGACAGAACAAAAATTAGACAGCTGCCTAACCTTCTCAGAGGATACGAACATATCAGGTTATCAACTCATAACGTTATCAATGCAGAAAATGCACCAGGAGGATCCTACAAAATTGGAACCTCAGGGTCTTGCCCTAACGTTACCAATGGAAACGGATTTTTCGCAACAATGGCTTGGGCCGTCCCAAAAAACGACAACAACAAAACAGCAACAAATTCATTAACAATAGAAGTACCATACATTTGTACAGAAGGAGAAGACCAAATTACCGTTTGGGGGTTCCACTCTGATAACGAAGCCCAAATGGCAAAGCTCTATGGGGACTCAAAGCCCCAGAAGTTCACCTCATCTGCCAACGGAGTGACCACACATTACGTTTCACAGATTGGTGGCTTCCCAAATCAAACAGAAGACGGAGGACTACCACAAAGTGGTAGAATTGTTGTTGATTACATGGTGCAAAAATCTGGGAAAACAGGAACAATTACCTATCAAAGAGGTATTTTATTGCCTCAAAAAGTGTGGTGCGCAAGTGGCAGGAGCAAGGTAATAAAAGGATCCTTGCCTTTAATTGGAGAAGCAGATTGCCTCCACGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGGGAACATGCAAAGGCCATAGGAAATTGCCCAATATGGGTGAAAACACCCTTGAAGCTGGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCTTAGAAGGAGGATGGGAAGGAATGATTGCAGGTTGGCACGGATACACATCCCATGGGGCACATGGAGTAGCGGTGGCAGCAGACCTTAAGAGCACTCAAGAGGCCATAAACAAGATAACAAAAAATCTCAACTCTTTGAGTGAGCTGGAAGTAAAGAATCTTCAAAGACTAAGCGGTGCCATGGATGAACTCCACAACGAAATACTAGAACTAGACGAGAAAGTGGATGATCTCAGAGCTGATACAATAAGCTCACAAATAGAACTCGCAGTCCTGCTTTCCAATGAAGGAATAATAAACAGTGAAGATGAGCATCTCTTGGCGCTTGAAAGAAAGCTGAAGAAAATGCTGGGCCCCTCTGCTGTAGAGATAGGGAATGGATGCTTTGAAACCAAACACAAGTGCAACCAGACCTGTCTCGACAGAATAGCTGCTGGTACCTTTGATGCAGGAGAATTTTCTCTCCCCACTTTTGATTCACTGAATATTACTGCTGCATCTTTAAATGACGATGGATTGGATAATCATACTATACTGCTTTACTACTCAACTGCTGCCTCCAGTTTGGCTGTAACATTGATGATAGCTATCTTTGTTGTTTATATGGTCTCCAGAGACAATGTTTCTTGCTCCATCTGTCTATAAGGAAAGTTAAACCCTGTATTTTCCTTTATTGTAGTGCTTGTTTGCTTGTTACCATTACAAAAAACGGTTATTGAAAAATGCTCTTGTTACTACTAATA',
			},
			{
				'strain':'B/Brisbane/60/2008',
				'isolate_id':'EPI_ISL_24365',
				'date': '2008-08-04',
				'lab':'Queensland Health Scientific Services',
				'region':'Oceania',
				'seq':'AGCAGAAGCAGAGCATTTTCTAATATCCACAAAATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCAGATCGAATCTGCACTGGGATAACATCGTCAAACTCACCACATGTCGTCAAAACTGCTACTCAAGGGGAGGTCAATGTGACTGGTGTAATACCACTGACAACAACACCCACCAAATCTCATTTTGCAAATCTCAAAGGAACAGAAACCAGGGGGAAACTATGCCCAAAATGCCTCAACTGCACAGATCTGGACGTAGCCTTGGGCAGACCAAAATGCACGGGGAAAATACCCTCGGCAAGAGTTTCAATACTCCATGAAGTCAGACCTGTTACATCTGGGTGCTTTCCTATAATGCACGACAGAACAAAAATTAGACAGCTGCCTAACCTTCTCCGAGGATACGAACATATCAGGTTATCAACCCATAACGTTATCAATGCAGAAAATGCACCAGGAGGACCCTACAAAATTGGAACCTCAGGGTCTTGCCCTAACATTACCAATGGAAACGGATTTTTCGCAACAATGGCTTGGGCCGTCCCAAAAAACGACAAAAACAAAACAGCAACAAATCCATTAACAATAGAAGTACCATACATTTGTACAGAAGGAGAAGACCAAATTACCGTTTGGGGGTTCCACTCTGACAACGAGACCCAAATGGCAAAGCTCTATGGGGACTCAAAGCCCCAGAAGTTCACCTCATCTGCCAACGGAGTGACCACACATTACGTTTCACAGATTGGTGGCTTCCCAAATCAAACAGAAGACGGAGGACTACCACAAAGTGGTAGAATTGTTGTTGATTACATGGTGCAAAAATCTGGGAAAACAGGAACAATTACCTATCAAAGGGGTATTTTATTGCCTCAAAAGGTGTGGTGCGCAAGTGGCAGGAGCAAGGTAATAAAAGGATCCTTGCCTTTAATTGGAGAAGCAGATTGCCTCCACGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGGGAACATGCAAAGGCCATAGGAAATTGCCCAATATGGGTGAAAACACCCTTGAAGCTGGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCTTAGAAGGAGGATGGGAAGGAATGATTGCAGGTTGGCACGGATACACATCCCATGGGGCACATGGAGTAGCGGTGGCAGCAGACCTTAAGAGCACTCAAGAGGCCATAAACAAGATAACAAAAAATCTCAACTCTTTGAGTGAGCTGGAAGTAAAGAATCTTCAAAGACTAAGCGGTGCCATGGATGAACTCCACAACGAAATACTAGAACTAGATGAGAAAGTGGATGATCTCAGAGCTGATACAATAAGCTCACAAATAGAACTCGCAGTCCTGCTTTCCAATGAAGGAATAATAAACAGTGAAGATGAACATCTCTTGGCGCTTGAAAGAAAGCTGAAGAAAATGCTGGGCCCCTCTGCTGTAGAGATAGGGAATGGATGCTTTGAAACCAAACACAAGTGCAACCAGACCTGTCTCGACAGAATAGCTGCTGGTACCTTTGATGCAGGAGAATTTTCTCTCCCCACCTTTGATTCACTGAATATTACTGCTGCATCTTTAAATGACGATGGATTGGATAATCATACTATACTGCTTTACTACTCAACTGCTGCCTCCAGT',
			}
		]
		tmp_outgroup = SeqIO.read('source-data/Vic_outgroup.gb', 'genbank')
		genome_annotation = tmp_outgroup.features
		self.cds = {x.qualifiers['gene'][0]:x for x in genome_annotation
				if 'gene' in x.qualifiers and x.type=='CDS' and
				x.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}
		self.outgroup = {
						'strain':'B/HongKong/02/1993',
						'region':'China',
						'isolate_id':'EPI_ISL_6617',
						'date':'1993-02-15', #(Month and day unknown)
						'seq': str(tmp_outgroup.seq).upper()
						}


class BVic_clean(virus_clean):
	def __init__(self,**kwargs):
		virus_clean.__init__(self, **kwargs)

	def clean_outbreaks(self):
		"""Remove duplicate strains, where the geographic location, date of sampling and sequence are identical"""
		virus_hashes = set()
		new_viruses = []
		for v in self.viruses:
			try:
				geo = re.search(r'B/([^/]+)/', v.strain).group(1)
				if geo:
					vhash = (geo, v.date, str(v.seq))
					if vhash not in virus_hashes:
						new_viruses.append(v)
						virus_hashes.add(vhash)
			except:
				print "Error parsing geo info for", v.strain

		self.viruses = MultipleSeqAlignment(new_viruses)
		return new_viruses

	def clean_outlier_strains(self):
		"""Remove single outlying viruses"""
		remove_viruses = []
		outlier_strains = ["B/Bangkok/SI17/2012", "B/Bangkok/SI58/2012", "B/Kol/2024/2008"]
		for outlier_strain in outlier_strains:
			for v in self.viruses:
				if (v.strain == outlier_strain):
					remove_viruses.append(v)
					if self.verbose > 1:
						print "\tremoving", v.strain
		self.viruses = MultipleSeqAlignment([v for v in self.viruses if v not in remove_viruses])

	def clean(self):
		self.clean_generic()
		self.clean_outbreaks()
		print "Number of viruses after outbreak filtering:",len(self.viruses)
		self.clean_outlier_strains()
		print "Number of viruses after outlier filtering:",len(self.viruses)

class BVic_process(process, BVic_filter, BVic_clean, BVic_refine, HI_tree):
	"""docstring for BVic_process, BVic_filter"""
	def __init__(self,verbose = 0, force_include = None,
				force_include_all = False, max_global= True, **kwargs):
		self.force_include = force_include
		self.force_include_all = force_include_all
		self.max_global = max_global
		process.__init__(self, **kwargs)
		BVic_filter.__init__(self,**kwargs)
		BVic_clean.__init__(self,**kwargs)
		BVic_refine.__init__(self,**kwargs)
		HI_tree.__init__(self,**kwargs)
		self.verbose = verbose

	def run(self, steps, viruses_per_month=50, raxml_time_limit = 1.0, lam_HI=2.0, lam_pot=0.3, lam_avi=2.0):
		if 'filter' in steps:
			print "--- Virus filtering at " + time.strftime("%H:%M:%S") + " ---"
			self.filter()
			if self.force_include is not None and os.path.isfile(self.force_include):
				with open(self.force_include) as infile:
					forced_strains = [fix_name(line.strip().split('\t')[0]).upper() for line in infile]
			else:
				forced_strains = []
			self.subsample(viruses_per_month,
				prioritize=forced_strains, all_priority=self.force_include_all,
				region_specific = self.max_global)
			self.add_older_vaccine_viruses(dt = 6)
			self.dump()
		else:
			self.load()
		if 'align' in steps:
			self.align()   	# -> self.viruses is an alignment object
		if 'clean' in steps:
			print "--- Clean at " + time.strftime("%H:%M:%S") + " ---"
			self.clean()   # -> every node as a numerical date
			self.dump()
		if 'tree' in steps:
			print "--- Tree	 infer at " + time.strftime("%H:%M:%S") + " ---"
			self.infer_tree(raxml_time_limit)  # -> self has a tree
			self.dump()
		if 'ancestral' in steps:
			print "--- Infer ancestral sequences " + time.strftime("%H:%M:%S") + " ---"
			self.infer_ancestral()  # -> every node has a sequence
			self.dump()
		if 'refine' in steps:
			print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"
			self.refine()
			self.dump()
		if 'frequencies' in steps:
			print "--- Estimating frequencies at " + time.strftime("%H:%M:%S") + " ---"
			self.determine_variable_positions()
			self.estimate_frequencies(tasks = ["mutations", "tree"])
			if 'genotype_frequencies' in steps:
					self.estimate_frequencies(tasks = ["genotypes"])
			self.dump()
		if 'HI' in steps:
			print "--- Adding HI titers to the tree " + time.strftime("%H:%M:%S") + " ---"
			try:
				self.determine_variable_positions()
				self.map_HI(training_fraction=1.0, method = 'nnl1reg',
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=True)
				self.map_HI(training_fraction=1.0, method = 'nnl1reg', force_redo=True,
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=False)
				self.dump()
			except:
				print("HI modeling failed!")
		if 'export' in steps:
			self.add_titers()
			self.temporal_regional_statistics()
			# exporting to json, including the BVic specific fields
			self.export_to_auspice(tree_fields = [
				'ep', 'ne', 'rb', 'aa_muts','accession','isolate_id', 'lab','db', 'country',
				'dHI', 'cHI', 'mean_HI_titers','HI_titers','HI_titers_raw', 'serum', 'HI_info',
				'avidity_tree','avidity_mut', 'potency_mut', 'potency_tree', 'mean_potency_mut', 'mean_potency_tree', 'autologous_titers'],
				annotations = ['1A', '1B', '117V'])
			if params.html:
				self.generate_indexHTML()
			self.export_HI_mutation_effects()


		if 'HIvalidate' in steps:
			print "--- generating validation figures " + time.strftime("%H:%M:%S") + " ---"
			self.generate_validation_figures()


if __name__=="__main__":
	all_steps = ['filter', 'align', 'clean', 'tree', 'ancestral', 'refine', 'frequencies','HI', 'export'] + ['HIvalidate']
	from process import parser
	params = parser.parse_args()

	lt = time.localtime()
	num_date = round(lt.tm_year+(lt.tm_yday-1.0)/365.0,2)
	params.time_interval = (num_date-params.years_back, num_date)
	if params.interval is not None and len(params.interval)==2 and params.interval[0]<params.interval[1]:
		params.time_interval = (params.interval[0], params.interval[1])
	dt= params.time_interval[1]-params.time_interval[0]
	params.pivots_per_year = 12.0 if dt<5 else 6.0
	steps = all_steps[all_steps.index(params.start):(all_steps.index(params.stop)+1)]
	if params.skip is not None:
		for tmp_step in params.skip:
			if tmp_step in steps:
				print "skipping",tmp_step
				steps.remove(tmp_step)

	# add all arguments to virus_config (possibly overriding)
	virus_config.update(params.__dict__)
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myBVic = BVic_process(**virus_config)
	if params.test:
		myBVic.load()
	else:
		myBVic.run(steps,viruses_per_month = virus_config['viruses_per_month'],
			raxml_time_limit = virus_config['raxml_time_limit'],
				   lam_HI = virus_config['lam_HI'],
				   lam_avi = virus_config['lam_avi'],
				   lam_pot = virus_config['lam_pot'],
				   )
