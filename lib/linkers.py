from __future__ import print_function
from __future__ import division
#from __future__ import unicode_literals
from tiedie_util import *

def min(vals):
	min = vals[0]
	for v in vals:
		if v < min:
			min = v

	return min

def getProduct(diffused):
	gene_scores = {}
	for file in diffused:
		# a hash of hashes: file is the index
		for (gene, heat) in diffused[file].iteritems():
			if gene not in gene_scores:
				gene_scores[gene] = []
			gene_scores[gene].append(heat)

	gene_products = {}
	for gene in gene_scores:
		product = 1
		for v in gene_scores[gene]:
			product *= v
		gene_products[gene] = product

	return gene_products


def getMinHeats(consider_top, diffused):
	"""
	Gets the minimum heats for all genes, from a number of diffused heat vectors.

	Input:
		diffused = { 'set':{'gene1':heat1, 'gene2':...}

	Returns:
		A minimum-heat vector over all genes

	"""

	gene_scores = {}
	for file in diffused:
		# a hash of hashes: file is the index
		for (gene, heat) in diffused[file].iteritems():
			if gene not in gene_scores:
				gene_scores[gene] = []
			gene_scores[gene].append(heat)


	min_gene_values = {}
	for gene in gene_scores:
		values = gene_scores[gene]
		# get the top X
		min_gene_values[gene] = min(sorted(values, reverse=True)[0:consider_top])

	return min_gene_values

def getMaxHeats(consider_top, diffused):
	"""
	Gets the maximum heats for all genes, from a number of diffused heat vectors.

	Input:
		diffused = { 'set':{'gene1':heat1, 'gene2':...}

	Returns:
		A minimum-heat vector over all genes

	"""

	gene_scores = {}
	for file in diffused:
		# a hash of hashes: file is the index
		for (gene, heat) in diffused[file].iteritems():
			if gene not in gene_scores:
				gene_scores[gene] = []
			gene_scores[gene].append(heat)


	max_gene_values = {}
	for gene in gene_scores:
		values = gene_scores[gene]
		# get the top X
		max_gene_values[gene] = max(sorted(values, reverse=True)[0:consider_top])

	return max_gene_values

def extractSubnetwork(network, input_heats, diffused_heats, size_control, opts):
	"""
		Generate a spanning subnetwork from the supplied inputs, diffused heats and
		size control cutoff

		Input:
			- input heats
			- diffused input heats
			- size control factor

		Output:
			- spanning network
			- list of nodes in that network

	>>> input_heats = {'source':{'s1':0.01,'s2':0.5}, 'target':{'e1':0.49}}
	>>> diffused_heats = {'source':{'s1':0.05,'s2':0.4}, 'target':{'e1':0.4,'t2':0.3,'t1':0.1}}
	>>> s1 = set()
	>>> s1.add(('-t>','t1'))
	>>> s2 = set()
	>>> s2.add(('-t>','t2'))
	>>> t2 = set()
	>>> t2.add(('-t>','e1'))
	>>> network = {'s1':s1, 's2':s2, 't2':t2}
	>>> extractSubnetwork(network, input_heats, diffused_heats, 0.25, {})
	({'s2': set([('-t>', 't2')]), 't2': set([('-t>', 'e1')])}, set(['s2', 't2', 'e1']), {'s2': 0.4, 's1': 0.05, 't2': 0.3, 'e1': 0.4, 't1': 0.1}, 0.2999)

	"""

	linker_cutoff = None
	linker_nodes = None
	linker_scores = None

	# get linker heats as a function of input sets
	consider = len(input_heats)

	linker_heats = getMinHeats(consider, diffused_heats)
	#linker_heats = getProduct(diffused_heats)

	EPSILON = 0.0001
	linker_cutoff = None
	score = None
	linkers = set()
	linker_scores = {}
	for (l,h) in sorted(linker_heats.iteritems(), key=operator.itemgetter(1), reverse=True):
		c = h-EPSILON
		score, size_frac = linkerScore(input_heats, linker_heats, c, size_control)
		linker_cutoff = c
		linkers.add(l)
		linker_scores[l] = h
		if size_frac > 1:
			break

	input_genes = set()
	for input in input_heats:
		input_genes = input_genes.union(input_heats[input].keys())
	# set of input heats
	ugraph = None
	# USE LINKER GENES AND INPUT GENES
	active_nodes = set(linkers)
	active_nodes = active_nodes.union(input_genes)
	ugraph = connectedSubnets(network, active_nodes)
	if len(ugraph) == 0:
		sys.stderr.write("Couldn't find any linking graph at this size setting!\n")
	subnet_soln = mapUGraphToNetwork(ugraph, network)

	subnet_soln_nodes = set()
	for s in subnet_soln:
		subnet_soln_nodes.add(s)
		for (i,t) in subnet_soln[s]:
			subnet_soln_nodes.add(t)

	return (subnet_soln, subnet_soln_nodes, linker_heats, linker_cutoff)

def linkerScore(input_heats, min_heats, cutoff, size):
	"""
		Get linkers greater than this cutoff according to reverse-sorted list.
		This version takes an arbitrary number of inputs.

		Inputs:
			input_heats: a dictionary of an arbitrary number of input heat sets.
			min_heats: pre-processed 'linker' heat values according to any particular
			linker function.
	"""

	# get the set of all input genes
	all_inputs = set()
	for name in input_heats:
		all_inputs = all_inputs.union(input_heats[name].keys())

	# generate the set of linker genes according to the supplied heat cutoff.
	all_linkers = set()
	for (gene, heat) in sorted(min_heats.iteritems(), key=operator.itemgetter(1), reverse=True):
		if heat < cutoff:
			break
		all_linkers.add(gene)

	# generate the union of input and linker sets, the exclusive 'connecting' set of linker/non-input genes
	# and score based on the fractional criterion
	all_genes = all_inputs.union(all_linkers)
	connecting = all_linkers.difference(all_inputs)
	score = len(connecting)/float(len(all_linkers))
	# the relative size of the connecting genes, compared to the input set sizes
	size_frac = (len(connecting)/float(len(all_inputs)))/float(size)

	return (score, size_frac)
