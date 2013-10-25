#!/usr/bin/env	python2.7

from tiedie_util import *
from linkers import *
import random
from collections import defaultdict
import copy

class ConsensusNetwork:
    # TODO write docstring

	def __init__(self, base_network, diffuser):
        # TODO write docstring

		# store nodes and edge counts
		self.nodes = defaultdict(int)
		self.edges = defaultdict(int)
		self.node_heats = defaultdict(list)
		# number of concensus networks
		self.data_points = 0

		self.base_network = base_network
		self.diffuser = diffuser

	def generate(self, input_sets, rounds, sample_rate, options):
        # TODO write docstring
			
		size_ranges = options['size']
		# copy to modify size options
		tiedie_opts = copy.copy(options)
		tiedie_opts['size'] = None

		# edges/nodes will be counted for each data subsample, and each 
		# size-cutoff option within those subsampled networks
		self.num_rounds = rounds*len(options['size'])
		subsample_set = options['subsample_which']

		for i in range(0, rounds):

			# generate random subsamples of the input sets
			subsampled_inputs = {}
			subsampled_diffused = {}

			for input in input_sets:

				input_set = input_sets[input]
				subsampled_inputs[input] = {}
			
				# get a random sample of the inputs, copy to new dictionary
				RATE = sample_rate
				# boundary case: obviously we can't subsample a single input, but it's still a valid input...
				if int(len(input_sets)*RATE) < 1:
					RATE = 1
				for key in random.sample(input_set.keys(), int(RATE*len(input_set))):
					subsampled_inputs[input][key] = input_set[key]

				# diffuse subsamples
				subsampled_diffused[input] = self.diffuser.diffuse(subsampled_inputs[input], reverse=False)


			first_size = True
			for network_size in size_ranges:
				# FIXME: this is quite inefficient, as the algorithm will repeat the same steps to find a 
				# proper heat-cutoff in each iteration. May be worth some re-engineering in the future. 
				#try:
					# extract network at this size cutoff
				subnet_soln, subnet_soln_nodes, linker_scores, cutoff = \
					extractSubnetwork(self.base_network, subsampled_inputs, subsampled_diffused, network_size, tiedie_opts)
				#except Exception, err:
				#	# just penalize with zero counts if we can't find a subnetwork at all
				#	sys.stderr.write(Exception.str()+"\t"+str(err)+'\n')
				#	continue


				# on the first round only: store node heats to generate distribution over the outer loop
				if first_size:
					first_size = False
					for (node, heat) in linker_scores.items():
						self.node_heats[node].append(heat)
	
				for s in subnet_soln:
					for (i, t) in subnet_soln[s]:
						self.edges[ (s, i, t) ] += 1

				for n in subnet_soln_nodes:
					self.nodes[n] += 1	

	def getStats(self):
		"""
		Return statistics pertaining to the consensus network.
        
        Returns: a tuple containing 3 fields:
        0. a dictionary of edge frequencies where edge_freqs[('from_node',
           'interaction', 'to_node'] is the frequency with which that edge
           appeared in the subsample diffusions
        1. a dictionary of node frequencies where node_freqs['node'] is the
           fraction of the time in which that node appeared in the resulting
           networks from subsample diffusions
		2. a heat distribution over subsampled inputs for each node in the
           network as a dictionary of lists where heat['node'][i] is the
           final heat at 'node' in subsample i
		"""	
		edge_fractions = defaultdict(float)
		node_fractions = defaultdict(float)

		for (edge, count) in self.edges.items():
			edge_fractions[edge] = count / float(self.num_rounds)

		for (node, count) in self.nodes.items():
			node_fractions[node] = count / float(self.num_rounds)

		return (edge_fractions, node_fractions, self.node_heats)

