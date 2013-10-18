#!/usr/bin/env	python2.7

from tiedie_util import *
import random
from collections import defaultdict
import copy

class ConcensusNetwork:

	def __init__(self, base_network, diffuser):

		# store nodes and edge counts
		self.nodes = defaultdict(int)
		self.edges = defaultdict(int)
		self.node_heats = defaultdict(list)
		# number of concensus networks
		self.data_points = 0

		self.base_network = base_network
		self.diffuser = diffuser

	def generate(self, input_set1, input_set2, rounds, sample_rate, options):

		# edges/nodes will be counted for each data subsample, and each 
		# size-cutoff option within those subsampled networks
		self.num_rounds = rounds*len(options['size'])
		subsample_set = options['subsample_which']

		for i in range(0, rounds):

			# generate random subsamples of the input sets
			subset1 = {}
			subset2 = {}

			# permute the upstream set unless 'd' is specified
			if subsample_set != "d":
				subset1 = {}
				for key in random.sample(input_set1.keys(), int(sample_rate*len(input_set1))):
					subset1[key] = input_set1[key]
			else:
				subset1 = input_set1

			# permute the downstream set unless 'u' is specified
			if subsample_set != "u":
				subset2 = {}
				for key in random.sample(input_set2.keys(), int(sample_rate*len(input_set2))):
					subset2[key] = input_set2[key]
			else:
				subset2 = input_set2

				
			# diffuse subsamples
			diff_subset1 = self.diffuser.diffuse(subset1, reverse=False)
			diff_subset2 = self.diffuser.diffuse(subset2, reverse=True)

			size_ranges = options['size']
			# copy to modify size options
			tiedie_opts = copy.copy(options)
			tiedie_opts['size'] = None
	
			first_size = True
			for network_size in size_ranges:
				# extract network at this size cutoff
				# FIXME: this is quite inefficient, as the algorithm will repeat the same steps to find a 
				# proper heat-cutoff in each iteration. May be worth some re-engineering in the future. 
				tiedie_opts['size'] = network_size
				subnet_soln, subnet_soln_nodes, alpha_score, linker_scores = \
					extractSubnetwork(subset1, subset2, diff_subset1, diff_subset2, self.base_network, tiedie_opts)

				# on the first round only: store node heats to generate distribution over the outer loop
				if first_size:
					first_size = False
					for (node, heat) in linker_scores.items():
						self.node_heats[node].append(heat)
	
				# count each additional edge, and node	
				for s in subnet_soln:
					for (i, t) in subnet_soln[s]:
						self.edges[ (s, i, t) ] += 1

				for n in subnet_soln_nodes:
					self.nodes[n] += 1	

	def getStats(self):
		"""
			Return frequencies for edges, nodes, and 
			also a heat distribution over subsampled inputs for each node in the network
			(a dictionary of lists)
		"""	
		edge_fractions = defaultdict(float)
		node_fractions = defaultdict(float)

		for (edge, count) in self.edges.items():
			edge_fractions[edge] = count / float(self.num_rounds)

		for (node, count) in self.nodes.items():
			node_fractions[node] = count / float(self.num_rounds)

		return (edge_fractions, node_fractions, self.node_heats)

