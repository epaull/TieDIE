from __future__ import print_function
from __future__ import division
#from __future__ import unicode_literals
import networkx as nx

class PPrDiffuser:

	def __init__(self, network):
		'''
			PPrDiffuser: object to perform the Personalized PageRank Algorithm
			This method creates the diffuser object from an networkx DiGraph() object,
			which can then be used to diffuse vectors over this

			Input:
				- network : a network hash object
		'''

		# create a directed graph with NetworkX
		self.G = nx.DiGraph()
		# create a reversed graph for diffusion from the 'target' set
		self.G_reversed = nx.DiGraph()
		# convert network format to networkX Graph object
		for source in network:
			for (i,t) in network[source]:
				self.G.add_edge(source, t)
				self.G_reversed.add_edge(t, source)

	def personal_page_rank(self, p_vector, reverse=False):
		'''
			Personal_Page_Rank: Get the personal pagerank of the supplied input vector

			Input:
				- p_vector: A hash-map of input values for a selection (or all) nodes
				(if supplied nodes aren't in the graph, they will be ignored)

			Output:
				- A vector of diffused heats in hash-map (key,value) format
		'''
		input_pvec = None
		#  without initializing this vector the initial probabilities will be flat
		# and this will be equivalent to standard page rank
		if p_vector:
			input_pvec = {}
			# doesn't seem to be necessary for a non-zero epsilon now, but
			# leave this as a place holder
			epsilon = 0.0
			for node in self.G.nodes(data=False):
				if node in p_vector:
					input_pvec[node] = p_vector[node]
				else:
					input_pvec[node] = epsilon

		if reverse:
			return nx.pagerank_numpy(self.G_reversed, 0.85, input_pvec)
		else:
			return nx.pagerank_numpy(self.G, 0.85, input_pvec)

	def diffuse(self, p_vector, reverse=False):
		'''
			Diffuse: perform generalized diffusion from the supplied input vector

			Input:
				- p_vector: A hash-map of input values for a selection (or all) nodes
				(if supplied nodes aren't in the graph, they will be ignored)

			Output:
				- A vector of diffused heats in hash-map (key,value) format
		'''
		return self.personal_page_rank(p_vector, reverse)
