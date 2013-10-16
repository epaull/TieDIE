#!/usr/bin/env  python2.7

import networkx as nx

# TODO perhaps Diffuser could be abstracted to an interface for better code
# clarity? - ESR
# TODO should Diffusers take an nx DiGraph as initialization arguments rather
# than create them? Storing graphss as hashes makes the code less clear IMO -ESR
class PPrDiffuser:

	def __init__(self, network):
		"""
		PPrDiffuser: object to perform the Personalized PageRank Algorithm
		This method creates the diffuser object from an networkx DiGraph() object, 
		which can then be used to diffuse vectors over this
		
		Input:
			- network : a network hash object
		"""

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
		"""
		Personal_Page_Rank: Get the personal pagerank of the supplied input vector

		Input: 
			- p_vector: A hash-map of input values for a selection (or all) nodes
			(if supplied nodes aren't in the graph, they will be ignored)

        Output:
			- A vector of diffused heats in hash-map (key,value) format
		"""
		input_pvec = None
		# without initializing this vector the initial probabilities will be flat
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
		"""
		Diffuse: perform generalized diffusion from the supplied input vector

		Input: 
			- p_vector: A dict of initial heat values for nodes. Keys are nodes'
              ID strings. Any node not given in the dictionary will receive no
              initial heat. Any supplied nodes not found in the graph will be
              ignored.

        Optional args:
            - reverse: reverse the direction of edges in the graph before
              applying heat. default=False

        Returns:
			- A vector of diffused heats in dict format where keys are the
              nodes' ID strings
    	"""
		return self.personal_page_rank(p_vector, reverse)
