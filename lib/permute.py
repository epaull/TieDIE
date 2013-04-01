
from copy import copy
from random import shuffle

class NetBalancedPermuter:
	"""
		Encapsulates the permutation logic for an input heat set. Permutes
		Node scores with other nodes of similar network degree. 
	"""


	def __init__(self, network, up_set):
		"""
			Input:
				network: net[source] = [(i, t)
				up_set: up_set[node] = score
		"""
	

		self.degrees = {}

		self.nodes = up_set.keys()
		# heuristic: block needs to be significantly larger than the input set size
		self.block_size = len(self.nodes)*10
		self.scores = {}
		for node in self.nodes:
			self.scores[(node, str(up_set[node]))] = 1

		# convert internal representation, and get degrees of each node
		for source in network:

			if source not in self.degrees:
				self.degrees[source] = 0

			for (i, target) in network[source]:	
				self.degrees[source] += 1

				if target not in self.degrees:
					self.degrees[target] = 0
				# add a degree for the incoming edge
				self.degrees[target] += 1	

		self.sorted_degrees = sorted(self.degrees.items(), key=lambda x:x[1], reverse=True)

			
	def permuteBlock(self, block):

		orig = copy(block)
		b = copy(block)
		map = {}
		shuffle(b)
		for i in range(0, len(b)):
			map[orig[i]] = b[i]

		return map

	def permuteOne(self):
	
		group_count = 0
		permuted_scores = {}
		block = []
		for (node, degree) in self.sorted_degrees:
			block.append(node)	
			group_count += 1
			if group_count % self.block_size == 0:
				map = self.permuteBlock(block)
				for (node, score) in self.scores:
					if node in map:
						permuted_scores[map[node]] = float(score)
				block = []	

		return permuted_scores	

	def permute(self, iterations):
		
		permuted = []	
		for i in range(0, iterations):
			permuted.append(self.permuteOne())

		return permuted
