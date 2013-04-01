
from permute import NetBalancedPermuter

class SizeSelector:
	"""
		Encapsulates the permutation logic for an input heat set. Permutes
		Node scores with other nodes of similar network degree. 
	"""


	def __init__(self, network, up_set, down_set, down_diffused, diffuser, no_permutations):
		"""
			Input:
				network: net[source] = [(i, t)
				up_set: up_set[node] = score
		"""
		self.network = network
		self.up_set = up_set
		self.down_set = down_set
		self.down_diffused = down_diffused
		self.diffuser = diffuser
		self.no_permutations = no_permutations
			
					
	def geometricSearch(self):
		
		# create permuted sets
		self.permuter = NetBalancedPermuter(self.network, self.up_set)
		self.permuted_heats = self.permuter.permute(self.no_permutations)

		real_diffused = self.diffuser.diffuse(self.up_set)	
		# for each size, sweep the set of permutations and create a score
		# self.background = {'size':[score1, score2...]}

		# diffuse each...
		diffused_heats = []
		for permuted_heat in self.permuted_heats:
			diffused_heats.append((permuted_heat, self.diffuser.diffuse(permuted_heat)))

		self.background = {}

		current_best_pval = 1

		starting = 0.0
		search_incr = 0.1

		self.graph = {}

		while True:

			size = starting + search_incr

			self.background[size] = {'real':None, 'permuted':[]}
			# get the real score first
			real_cutoff, real_score = findLinkerCutoff(self.up_set, self.down_set, real_diffused, self.down_diffused, size)
			self.background[size]['real'] = real_score

			for (permuted_heat, diffused_heat) in diffused_heats:
				permuted_cutoff, permuted_score = findLinkerCutoff(permuted_heat, self.down_set, diffused_heat, self.down_diffused, size)
				self.background[size]['permuted'].append(permuted_score)

			pval = self.computePval(size)
			self.graph[size] = (real_score, pval)

			if pval < current_best_pval:
				current_best_pval = pval
				# double the search interval
				search_incr = search_incr*2
			else:
				# pval greater: find something in this	
				starting = starting + search_incr/2
				search_incr = 0.1
			

	def generateBackground(self, size_options):

		# create permuted sets
		self.permuter = NetBalancedPermuter(self.network, self.up_set)
		self.permuted_heats = self.permuter.permute(self.no_permutations)

		# 
		self.size_opts = size_options

		# store size/pval info
		self.graph = {}
	
		real_diffused = self.diffuser.diffuse(self.up_set)	
		# for each size, sweep the set of permutations and create a score
		# self.background = {'size':[score1, score2...]}

		# diffuse each...
		diffused_heats = []
		for permuted_heat in self.permuted_heats:
			diffused_heats.append((permuted_heat, self.diffuser.diffuse(permuted_heat)))

		self.background = {}
		for size in size_options:

			self.background[size] = {'real':None, 'permuted':[]}
			# get the real score first
			# real_cutoff = alpha linker cutoff 
			real_cutoff, real_score = findLinkerCutoff(self.up_set, self.down_set, real_diffused, self.down_diffused, size)
			self.background[size]['real'] = real_score

			for (permuted_heat, diffused_heat) in diffused_heats:
				permuted_cutoff, permuted_score = findLinkerCutoff(permuted_heat, self.down_set, diffused_heat, self.down_diffused, size)
				self.background[size]['permuted'].append(permuted_score)

			pval = self.computePval(size)
			self.graph[size] = (real_score, pval)

			print str(size)+"\t"+str(pval)+"\t"+str(real_score)+"\t"+str(real_cutoff)
			sys.stdout.flush()	


	def computePval(self, size):

		real = self.background[size]['real']
		gt = 0
		for permuted_score in sorted(self.background[size]['permuted'], reverse=True):
			if permuted_score > real:
				gt += 1	
			else:
				break
		pval = float(gt+1)/float(len(self.background[size]['permuted'])+1)

		return pval

	def getGraph(self):

		return self.graph
