#!/usr/bin/env	python2.7

from tiedie_util import *
from linkers import *
import random
from collections import defaultdict
import copy
#from sklearn.neighbors import NearestNeighbors
from scipy import stats

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
				try:
					# extract network at this size cutoff
					subnet_soln, subnet_soln_nodes, linker_scores, cutoff = \
						extractSubnetwork(self.base_network, subsampled_inputs, subsampled_diffused, network_size, tiedie_opts)
				except Exception, err:
					# just penalize with zero counts if we can't find a subnetwork at all
					continue


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

class ConsensusCentroid:
	"""
	Instantiate using a set of tiedie vector scores. Compute the mean
	vector and find the correlation between each vector to that mean. 

		Modes are:
			- centroid: generate a single centroid and compute distances to
			that centroid. 
			- knn: compute the distance only to the k-nearest neighbors

		Inputs:
			vectors: a hash indexed by sample_id; each contains a hash of gene by score
			values
			

	"""

	def __init__(self, vectors, mode='knn', nbrs=5):

		self.mode = mode
		if self.mode == 'centroid':
			self.nbrs = None
		elif self.mode == 'knn':
			self.num_nbrs = int(nbrs)
		else:
			raise Exception("Error: unsupported mode!")

		# values: indexed by sample, then by gene
		self.vectors = vectors
		# Each sample's correlation to it's nearest centroid
		self.scores = {}
		# index the centroid vectors for each sample
		self.centroids = {}

		# genes must be in the same order: verify this by first creating a 
		# single index
		self.gene_idx = []
		for sample in self.vectors:
			for gene in self.vectors[sample]:
				self.gene_idx.append(gene)
			break

		if self.mode == 'centroid':
			# get a single centroid using all samples
			centroid = self.getCentroid(self.vectors.keys())
			self.scores = self.scoreSamples(self.vectors.keys(), centroid)
			self.centroids = {}
			for sample in self.vectors.keys():
				self.centroids[sample] = centroid
			return 

		# otherwise score based on KNN
		self.knnScore()

	def getCentroids(self):
		return self.centroids

	def getScores(self):
		return self.scores

	def knnScoreExternal(self, gene_scores_by_sample):
		"""
		Score external sets against the tiedie vectors stored here
		WARNING: Must call knnScore() before this method is run!

		Input:
			gene_vector: a gene-indexed vector
		"""
		scores = {}
		for sample in gene_scores_by_sample: 

			gene_vector = gene_scores_by_sample[sample]	
	
			if self.mode == 'centroid':
				# compute the centroid and score this sample
				centroid = self.getCentroid(self.scores.keys())
				# score just this sample against this centroid
				spearmanRHO, pval = self.scoreSample(gene_vector, centroid)
				scores[sample] = spearmanRHO
				continue

			input_vector = []
			for gene in self.gene_idx:

				if gene not in gene_vector:
					raise Exception("Error: gene not found in input vector!")

				input_vector.append(float(gene_vector[gene]))
			distances, indices = self.neigh.kneighbors([input_vector])

			# get the k-nearest neighbors in our current dataset, for this outside sample
			knn_samples = []
			for j in indices[0]:
				knn_samples.append(self.sample_idx[j])

			# compute the centroid and score this sample
			centroid = self.getCentroid(knn_samples)
			# score just this sample against this centroid
			spearmanRHO, pval = self.scoreSample(gene_vector, centroid)
			scores[sample] = spearmanRHO

		return scores

	def knnScore(self):
		"""	
			Build self.scores by doing a KNN search for each sample, and 
			computing the closeness to each

		"""
		# we need to convert to list vectors first
		list_vectors = {}

		for sample in self.vectors:
			list_vectors[sample] = []
			for gene in self.gene_idx:
				list_vectors[sample].append(float(self.vectors[sample][gene]))

		# now create a data matrix: columns are sample vectors
		data = []
		self.sample_idx = list_vectors.keys()
		for sample in self.sample_idx:
			data.append(list_vectors[sample])

		self.neigh = NearestNeighbors(self.num_nbrs)
		self.neigh.fit(data)
		distances, indices = self.neigh.kneighbors(data)
		# indexed by sample, then gene
		self.neighbors_gene_distributions = {}

		# store the list of neighbors for each sample	
		self.neighbors = {}	
		for i in range(0, len(indices)):
			# for each sample, get the KNN, and compute it's distance to those			
			# nbrs is an array of the sample indexes
			target_sample = self.sample_idx[i]

			knn_samples = []
			for j in indices[i]:
				knn_samples.append(self.sample_idx[j])

			# compute the centroid and score this sample
			centroid = self.getCentroid(knn_samples)
			# score just this sample against this centroid
			target_scores = self.scoreSamples([target_sample], centroid)
			self.scores[target_sample] = target_scores[target_sample]
			self.centroids[target_sample] = centroid

			# FIXME: store a index of 
			self.neighbors[target_sample] = knn_samples
			# indexed by sample, then gene
			#self.neighbors_gene_distributions[target_sample] = {}
			#for gene in self.vectors[target_sample]:
			#	self.neighbors_gene_distributions[target_sample][gene] = []
			#	self.neighbors_gene_distributions[target_sample][gene].append(self.vectors[target_sample][gene])
			#	for s in knn_samples:	
			#		self.neighbors_gene_distributions[target_sample][gene].append(self.vectors[s][gene])

	def getNeighbors(self, target_sample):
		return self.neighbors[target_sample]

	def getGeneDistributions(self):
		return self.neighbors_gene_distributions
			
	@staticmethod
	def vectorCorr(hash1, hash2, method='spearman'):
		"""
			Input:
			   2 hashes with the same set of keys indexing floating point values
			Output:
			   The spearman correlation and the significance of that correlation
		"""

		if method != 'spearman' and method != 'pearson':
			raise Exception("Error: only Spearman and Pearson correlation supported right now!")
	
		vals1 = []
		vals2 = []
	
		for node in hash1:
	
			if node not in hash2:
				continue
	
			vals1.append(float(hash1[node]))
			vals2.append(float(hash2[node]))

		corr = None
		pval = None
		if method == 'spearman':	
			corr, pval = stats.spearmanr(vals1, vals2)	
		elif method == 'pearson':	
			corr, pval = stats.pearsonr(vals1, vals2)	
			
		return (corr, pval)
	
	def scoreSample(self, sample_vector, centroid, method='spearman', sig_filter=None):
		if method != 'spearman':
			raise Exception("Error: only spearman correlation supported right now!")

		spearmanRHO, pval = ConsensusCentroid.vectorCorr(centroid, sample_vector, "pearson")
		#if sig_filter is not None and pval > float(sig_filter):

		return (spearmanRHO, pval)	


	def scoreSamples(self, sample_ids, centroid, method='spearman', sig_filter=None):
	
		if method != 'spearman':
			raise Exception("Error: only spearman correlation supported right now!")
	
		scores = {}	
		for sample_id in sample_ids:
			spearmanRHO, pval = self.scoreSample(self.vectors[sample_id], centroid, method)
			if sig_filter and pval > float(sig_filter):
				continue
			scores[sample_id] = spearmanRHO
	
		return scores	

	def getCentroid(self, sample_ids):
		"""
		Compute the centroid for the sample ids listed here
		"""

		# indexed by gene: get the average over all sample-specific networks
		all_scores = defaultdict(list)
		for sample_id in sample_ids:

			if sample_id not in self.vectors:
				raise Exception("Error: sample id: "+sample_id+" not in internal list")

			for gene in self.vectors[sample_id]:
				all_scores[gene].append(float(self.vectors[sample_id][gene]))

		# average the scores
		centroid_scores = {}
		for gene in all_scores:
			centroid_scores[gene] = Mean(all_scores[gene])
	
		return centroid_scores
		

