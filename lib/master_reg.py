
from tiedie_util import *
import operator
import math
import random
from scipy import stats
import numpy as np
from distributions import Dist
from itertools import combinations
import operator

class ActivityScores: 
	"""
		Uses the supplied pathway to find activity scores based on the expression of each
		candidate regulators downstream targets. This is similar to the 'Master Regulator'
		method developed by Andrea Califano's lab, but uses a simpler background model that
		swaps scores between genes rather than shuffling samples (in the case that you don't
		supply reference/null samples). The advantage of this over the Master Regulator method 
		is you don't need the original sample by expression matrix and you can use a signature 
		computed by a data-specific method such as DESeq, edgeR, etc. 
		The disadvantage is that you break the gene-gene correlations present
		in the dataset; I haven't done extensive testing to compare these approaches, but the 
		recommendation is to either provide normal or other reference samples, or use the official 
		Master Regulator method (see run-viper.R included here) as a default.

	"""

	def __init__(self, network, scores, min_hub=10, regulation_type="transcriptional", p=1):
		"""
			Input:
				network: net[source] = [(i, t)]
				scores: hash map of differential gene expression (think D-statistics from SAM)
				min_hub: minimum number of genes regulated transcriptionally required 
				to be considered as a potential 'master regulator'
				regulation_type: (post) transcriptional 
				p: the power to raise each element to when computing the running sum
		"""

		#FIXME: take in 
	 	considered_modes = set()
		if regulation_type == "transcriptional":
			considered_modes.add('t') 
		elif regulation_type == "post":
			considered_modes.add('a') 
		# build a list of candidate regulators
		self.candidates = {}	
		for source in network:

			positive_regulon = set()
			negative_regulon = set()
			for (i, t) in network[source]:
				type, mode = classifyInteraction(i)	
				# only consider transcriptional regulation
				if mode not in considered_modes:
					continue

				if type == 1:
					positive_regulon.add(t)
				elif type == -1:
					negative_regulon.add(t)

			if (len(positive_regulon) + len(negative_regulon)) >= int(min_hub):
				self.candidates[source] = (positive_regulon, negative_regulon)

		self.generateRankings(scores)
		# for chisquare test
		self.generateCategories(scores)

	@staticmethod
	def	getEnrichmentScore(network, scores, test_set, nperms=1000):
		mrObj = ActivityScores(network, scores, min_hub=10)

		# index all network nodes
		network_nodes = set()
		for s in network:
			network_nodes.add(s)
			for (i, t) in network[s]:
				network_nodes.add(t)

		# generate GSEA score
		score = mrObj.scoreReg(test_set, set())
		# perform random permutations, get background scores
		no_gte = 0.0
		for i in range(0, nperms):
			permuted_set = random.sample(network_nodes, len(test_set))
			p_score = mrObj.scoreReg(permuted_set, set())
			if p_score >= score:
				no_gte += 1.0

		pval = (no_gte+1)/(nperms+1)
		return (score, pval)

	@staticmethod
	def findRegulators(network, de_file, min_hub=10, nperms=1000):
		"""
		Input:
			file with differential expression (or otherwise scored) values 
		
		Returns:
			A hash of master regulators, with signed, weighted scores normalized
			so that absolute values sum to 1.
		"""
		scores, signs = parseHeats(de_file)
		mrObj = ActivityScores(network, scores, min_hub=min_hub)
		# perform 1000 random permutations of the data to get significance scores for each
		result = mrObj.scoreCandidates(nperms)
		tfs_heats = {}
		for (tf, result) in sorted(result.items(), key=lambda t: t[1][0]):
			# filter on p-value
			if result[1] > 0.05:
				continue
			tfs_heats[tf] = float(result[0])

		if len(tfs_heats) == 0:
			raise Exception("No Significant Regulators Active!")
	
		t_total = 0
		for (g, h) in tfs_heats.items():
			t_total += abs(float(h))

		# normalize abs values to sum to 1	
		norm_factor = 1.0/t_total
	
		for (g, h) in tfs_heats.items():
			tfs_heats[g] = h*norm_factor

		return tfs_heats

	@staticmethod
	def findRegulatorsReport(network, de_file, min_hub=10, nperms=1000, reg_type="transcriptional"):
		"""
		Input:
			file with differential expression (or otherwise scored) values 
		
		Returns:
			A hash of master regulators, with signed, weighted scores normalized
			so that absolute values sum to 1.
		"""
		scores, signs = parseHeats(de_file)
		mrObj = ActivityScores(network, scores, regulation_type=reg_type, min_hub=min_hub)
		if len(mrObj.candidates) == 0:
			raise Exception("Error: no canidate hubs at this min cutoff")
		# perform 1000 random permutations of the data to get significance scores for each
		result = mrObj.scoreCandidates(nperms)

		tfs_heats = {}
		all_pvals = []
		for (tf, score) in sorted(result.items(), key=lambda t: t[1][0]):
			pval = score[1]
			heat = score[0]
			tfs_heats[tf] = (pval, heat)
			all_pvals.append(pval)

		t_total = 0
		for (g, h) in tfs_heats.items():
			t_total += abs(float(h[1]))

		# normalize abs values to sum to 1	
		norm_factor = 1.0/t_total
	
		for (g, h) in tfs_heats.items():
			tfs_heats[g] = (h[0], h[1]*norm_factor)

		# get BHFDRs
		corrected_p = ActivityScores.getBHYfdr(all_pvals)
		for (tf, heat) in sorted(tfs_heats.items(), key=lambda t: t[1][0]):
			pval = heat[0]	
			corrected = corrected_p[pval]
			tfs_heats[tf] = (heat[0], corrected, heat[1])
		
		return tfs_heats

		bfdr_neg = {}

	@staticmethod
	def getBHYfdr(pvals):
		# BHY-FDR procedure:
		m = float(len(pvals)*len(pvals))
		# try a 0.1 fdr
		pvals = sorted(pvals)
		corrected_vals = [1 for i in range(0, len(pvals))]
		for x in range(0, 20):

			alpha = x / 20.0
			# sort by p-value, then step through k=1... such that the equation is satisfied
			# Doing the Yekutieli procedure, where positive dependence is assumed
			k = 1
			i = 0
			for pval in pvals:

				if float(pval) > (k*alpha)/m:
					break

				if alpha < corrected_vals[i]:
					corrected_vals[i] = alpha

				k += 1
				i += 1
		
		ch = {}
		for i in range(0, len(pvals)):
			ch[pvals[i]] = corrected_vals[i]

		return ch

	@staticmethod
	def getPval(real, background):

		count = 0.0
		empirical_pval = None
		if real >= 0:
			# sort in descending order
			for val in sorted(background, reverse=True):
				if val >= real:
					count += 1
				else:
					break
			empirical_pval = (count+1)/(len(background)+1)
		else:
			# ascending order	
			for val in sorted(background, reverse=False):
				if val <= real:
					count += 1
				else:
					break
			empirical_pval = (count+1)/(len(background)+1)

		return empirical_pval

		
	def scoreCandidatesFISHER(self):
	
		scores = {}	
		for c in self.candidates:
			pos, neg = self.candidates[c]
			pval = self.scoreCHISQ(pos, neg)
			scores[c] = pval

		return scores	

	def scoreCandidates(self, nperms=1000):
	
		scores = {}	
		for c in self.candidates:
			pos, neg = self.candidates[c]
			score = self.scoreReg(pos, neg)
			bg = self.generateBackground(c, int(nperms))
			pval = ActivityScores.getPval(score, bg)
			scores[c] = (score, pval)
			# filter first by p-value, then weight by the score
		#	if pval < threshold:
		#		scores[c] = (score, pval)

		return scores	

	def generateBackground(self, candidate, nperms):
	
		pos, neg = self.candidates[candidate]
		# sample of this set size
		# of random genes to generate each permutation
		background_scores = []
		for i in range(0, nperms):
			sampled_pos = None
			sampled_neg = None
			if len(pos) < len(self.gene_list):
				sampled_pos = set(random.sample(self.gene_list, len(pos)))
			else:
				sampled_pos = self.gene_list

			if len(neg) < len(self.gene_list):
				sampled_neg = set(random.sample(self.gene_list, len(neg)))
			else:
				sampled_neg = self.gene_list

			score = self.scoreReg(sampled_pos, sampled_neg)
			background_scores.append(score)	

		return background_scores	

	def generateCategories(self, scores):
		"""
		Used for fisher's exact test: bin positive and negative sets 
		of each
		"""

		# create two sets: genes are either significantly up or down in either
		self.pos_de_set = set()
		self.neg_de_set = set()

		for (gene, score) in sorted(scores.iteritems(), key=operator.itemgetter(1), reverse=True):
			if score > 0:
				self.pos_de_set.add(gene)
			else:
				self.neg_de_set.add(gene)
	

	def generateRankings(self, scores):

		"""
			scores: scores of differential gene expression. These canonically are 
			d-statistic values output from Significance of Microarrays (SAM, Tishirani 2003).
			Input as a hash-map.
			Store the results in the internal index
		"""

		# invert the list, and then merge the postive and negative lists
		# descending order

		# save this data	
		self.gene_list = []
		self.scores = scores

		forward_genes = []
		forward_scores = []
		for (gene, score) in sorted(scores.iteritems(), key=operator.itemgetter(1), reverse=True):
			self.gene_list.append(gene)
			forward_genes.append(gene)
			forward_scores.append(score)
		# ascending order
		reverse_genes = []
		reverse_scores = []
		for (gene, score) in sorted(scores.iteritems(), key=operator.itemgetter(1), reverse=False):
			reverse_genes.append(gene)
			reverse_scores.append(score)

		# maintain two indexes
		indexF = 0
		indexR = 0
		# index by (gene, class (positive or negative))
		R_c = []
		# scores are combined
		R_c_SCORES = []
		while True:

			# termination conditions
			if indexF >= len(forward_genes) and indexR >= len(reverse_genes):
				break
			# append from the other list if one is finished
			elif indexF >= len(forward_genes):
				# the gene name and set are indexed
				R_c.append( (reverse_genes[indexR], '-') )
				R_c_SCORES.append( -reverse_scores[indexR] )
				indexR += 1
				continue
			elif indexR >= len(reverse_genes):
				R_c.append( (forward_genes[indexF], '+') )
				R_c_SCORES.append( forward_scores[indexF] )
				indexF += 1
				continue
					
			f_score = forward_scores[indexF]
			# inverse score...
			r_score = -reverse_scores[indexR]

			if f_score > r_score:
				R_c.append( (forward_genes[indexF], '+') )
				R_c_SCORES.append( f_score )
				indexF += 1
			else:
				R_c.append( (reverse_genes[indexR], '-') )
				R_c_SCORES.append( r_score )
				indexR += 1

		self.scores = R_c_SCORES
		self.list = R_c	


	def scoreCHISQ(self, pos_query_set, neg_query_set):
		"""
		Use chisquare approximation to fisher's exact test
		to calculate p-values for each
		"""

		# compute probability of random distribution in each
		# category by simple combinatorics
		s1 = len(self.pos_de_set)
		s2 = len(self.neg_de_set)
		norm = float(s1+s2)
		s1 = s1/norm
		s2 = s2/norm
		# expected frequencies for each set
		expected = np.array([s1, s2])

		up_AGREE = float(len(pos_query_set.intersection(self.pos_de_set)))
		up_DISAGREE = float(len(pos_query_set.intersection(self.neg_de_set)))
		observed = np.array([up_AGREE, up_DISAGREE])
		UP_chisq, UP_pval = stats.chisquare(observed, expected)
		
		down_AGREE = float(len(neg_query_set.intersection(self.neg_de_set)))
		down_DISAGREE = float(len(neg_query_set.intersection(self.pos_de_set)))
		observed = np.array([down_DISAGREE, down_AGREE])
		DOWN_chisq, DOWN_pval = stats.chisquare(observed, expected)

		combined_p = UP_pval*DOWN_pval
		return combined_p	

	def scoreReg(self, pos_query_set, neg_query_set):
		"""
			
		"""

		# from Lim et al., 2009 PSB
		rs_const = float(2.0*len(self.scores)-(len(pos_query_set)+len(neg_query_set)))
		running_sum = 0.0

		# -- norm const
		sum_norm_const = 0.0
		for i in range(0, len(self.list)):
			gene, set = self.list[i]
			if (set == '-' and gene in neg_query_set) or (set == '+' and gene in pos_query_set):
				# compute the sum of abs values of all scores in this set to get a normalization 
				# constant at the end
				sum_norm_const += abs(self.scores[i])

		# FIXME: handle this with a better error indication
		if sum_norm_const == 0:
			return 0

		running_sum = 0.0
		max_rs = 0
		min_rs = 0
		for i in range(0, len(self.list)):

			gene, set = self.list[i]
			if (set == '-' and gene in neg_query_set) or (set == '+' and gene in pos_query_set):
				running_sum += self.scores[i]/sum_norm_const
			else:
				# score decreases in this case
				running_sum -= 1/rs_const	

			if running_sum > max_rs:
				max_rs = running_sum	
			elif running_sum < min_rs:
				min_rs = running_sum	

		return max_rs+min_rs


class SSActivityScores: 
	"""
		regulon: 
			{ tf/regulon: (pos_query_set), (neg_query_set), ...}

		phenotypes:
			{ sample: 1 (test sample), 0 (reference sample) }

		:
			{ data[sample][expression] }
	"""
	def __init__(self, regulon, data, test_samples, ref_samples, min_hub=10, nperms=1000):

		# indexed by sample, then by gene	
		self.data = data

		# get the genes
		self.genes = []
		for sample in self.data:
			for gene in self.data[sample]:
				self.genes.append(gene)
			break	

	
		# get the true signature
		#self.signature = self.getSignature(test_samples, ref_samples)
		#
		#self.nes = self.dsGSEA(regulon, self.signature)

		# compute null model with random split of reference samples compared against each 
		# fit distributions to each regulators scores. 
		print 'Generating null model with '+str(len(ref_samples))+' reference samples and '+str(nperms)+' permutations...'
		k = int(math.ceil(len(ref_samples)/2.0))
		self.null_nes_scores = defaultdict(list)

		pool = []
		for c in combinations(ref_samples, k):
			pool.append(c)

		nperms = int(nperms)
		bg_samples = None
		if nperms > len(pool):
			bg_samples = pool
		else:
			bg_samples = random.sample(pool, nperms)
		
		for ref_samples_A in bg_samples:
			ref_samples_B = ref_samples.difference(ref_samples_A)
			null_signature = self.getSignature(ref_samples_A, ref_samples_B)
			null_nes = self.dsGSEA(regulon, null_signature)
			for gene in null_nes:
				self.null_nes_scores[gene].append(null_nes[gene]) 
			print '.'
		print "done."

		# indexed by sample, then gene
		self.nes = {}
		print "generating signatures for each sample"
		for sample in self.data:
			sig = self.getSignature(set([sample]), ref_samples_B)
			nes = self.dsGSEA(regulon, sig)
			self.nes[sample] = nes
			print '.'
		print "done."

	def getNullScores(self):

		# indexed by gene
		return self.null_nes_scores

	def getScores(self):

		scores = {}
		for sample in self.nes:
			scores[sample] = {}
			for gene in self.nes[sample]:
				pval = SSActivityScores.getPVAL(self.nes[sample][gene], self.null_nes_scores[gene])
				z_score = SSActivityScores.scoreAgainstBG(self.nes[sample][gene], self.null_nes_scores[gene])
				scores[sample][gene] = (z_score, pval)
		
		return scores

	@staticmethod	
	def scoreAgainstBG(real_score, bg_scores):

		z_score = (real_score - np.mean(bg_scores)) / np.std(bg_scores)
		return z_score

	@staticmethod	
	def getPVAL(real_score, bg_scores):

		bg_size = len(bg_scores)

		no_gte = 0.0
		for val in sorted(bg_scores, reverse=True):
			if val >= real_score:
				no_gte += 1
			else:
				break

		empirical_pval = (no_gte+1)/(bg_size+1)

		return empirical_pval


	def rowMeans(self, sample_set):

		row_means = {}
		std = {}
		EPSILON = 0.001
		for gene in self.genes:
			sum = 0.0
			values = []
			for sample in sample_set:
				if sample not in self.data:
					raise Exception("Error: sample not found in data input:"+sample)
				elif gene not in self.data[sample]:
					raise Exception("Error: gene data not found for sample:"+sample+'\t'+gene)

				sum += self.data[sample][gene]
				values.append(self.data[sample][gene])

			row_means[gene] = sum/len(sample_set)
			std[gene] = np.std(values)+EPSILON

		return (row_means, std)

	def getSignature(self, test_samples, ref_samples, method='zscore'):

		# difference in row-means (or z-scores) between sample sets

		test_means, test_std = self.rowMeans(test_samples)	
		ref_means, ref_std = self.rowMeans(ref_samples)	

		signature = {}
		if method == 'zscore':
			for gene in test_means:
				# fixme: scale by difference in standard deviation of the gene in each group. 
				signature[gene] = (test_means[gene] - ref_means[gene])/(test_std[gene] + ref_std[gene])
			
		return signature

	@staticmethod	
	def dsGSEA(regulon, signature):
		"""
			signature: scores of differential gene expression. These canonically are 
			d-statistic values output from Significance of Microarrays (SAM, Tishirani 2003).
			Input as a hash-map.

			regulon:
				{key: (pos_query_set), (neg_query_set), ...}

			Rank the signature: both the positive and negative lists, combine into a single ranked
			list. Then match each regulon, compute a enrichment score, return the enrichment scores.

			Returns:
				a hash indexed by each tf/regulator name and corresponding enrichment scores 

		"""

		# invert the list, and then merge the postive and negative lists
		# descending order

		# save this data	
		gene_list = []

		forward_genes = []
		forward_scores = []
		forward_sorted = sorted(signature.iteritems(), key=operator.itemgetter(1), reverse=True)
		for (gene, score) in forward_sorted:
			gene_list.append(gene)
			forward_genes.append(gene)
			forward_scores.append(score)
		# ascending order
		reverse_genes = []
		reverse_scores = []
		for (gene, score) in reversed(forward_sorted):
			reverse_genes.append(gene)
			reverse_scores.append(score)

		# maintain two indexes
		indexF = 0
		indexR = 0
		# index by (gene, class (positive or negative))
		# tuple comtaining gene, 'class' sign
		R_c = []
		# scores are combined
		# corresponding scores
		R_c_SCORES = []
		while True:

			# termination conditions
			if indexF >= len(forward_genes) and indexR >= len(reverse_genes):
				break
			# append from the other list if one is finished
			elif indexF >= len(forward_genes):
				# the gene name and set are indexed
				R_c.append( (reverse_genes[indexR], '-') )
				R_c_SCORES.append( -reverse_scores[indexR] )
				indexR += 1
				continue
			elif indexR >= len(reverse_genes):
				R_c.append( (forward_genes[indexF], '+') )
				R_c_SCORES.append( forward_scores[indexF] )
				indexF += 1
				continue
					
			f_score = forward_scores[indexF]
			# inverse score...
			r_score = -reverse_scores[indexR]

			if f_score > r_score:
				R_c.append( (forward_genes[indexF], '+') )
				R_c_SCORES.append( f_score )
				indexF += 1
			else:
				R_c.append( (reverse_genes[indexR], '-') )
				R_c_SCORES.append( r_score )
				indexR += 1

		#.scores = R_c_SCORES
		#list = R_c	

		nes_scores = {}
		for regulator in regulon:
			pos_query_set, neg_query_set = regulon[regulator]

			# from Lim et al., 2009 PSB
			rs_const = float(2.0*len(signature)-(len(pos_query_set)+len(neg_query_set)))
			running_sum = 0.0
	
			# -- norm const
			sum_norm_const = 0.0
			for i in range(0, len(R_c)):
				gene, set = R_c[i]
				if (set == '-' and gene in neg_query_set) or (set == '+' and gene in pos_query_set):
					# compute the sum of abs values of all scores in this set to get a normalization 
					# constant at the end
					sum_norm_const += abs(R_c_SCORES[i])
	
			# FIXME: handle this with a better error indication
			if sum_norm_const == 0:
				return 0
	
			running_sum = 0.0
			max_rs = 0
			min_rs = 0
			for i in range(0, len(R_c)):
	
				gene, set = R_c[i]
				if (set == '-' and gene in neg_query_set) or (set == '+' and gene in pos_query_set):
					running_sum += R_c_SCORES[i]/sum_norm_const
				else:
					# score decreases in this case
					running_sum -= 1/rs_const	
	
				if running_sum > max_rs:
					max_rs = running_sum	
				elif running_sum < min_rs:
					min_rs = running_sum	

			# Again, Lim 2009 defined this..
			#nes = abs(max_rs)+abs(min_rs)
			nes = abs(max_rs)
			nes_scores[regulator] = nes
			#print regulator+'\t'+str(nes)

		return nes_scores	

class SSActivityScores_Naive:
	"""
	Sample-specific activity scores module. The approach is to compute a 'activity' score for each 
	hub and sample, by averaging the expression (vs normal) of it's regulon. A background distribution
	for each is then built from permuted data, and z-scores for each sample/hub pair are returned. 
	"""
	def __init__(self, network, expression_matrix, normal_expression_matrix=None, min_hub=10, restrict_samples=None, z_threshold=1.5, predefined_activities=None):
		"""
		predefined_activities: a tuple with two references--a sample-gene indexed double hash of predefined 
		activity scores (i.e. z-scores) that have been pre-computed or filtered by some other algorithm. 
		If we have these, don't infer activities. The second is a regulon definition, indexed by each regulon
		and referencing the children of each. 
		
		"""

		NET_parents, NET_children = getTFparents(network)
		tf_candidates = set()
		for tf in NET_children:
			if len(NET_children[tf]) >= min_hub:
				tf_candidates.add(tf)

		if len(tf_candidates) < 1 and not predefined_activities[0]:
			raise Exception("Error: couldn't find any TF candidate hubs at this threshold: "+str(min_hub))
	
		# get expression activity scores
		self.sample_expr = parseMatrix(expression_matrix, restrict_samples=restrict_samples)

		# split into normal and 
		self.reference_expr = None
		if normal_expression_matrix:
			self.reference_expr = parseMatrix(normal_expression_matrix, restrict_samples=restrict_samples)

		self.sample_activities = None
		self.regulons = None
		# Either compute activity scores from the network provided, or use the pre-defined scores and		# define the  regulon.
		if not predefined_activities or not predefined_activities[0]:
			self.sample_activities = getActivityScores(self.sample_expr, tf_candidates, NET_parents)
		else:
			self.sample_activities = predefined_activities[0]
			self.regulons = {}
			for sample in self.sample_activities:
				self.regulons[sample] = predefined_activities[1]

		# fitted distributions for the activity scores
		bg_stats = None
		background_expression = None
		if self.reference_expr and len(self.reference_expr) > 0:	
			# If we can compare to normal or some kind of reference expression data, then 
			# use that as the background
			background_expression = self.reference_expr
		else:
			# generate a background distribution based on permuted (real) data, 
			# then fit to a gaussian distribution 
			background_expression = self.permuteLabels(30)	

		# background distributions for activity scores
		bg_act = None
		bg_stats = None
		if not predefined_activities or not predefined_activities[0]:
			bg_act = self.getActivityScores(background_expression, tf_candidates, NET_parents)
			bg_stats, bg_activity = self.fitGeneData(bg_act)
			self.z_scores = self.getZScores(self.sample_activities, bg_stats, z_threshold)
		else:
			self.z_scores = self.sample_activities
		# background distributions for the gene expression data
		cis_expression_background, cis_activity = self.fitGeneData(background_expression)
	
		# z-scores based on fitted background data
		# indexed by sample, then regulator
		self.downstream_expr_z_scores = self.getZScores(self.sample_expr, cis_expression_background, z_threshold)

		# FIXME: if we have pre-defined regulons, remove anything without the gene expression
		# above the threshold
		if not self.regulons:
			self.regulons = {}
			# only look at genes downstream of active regulators in this patient. 
			# create a hash indexed by regulators, and pointing to the set of genes regulated. 
			for sample in self.z_scores:
				self.regulons[sample] = {}
				for regulator in self.z_scores[sample]:
					self.regulons[sample][regulator] = set()
					for child in NET_children[regulator]:
						if child in self.downstream_expr_z_scores[sample]:	
							self.regulons[sample][regulator].add(child)


		# basic sanity checks
		if self.z_scores is None:
			raise Exception("Error: no sample activity found on this network!")

		non_trivial_values = False
		for sample in self.z_scores:
			if self.z_scores[sample] is not None:
				non_trivial_values = True
				break

		if not non_trivial_values:
			raise Exception("Error: no sample activity found on this network!")

	def getDownstreamRegulon(self, z_threshold):
		return self.regulons

	def getDownstreamExpression(self, z_threshold):
		return self.downstream_expr_z_scores

	def getActivities(self):
		#FIXME: at a supplied FDR cutoff, return the z-cutoff needed to acheive that. 
		# have to get z-scores for normal samples against this distribution
		return self.z_scores

	def fitGeneData(self, bg_activity_scores):
		"""
		Fit a normal distribution to the activity scores for each gene.
	
			Input
				activity_scores: permuted scores by sample
			Output
				summary: Hash indexed by gene, with mean/sd summary statistics
	
		"""

		activity = {}
		for sample in bg_activity_scores:
			for gene in bg_activity_scores[sample]:
				if gene not in activity:
					activity[gene] = []
				activity[gene].append(bg_activity_scores[sample][gene])
	
		distributions = {} 
		for gene in activity:
			distributions[gene] = Dist(activity[gene], method="gaussian")
			
		return (distributions, activity)

	def getActivityScores(self, expr_data, tf_genes, tf_parents, binary_threshold=0): 
	
		# the number of counts per gene, per sample 
		counts = {}
		# store the mean scores of each gene, for each sample
		activities = {}
	
		for sample in expr_data:
			activities[sample] = defaultdict(float)
			counts[sample] = defaultdict(int)
			for gene in expr_data[sample]:	
	
				val = expr_data[sample][gene]
	
				if abs(val) >= binary_threshold:
	
					if gene not in tf_parents:
						continue
	
					# check, is this downstream of a TF of interest? 
					# if so, add the TF, not the gene
					parents, activation_type = tf_parents[gene]
					for parent in tf_genes.intersection(parents):
					
						act = activation_type[parent]	
						tf_act = None
						# is this TF active? 
						if act == 'a':
							tf_act = val
						elif act == 'i':	
							tf_act = -1*val
	
						# add the activity, and the count
						activities[sample][parent] += tf_act
						counts[sample][parent] += 1
	
		# convert sums to means
		for sample in activities:
			for gene in activities[sample]:
				activities[sample][gene] = activities[sample][gene]/float(counts[sample][gene])
	
		return activities
	
	def getZScores(self, test_data, stats_summary, z_threshold, num_tests=None):
		"""
		Return z-scores for activity against the supplied background
		distribution (summary statistics)	
		Input
			scores: a hash indexed by genes, values are activity scores
			stats_summary: summary statistics of the per gene distributions
		
		Returns:
			activity_z_scores: z-scores for each gene
		"""
		
		z_scores = {}

		for sample in test_data:
			z_scores[sample] = {}
			for gene in test_data[sample]:

				if gene not in stats_summary:
					z_scores[sample][gene] = 0
					continue
	
				val = test_data[sample][gene]
				
				p = stats_summary[gene].getP(val)
				z = stats_summary[gene].getZ(val)

				# only save those that pass the threshold
				if abs(z) < z_threshold:
					continue
	
				#activity_z_scores[sample][gene] = (z, p, corrected_p)
				z_scores[sample][gene] = z
	
		return z_scores
	
	def permuteLabels(self, num_permuted_samples, by_gene=False):
		"""
		Permute gene-labels
		warning: this breaks the correlation between genes. 
		"""	
		data_by_gene = {}
		all_data = []
		for sample in self.sample_expr:
			for gene in self.sample_expr[sample]:
				if gene not in data_by_gene:
					data_by_gene[gene] = []
				data_by_gene[gene].append(self.sample_expr[sample][gene])
				all_data.append(self.sample_expr[sample][gene])
	
		permuted = {}
		for i in range(0, num_permuted_samples):
			permuted[i] = {}
			for gene in data_by_gene:
				vals = None
				if not by_gene:
					permuted[i][gene] = random.sample(all_data,1)[0]
				else:
					vals = data_by_gene[gene]
					# sample with replacement from that gene's data
					permuted[i][gene] = random.sample(vals,1)[0]
	
		return permuted
		
