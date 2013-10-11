
from tiedie_util import classifyInteraction
import operator

class ActivityScores: 
	"""
		Uses the supplied pathway to find 	

	"""

	def __init__(self, network, scores, min_hub=10, p=1):
		"""
			Input:
				network: net[source] = [(i, t)]
				scores: hash map of differential gene expression (think D-statistics from SAM)
				min_hub: minimum number of genes regulated transcriptionally required 
				to be considered as a potential 'master regulator'
				p: the power to raise each element to when computing the running sum
		"""

		# build a list of candidate regulators
		self.candidates = {}	

		for source in network:

			positive_regulon = set()
			negative_regulon = set()
			for (i, t) in network[source]:
				type, mode = classifyInteraction(i)	
				# only consider transcriptional regulation
				if mode != 't':
					continue

				if type == 1:
					positive_regulon.add(t)
				elif type == -1:
					negative_regulon.add(t)

			if (len(positive_regulon) + len(negative_regulon)) >= min_hub:
				self.candidates[source] = (positive_regulon, negative_regulon)

		self.generateRankings(scores)

	def scoreCandidates(self):
	
		scores = {}	
		for c in self.candidates:
			pos, neg = self.candidates[c]
			score = self.scoreReg(pos, neg)
			scores[c] = score

		return scores	

	def generateRankings(self, scores):

		"""
			scores: scores of differential gene expression. These canonically are 
			d-statistic values output from Significance of Microarrays (SAM, Tishirani 2003).
			Input as a hash-map.
			Store the results in the internal index
		"""

		# invert the list, and then merge the postive and negative lists
		# descending order
		forward_genes = []
		forward_scores = []
		for (gene, score) in sorted(scores.iteritems(), key=operator.itemgetter(1), reverse=True):
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

		# build the total sum
		sum = 0
		for score in self.scores:
			sum += abs(score)

		self.norm_const = sum

	def scoreReg(self, pos_query_set, neg_query_set):
		"""
			
		"""

		# from Lim et al., 2009 PSB
		self.rs_const = float(2*len(self.scores)-(len(pos_query_set)+len(neg_query_set)))
		running_sum = 0.0
		max_rs = 0
		min_rs = 0
		for i in range(0, len(self.list)):

			gene, set = self.list[i]
			if (set == '-' and gene in neg_query_set) or (set == '+' and gene in pos_query_set):
				running_sum += self.scores[i]/self.norm_const
			else:
				# score decreases in this case
				running_sum -= 1/self.rs_const	

			if running_sum > max_rs:
				max_rs = running_sum	
			elif running_sum < min_rs:
				min_rs = running_sum	
	
		return max_rs+min_rs
