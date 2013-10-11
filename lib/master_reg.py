
from tiedie_util import classifyInteraction

class ActivityScores: 
	"""
		Uses the supplied pathway to find 	

	"""

	def __init__(self, network, min_hub=10):
        """
            Input:
                network: net[source] = [(i, t)]
				ranked_list: ranked list of differential gene expression
				min_hub: minimum number of genes regulated transcriptionally required 
				to be considered as a potential 'master regulator'
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

			if len(positive_regulon) + len(negative_regulon) >= min_hub:
				self.candidates[source] = (positive_regulon, negative_regulon)


	def scoreRegulators(self, scores):

		"""
			scores: scores of differential gene expression. These canonically are 
			d-statistic values output from Significance of Microarrays (SAM, Tishirani 2003).
			Input as a hash-map.
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
		reverse_sorted = sorted(scores.iteritems(), key=operator.itemgetter(1), reverse=False)
			reverse_genes.append(gene)
			reverse_scores.append(score)

		# maintain two indexes
		indexF = 0
		indexR = 0
		# index by (gene, class (positive or negative))
		R_c = []
		# scores are combined
		R_c_SCORES = []
		R_c_index = 0	
		while True:

			# termination conditions
			if indexF >= len(forward_genes) and indexR >= len(reverse_genes):
				break
			# append from the other list if one is finished
			elif indexF >= len(forward_genes):
				R_c.append( (reverse_genes[indexR], '-') )
				R_c_SCORES.append( reverse_scores[indexR] )
				indexR += 1
			elif indexR >= len(reverse_genes):
				R_c.append( (forward_genes[indexF], '+') )
				R_c_SCORES.append( forward_scores[indexF] )
				indexF += 1
					
			f_score = forward_scores[indexF]
			# inverse score...
			r_score = -forward_scores[indexR]

			if f_score > r_score:
				R_c.append( (forward_genes[indexF], '+') )
				R_c_SCORES.append( f_score )
				indexF += 1
			else:
				R_c.append( (reverse_genes[indexR], '-') )
				R_c_SCORES.append( r_score )
				indexR += 1
	

class ExtendedGSEA:


	@staticmethod
	def score(ranked, pos_query_set, neg_query_set):
		"""
			
		"""





