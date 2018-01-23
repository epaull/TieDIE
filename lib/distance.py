from __future__ import print_function
from __future__ import division
#from __future__ import unicode_literals
import numpy as np

class ProbDistance:
	"""
		Holds functions for various distance metrics between probability vectors
		including symmetric and non-symmetric measures.
	"""

	@staticmethod
	def getSymmetricMeasure(v1, v2, measure="kl_div", symmetric=True):
		"""
		Inputs:
			v1: first hash-key/value probability vector
			v2: second hash-key/value probability vector
			measure: distance metric
			symmetric: True if the symmetric version of a measure is required

		Outputs:
			Floating-point value of distance

		Kullback Lieber Divergence is the only supported measure at this point
		no immediate plans for future measures
		"""
		if measure == "kl_div" and symmetric:
			return Distance.getSYMKLDiv(v1, v2)
		elif measure == "kl_div" and not symmetric:
			return Distance.getKLDiv(v1, v2)
		else:
			raise Exception("Error, method "+measure+" not yet implemented")

	@staticmethod
    def getSYMKLDiv(v1, v2):
        """
        <Development module> Get the symmetric Kullback-Leibler divergence between input vectors.
        """
        return (Kernel.getKLDiv(v1, v2) + Kernel.getKLDiv(v2, v1) )/ 2

    @staticmethod
    def getKLDiv(v1, v2):
        """
        <Development module> Get the Kullback-Leibler divergence between input vectors,
        which typically represent diffused heat values.
        Input:
            2 diffused heat vectors
        Output:
            (float) The KL-divergence metric between vectors
        """

        # convert values to ordered array
        arry1 = []
        arry2 = []
        for key in v1:
            arry1.append(float(v1[key]))
            arry2.append(float(v2[key]))

        # normalize both vectors, adding an epsilon (pseudo-count) value for zero-values
        EPSILON = 0.00001
        norm_arry1 = [ a+EPSILON/(sum(arry1)+len(arry1)*EPSILON) for a in arry1 ]
        norm_arry2 = [ a+EPSILON/(sum(arry2)+len(arry2)*EPSILON) for a in arry2 ]

        # calculate the KL-Div
        div_sum = 0
        for i in range(0, len(norm_arry1)):
            div_sum += np.log( (norm_arry1[i]/norm_arry2[i]) )*norm_arry1[i]

        return div_sum
