
from scipy.stats import norm
from math import log
   
class Dist: 

	"""
	This class holds various parameterized distributions used to fit empirical data
	generated by background models used by the TieDIE algorithm. 
	"""

	@staticmethod
	def fitLogNorm(vector, test_value):
		"""
			Fit a log-normal to the background distrubtion supplied. 
			Get the p-value based on the log of the test value. 

			Input:
				vector: background distribution to fit
				test_value: value to test against the fitted background distribution

			Output:
				A p-value based on that distribution	
		"""
		EPSILON = 0.001
		mean, sd = norm.fit([log(v) for v in vector])		
		# just the cdf: this value should be smaller 
		p_val = norm.cdf(log(test_value+EPSILON), loc=mean,scale=sd)		

		return p_val
