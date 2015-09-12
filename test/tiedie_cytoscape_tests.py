#!/usr/bin/env python

import unittest
import sys, os

sys.path.append(os.path.dirname(sys.argv[0])+'/../lib')
from kernel import Kernel
from ppr import PPrDiffuser
from permute import NetBalancedPermuter
from tiedie_util import *
from kernel_scipy import SciPYKernel

TEST_PATHWAY = "test_files/cytoscape/pathway.sif"
#TEST_KERNEL = "test_files/cytoscape/scipy_results/kernel.tab"
TEST_INPUT = "test_files/cytoscape/upstream.input"
TEST_DIFFUSED_SOLN = "test_files/cytoscape/upstream.diffused"

JAVA_KERNEL = "test_files/cytoscape/java_results/exp.txt"

class TestSequenceFunctions(unittest.TestCase):

	def matrixComp(self, data1, data2):

		for sample in data1:
			self.assertTrue(sample in data2)
			for gene in data1[sample]:
				self.assertAlmostEqual(data1[sample][gene], data2[sample][gene], places=5)

		for sample in data2:
			self.assertTrue(sample in data1)
			for gene in data2[sample]:
				self.assertAlmostEqual(data1[sample][gene], data2[sample][gene], places=7)
	
	def testKernel(self):

		# make kernel out of the pathway
		diffuser = SciPYKernel(TEST_PATHWAY)
		scipy_kernel = diffuser.getKernelMatrix()

		# parse java pre-computed kernel and compare...
		java_kernel = parseMatrix(JAVA_KERNEL)
		self.matrixComp(scipy_kernel, java_kernel)

if __name__ == '__main__':
    unittest.main()

