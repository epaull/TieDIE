#!/usr/bin/env python

import unittest
import sys, os

sys.path.append(os.path.dirname(sys.argv[0])+'/../lib')
from kernel import Kernel
from ppr import PPrDiffuser
from permute import NetBalancedPermuter
from tiedie_util import *
from kernel_scipy import SciPYKernel

TEST_PATHWAY = "test_files/test.pathway.sif"
TEST_KERNEL = "test_files/kernel.tab"
TEST_INPUT = "test_files/upstream.input"
TEST_DIFFUSED_SOLN = "test_files/upstream.diffused"

class TestSequenceFunctions(unittest.TestCase):
	
	def testDiffuseSciPY(self):

		correct_heats, s = parseHeats(TEST_DIFFUSED_SOLN)	
		input_heats, input_signs = parseHeats(TEST_INPUT)
		diffuser = SciPYKernel(TEST_PATHWAY)
		diffused = diffuser.diffuse(input_heats, reverse=False)
		for (key, val) in diffused.iteritems():
			self.assertAlmostEqual(val, correct_heats[key], places=10)

	def testDiffuseKernel(self):

		correct_heats, s = parseHeats(TEST_DIFFUSED_SOLN)	
		input_heats, input_signs = parseHeats(TEST_INPUT)
		diffuser = Kernel(TEST_KERNEL)
		diffused = diffuser.diffuse(input_heats, reverse=False)
		for (key, val) in diffused.iteritems():
			if val < 1:
				self.assertAlmostEqual(val, correct_heats[key], places=5)
			else:
				self.assertAlmostEqual(val, correct_heats[key], places=3)

if __name__ == '__main__':
    unittest.main()

