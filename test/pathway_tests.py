#!/usr/bin/env python

import unittest
import sys, os

sys.path.append(os.path.dirname(sys.argv[0])+'/../lib')
from pathway import *

TEST_PATHWAY = "test_files/test.pathway.sif"

class TestSequenceFunctions(unittest.TestCase):
	
	def testDiffuseKernel(self):

		validator = TrivialPathValidator()
		pathwayObj = Pathway(validator=validator)
		pathwayObj.parseNet(TEST_PATHWAY)
		self.assertEqual(0, len(pathwayObj.allPaths(set(['CDK4', 'PI3KCA', 'PIP3', 'PTEN']), set(['TP53', 'cell_cycle_progression', 'INK4']), 1)))
		self.assertEqual(2, len(pathwayObj.allPaths(set(['CDK4', 'PI3KCA', 'PIP3', 'PTEN']), set(['TP53', 'cell_cycle_progression', 'INK4']), 2)))
		self.assertEqual(4, len(pathwayObj.allPaths(set(['CDK4', 'PI3KCA', 'PIP3', 'PTEN']), set(['TP53', 'cell_cycle_progression', 'INK4']), 3)))
		self.assertEqual(8, len(pathwayObj.allPaths(set(['CDK4', 'PI3KCA', 'PIP3', 'PTEN']), set(['TP53', 'cell_cycle_progression', 'INK4']), 4)))
		self.assertEqual(9,len(pathwayObj.allPaths(set(['CDK4', 'PI3KCA', 'PIP3', 'PTEN']), set(['TP53', 'cell_cycle_progression', 'INK4']), 5)))

if __name__ == '__main__':
    unittest.main()

