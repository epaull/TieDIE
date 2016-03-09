#!/usr/bin/env python

import unittest
import os

TEST_PATHWAY = "test_files/test.pathway.sif"
TEST_KERNEL = "test_files/kernel.tab"
TEST_SOURCE = "test_files/upstream.input"
TEST_TARGET = "test_files/downstream.input"
TEST_DIFFUSED_SOLN = "test_files/upstream.diffused"

class TestSequenceFunctions(unittest.TestCase):

	def readFile(self, name):
		lines = []
		for line in open(name, 'r'):
			lines.append(line.rstrip())

		return lines

	def filesEqual(self, file1, file2):

		f1_lines = self.readFile(file1)
		f2_lines = self.readFile(file2)

		if len(f1_lines) != len(f2_lines):
			return False
	
		for i in range(0, len(f1_lines)):
			if f1_lines[i] != f2_lines[i]:
				return False

		return True

	def cleanup(self, out_dir):
		os.system("rm -rf "+out_dir)

	def testRUN(self):

		test_dir = '/tmp'		
		cmd = "../bin/tiedie -n "+TEST_PATHWAY+" -k "+TEST_KERNEL+" "+TEST_SOURCE+" "+TEST_TARGET+" -s 1.0 --output "+test_dir
		print cmd
		os.system(cmd)
	

		# file output
		reg_dir = "test_files/REGRESSION/"
		self.assertTrue( self.filesEqual(test_dir+"/TieDIE.sif", reg_dir+"TieDIE.sif") )
		self.assertTrue( self.filesEqual(test_dir+"/heats.tab", reg_dir+"heats.tab") )

		# these are stochastic based on random seed: need to add method to set seed:

		#self.assertTrue( self.filesEqual(test_dir+"/score.txt", reg_dir+"score.txt") )
		#self.assertTrue( self.filesEqual(test_dir+"/edge_frequencies.txt", reg_dir+"edge_frequencies.txt") )
		#self.assertTrue( self.filesEqual(test_dir+"/node_frequencies.txt", reg_dir+"node_frequencies.txt") )
	
if __name__ == '__main__':
    unittest.main()

