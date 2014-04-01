#!/usr/bin/env python

import unittest
import os

TEST_PATHWAY = "test_files/test.pathway.sif"
TEST_KERNEL = "test_files/kernel.tab"
TEST_SOURCE = "test_files/upstream.input"
TEST_TARGET = "test_files/downstream.input"
TEST_DE = "test_files/de.input"
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

		self.testIN1()
		self.testIN2()

	def testIN1(self):
		cmd = "../bin/tiedie -n "+TEST_PATHWAY+" -k "+TEST_KERNEL+" -u "+TEST_SOURCE+" -d "+TEST_TARGET+" -s 1.0"
		print cmd
		os.system(cmd)
	
		output_dir = "test_files/TieDIE_RESULT_size=1.0_depth=3/"
		reg_dir = "test_files/REGRESSION/"
		self.validateOutput(output_dir, reg_dir)

	def testIN2(self):
		cmd = "../bin/tiedie -n "+TEST_PATHWAY+" -k "+TEST_KERNEL+" -u "+TEST_SOURCE+" --d_expr "+TEST_DE+" -s 1.0 --min_hub 3"
		print cmd
		os.system(cmd)
	
		output_dir = "test_files/TieDIE_RESULT_size=1.0_depth=3/"
		reg_dir = "test_files/REGRESSION/"
		self.validateOutput(output_dir, reg_dir)

	def validateOutput(self, output_dir, reg_dir):
		# test report output
		report_text = self.readFile(output_dir+"report.txt")
		self.assertEqual(report_text[0].split("\t")[0], "0.375")
		self.assertEqual(report_text[1].split("\t")[0], "0.5")
		self.assertEqual(report_text[2], "And 9 connecting nodes")
		self.assertEqual(report_text[3], "Compactness Score:0.3475")

		# file output
		self.assertTrue( self.filesEqual(output_dir+"tiedie.sif", reg_dir+"tiedie.sif") )
		self.assertTrue( self.filesEqual(output_dir+"heats.NA", reg_dir+"heats.NA") )
		self.assertTrue( self.filesEqual(output_dir+"node.stats", reg_dir+"node.stats") )
		self.assertTrue( self.filesEqual(output_dir+"tiedie.cn.sif", reg_dir+"tiedie.cn.sif") )
		self.assertTrue( self.filesEqual(output_dir+"score.txt", reg_dir+"score.txt") )
	
		self.cleanup(output_dir)

if __name__ == '__main__':
    unittest.main()

