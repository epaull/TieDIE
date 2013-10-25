#!/usr/bin/env python

import sys, os
sys.path.append(os.path.dirname(sys.argv[0])+'/../lib')
import unittest
from master_reg import ActivityScores
from tiedie_util import *
import time

TEST_PATHWAY = "test_files/PSN/pathway.sif"
TEST_KERNEL = "test_files/PSN/kernel.tab"
TEST_DE = "test_files/PSN/expr.ranked.tab"
TEST_EXPR = "test_files/PSN/expr.data"
TEST_MUT = "test_files/PSN/mut.data"
TEST_SOURCE = "test_files/PSN/upstream.input"

GENERATED_HEATS = "TieDIE/heats.tab"

class TestSequenceFunctions(unittest.TestCase):

	def setUp(self):
		self.startTime = time.time()

	def tearDown(self):
		t = time.time() - self.startTime
		print "%s: %.3f" % (self.id(), t)

	def testRUN(self):
		# build the concensus network first
		cmd = "../bin/tiedie "+\
			" -n "+TEST_PATHWAY+\
			" -k "+TEST_KERNEL+\
			" --d_expr "+TEST_DE+\
			" -m "+"3"+\
			" "+TEST_SOURCE
		print "Running command:"
		print cmd
		os.system(cmd)

		print "Testing Patient Networks against Cohort"
		cmd = "../bin/tiedie.PSN "+\
			" -n "+TEST_PATHWAY+\
			" -k "+TEST_KERNEL+\
			" --p_expr "+TEST_EXPR+\
			" --p_mut "+TEST_MUT+\
			" --node_ranks "+GENERATED_HEATS+\
			" -m "+"3"
		print "Running command:"
		print cmd
		os.system(cmd)


if __name__ == '__main__':
	unittest.main()


