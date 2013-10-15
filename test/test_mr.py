#!/usr/bin/env python

import sys, os
sys.path.append(os.path.dirname(sys.argv[0])+'/../lib')
import unittest
from master_reg import ActivityScores
from tiedie_util import *
import time

TEST_PATHWAY = "test_files/test.tfnet.sif"
TEST_DATA = "test_files/test.tfnet.data.tab"

P1 = "test_files/test.tfnetbig.sif"
P1_D = "test.de"
#P1_D = "test_files/test.thcaBRAF.vs.RAS.tab"

class TestSequenceFunctions(unittest.TestCase):

	def setUp(self):
		self.startTime = time.time()

	def tearDown(self):
		t = time.time() - self.startTime
		print "%s: %.3f" % (self.id(), t)

	def testRUN(self):
		network = parseNet(TEST_PATHWAY)
		# signs is empty here
		scores, signs = parseHeats(TEST_DATA)
		mrObj = ActivityScores(network, scores, min_hub=3)
		valid_scores = [10.0, 8.0, 5.0, 4.05, 3.05, 1.05, -1.05, -3.05, -4.05, -5.0, -8.0, -10.0]
		valid_indexes = [('pos1', '+'), ('pos2', '+'), ('pos3', '+'), ('neg3', '-'), ('neg2', '-'), ('neg1', '-'), ('neg1', '+'), ('neg2', '+'), ('neg3', '+'), ('pos3', '-'), ('pos2', '-'), ('pos1', '-')]
		i = 0
		for i in range(0, len(valid_scores)):
			self.assertEqual(mrObj.scores[i], valid_scores[i])
			self.assertEqual(mrObj.list[i], valid_indexes[i])

		#print "Tiny Example:"
		#for (tf, score) in mrObj.scoreCandidates(nperms=100).items():
		#	print tf+"\t"+str(score)

		network = parseNet(P1)
		# signs is empty here
		scores, signs = parseHeats(P1_D)
		mrObj = ActivityScores(network, scores, min_hub=100)
		result = mrObj.scoreCandidates(nperms=100)
		for (tf, result) in sorted(result.items(), key=lambda t: t[1][0]):
			#if result > 0.05:
			#	continue
			print tf+"\t"+"\t".join([str(float(v)) for v in result])


if __name__ == '__main__':
	unittest.main()


