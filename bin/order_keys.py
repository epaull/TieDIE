#!/usr/bin/env	python

import sys

pairs = set()
for line in sys.stdin:
	A, B = line.rstrip().split('\t')

	if (B, A) in pairs:
		continue

	pairs.add( (A, B) )

for pair in sorted(pairs):
	print ('\t'.join(pair))
