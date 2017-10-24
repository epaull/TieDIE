#!/usr/bin/env	python

from optparse import OptionParser
parser = OptionParser()
(opts, args) = parser.parse_args()

from collections import defaultdict

tissue = args[0]

upstream_input_heats = defaultdict(float)

amps = set()
dels = set()
muts = set()

for line in open(tissue+'/DIFF_A.cnv.amps.txt', 'r'):
	protein = line.rstrip()
	amps.add(protein)	

for line in open(tissue+'/DIFF_A.cnv.dels.txt', 'r'):
	protein = line.rstrip()
	dels.add(protein)	

for line in open(tissue+'/DIFF_A.muts.txt', 'r'):
	protein = line.rstrip()
	muts.add(protein)	

total_heat = 1000

for amp in amps:
	upstream_input_heats[amp] += (total_heat/3.0) * (1.0/len(amps))

for delet in dels:
	upstream_input_heats[delet] += (total_heat/3.0) * (1.0/len(dels))

for mut in muts:
	upstream_input_heats[mut] += (total_heat/3.0) * (1.0/len(muts))

for protein in upstream_input_heats:
	print (protein+'\t'+str(upstream_input_heats[protein]))
