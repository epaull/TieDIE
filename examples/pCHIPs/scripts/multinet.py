#!/usr/bin/env	python

head = None
for line in open('pathways/multinet.full.txt', 'r'):

	data = line.rstrip().split('\t')
	if not head:
		head = data
		continue


	data = dict(zip(head, data))

	geneA, geneB = data['INTERACTION_NAME'].split('_')
	for interaction in data.keys()[1:len(data)]:
		if interaction=='NUM_NETWORKS':
			continue
		if data[interaction]=='1':
			print '\t'.join([geneA, interaction, geneB])
		
