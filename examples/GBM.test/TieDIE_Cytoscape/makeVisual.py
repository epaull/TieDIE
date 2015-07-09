#!/usr/bin/env	python

import sys


from optparse import OptionParser
parser = OptionParser()
parser.add_option("-o","--output",dest="output",action="store")
parser.add_option("-i","--input_heats",dest="input_heats",action="store")
(opts, args) = parser.parse_args()

def parseInputHeats(file):

	input_heats = {}

	fh = open(file, 'r')
	for line in fh:
		node, heat, sign = line.rstrip().split('\t')
		input_heats[node] = heat
		
	return input_heats

def parseNet(file):

	edges = {}
	nodes = set()

	fh = open(file, 'r')
	for line in fh:
		a, I, b = line.rstrip().split('\t')
		edges[(a, b)] = I
		nodes.add(a)
		nodes.add(b)
		
	return (edges, nodes)	

inputs = parseInputHeats(opts.input_heats)

full_network, network_nodes = parseNet(args[0])
subnetwork, subnet_nodes = parseNet(args[1])


fh = open(opts.output+'/node_attributes.txt', 'w')
fh.write('name\tin_subnet\tinput_heat\n')
for node in network_nodes:
	printstr = node
	if node in subnet_nodes:
		printstr += '\t1'
	else:
		printstr += '\t0'

	if node in inputs:
		printstr += '\t'+inputs[node]

	fh.write(printstr+'\n')

fh.close()

fh = open(opts.output+'/edge_attributes.txt', 'w')
fh.write('\t'.join(['source', 'interaction', 'target', 'in_subnet'])+'\n')
for edge in full_network:
	i = full_network[edge]

	if edge in subnetwork:
		fh.write('\t'.join([edge[0], i, edge[1]])+'\t'+'1'+'\n')
	else:
		fh.write('\t'.join([edge[0], i, edge[1]])+'\t'+'0'+'\n')

fh.close()
