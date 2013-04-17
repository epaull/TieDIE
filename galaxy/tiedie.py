#!/usr/bin/env python

import os, sys
from collections import defaultdict
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-u","--up_heats",dest="up_heats",action="store",default=None,help="\
File with upstream heats: <gene>    <input heat (0-100)>    <sign (+/-)>")
parser.add_option("-d","--down_heats",dest="down_heats",action="store",default=None,help="\
File with downstream heats: <gene>    <input heat (0-100)>    <sign (+/-)>")
parser.add_option("-k","--kernel",dest="kernel",action="store",default="superpathway.tab",help="\
Pre-computed heat diffusion kernel in tab-delimited form. Should have both a header and row labels")
parser.add_option("--pcst",dest="pcst",action="store_true",default=False,help="\
Use the Prize-Collecting Steiner Tree Formulation to Generate a Connecting Subnetwork (Bionet must be installed)")
parser.add_option("-n","--network",dest="network",action="store",default="superpathway.sif",help="\
.sif network file for the curated pathway to search. <source>   <(-a>,-a|,-t>,-t|,-component>)> <target>")
parser.add_option("-s","--size",dest="size",action="store",default="1",help="\
Network size control factor (default 1)")
parser.add_option("-a","--alpha",dest="alpha",action="store",default=None,help="Linker Cutoff (overrides the Size factor)")
parser.add_option("-c","--depth",dest="depth",action="store",default="3",help="\
Search depth for causal paths (default 3)")
parser.add_option("-p","--permute",dest="permute",action="store",default=1000,help="\
Number of random permutations performed for significance analysis (default 1000)")
parser.add_option("--pagerank",dest="pagerank",action="store_true",default=False,help="Use Personalized PageRank to Diffuse")
parser.add_option("--all_paths",dest="all_paths",action="store_true",default=False)
parser.add_option("--select_size",dest="size_select",action="store",default=None)
(opts, args) = parser.parse_args()

# local imports
from kernel import Kernel
from ppr import PPrDiffuser
from permute import NetBalancedPermuter
from size_selector import SizeSelector
from tiedie_util import *

# Program Constants
SCORE_MU = 0.1

def extractSubnetwork(up_heats, down_heats, up_heats_diffused, down_heats_diffused, size_control, set_alpha):
	"""
		Generate a spanning subnetwork from the supplied inputs, diffused heats and 
		size control cutoff

		Input:
			- upstream heats
			- downstream heats
			- diffused upstream heats
			- diffused downstream heats	
			- size control factor

		Output:
			- spanning network
			- list of nodes in that network
	"""

	# find linker cutoff
	linker_cutoff = None
	linker_nodes = None
	linker_scores = None
	alpha_score = None
	# if explicitly set, use the supplied linker cutoff
	if set_alpha:
		alpha_score = None
		linker_cutoff = float(set_alpha)
	else:
		linker_cutoff, alpha_score = findLinkerCutoff(up_heats, down_heats, up_heats_diffused, down_heats_diffused, size_control)

	linker_nodes, linker_scores = filterLinkers(up_heats_diffused, down_heats_diffused, linker_cutoff)

	ugraph = None
	if opts.pcst:
		ugraph = runPCST(up_heats, down_heats, linker_nodes, opts.network)
	else:
		nodes = set(up_heats).union(set(down_heats)).union(set(linker_nodes))
		ugraph = connectedSubnets(network, nodes)

	if len(ugraph) == 0:
		print "Couldn't find any linking graph at this size setting!"
		return (None, None, None, None)
	subnet_soln = mapUGraphToNetwork(ugraph, network)
	
	subnet_soln_nodes = set()
	for s in subnet_soln:
		subnet_soln_nodes.add(s)
		for (i,t) in subnet_soln[s]:
			subnet_soln_nodes.add(t)

	return (subnet_soln, subnet_soln_nodes, alpha_score, linker_scores)

def findConsistentPaths(up_signs, down_signs, searchNetwork, output_folder, output):

	"""
		Input:
			- hash with up/down signs for each upstream node
			- hash with up/down signs for each downstream node
			- subnetwork to search over 
			
		Options:
			- output: flag to write output 

	"""
	# connect all sources
	gene_states, t_states = classifyState(up_signs, down_signs)
	validated = set()
	down_set = set(down_signs.keys())
	# if we're doing a randomized link analysis, keep track of TP and FP scores
	TP = 0
	FP = 0
	for source in up_signs:
		action = gene_states[source]
		falsePaths = []
		truePaths = []	
		edges_this_source = set()
		searchDFS(source, action, edges_this_source, set(), down_set, searchNetwork, gene_states, t_states, search_depth, truePaths, falsePaths, False)
		for edge in edges_this_source:
			validated.add(edge)
		TP += len(truePaths)	
		FP += len(falsePaths)	

		# write it
		if output:
			out_file = output_folder+"/"+source+".cn.sif"
			print "Writing Single Causal Neighborhood to "+out_file
			writeEL(edges_this_source, source, down_set, out_file)

	if output:	
		out_file = output_folder+"/tiedie.cn.sif"
		print "Writing Full Causal Neighborhood to "+out_file
		writeEL(validated, "ALL", down_set, out_file)

	return (TP, FP, validated)

def scoreSubnet(subnet_soln_nodes, up_heats, down_heats):
	# Score Sets According to Compactness Score
	S = set(up_heats.keys())
	T = set(down_heats.keys())
	C = subnet_soln_nodes.difference(S).difference(T)
	U = S.union(T)
	Sr = S.intersection(subnet_soln_nodes)
	Tr = T.intersection(subnet_soln_nodes)
	penalty = float(len(C))/len(U)
	score = float(len(Sr))/(len(S)*2) + float(len(Tr))/(len(T)*2) - SCORE_MU*penalty

	return score

def sampleDiffuse(input_heats, input_signs, diffuser, rate, reverse):
	"""
		Sample without replacement at the specified rate. Use the subsampled 
		values to diffuse and return both the subsampled set and the diffused heats
 		from it
	"""
	sample_size	= int(len(input_heats)*float(rate))

	subsampled_heats = {}
	subsampled_signs = {}
	for node in random.sample(input_heats, sample_size):
		subsampled_heats[node] = input_heats[node]
		subsampled_signs[node] = input_signs[node]

	diffused_heats = diffuser.diffuse(subsampled_heats, reverse)

	return (subsampled_heats, subsampled_signs, diffused_heats)

# parse inputs, set up output folder:
PERMUTE = int(opts.permute)
opts.rate = 1.0

up_heats, up_signs = parseHeats(opts.up_heats)
down_heats, down_signs = parseHeats(opts.down_heats)
size_control = float(opts.size)
search_depth = int(opts.depth)

# output options for Galaxy

out_prefix = os.path.dirname(opts.up_heats)
if out_prefix == "":
	out_prefix = "."
output_folder = out_prefix+"/TieDIE_RESULT_size="+opts.size+"_depth="+opts.depth
if opts.alpha:
	output_folder = out_prefix+"/TieDIE_RESULT_"+str(opts.alpha)
if opts.pagerank:
	output_folder += "_PAGERANK"
		
# create directory
if not os.path.exists(output_folder):
	os.mkdir(output_folder)

#
# Diffusion Step:
#	Load the heat diffusion kernel and perform a kernel-multiply, or alternatively use supplied
# 	page-rank diffused vectors
#

print "Parsing Network File..."
network = parseNet(opts.network)
if opts.pagerank:
	diffuser = PPrDiffuser(network)
else:
	print "Loading Heat Diffusion Kernel.."
	diffuser = Kernel(opts.kernel)

print "Diffusing Heats..."
up_heats_diffused = diffuser.diffuse(up_heats, reverse=False)
down_heats_diffused = diffuser.diffuse(down_heats, reverse=True)

if opts.size_select:

	# Search a number of network sizes to get a range with the highest significance. 
	# Somewhat slow. 
	search_space = []
	start, end = opts.size_select.split(":")
	start = int(start)
	end = int(end)
	for size in range(1, 20):
		search_space.append(float(size)/20)
	for size in range(11, 20):
		search_space.append(float(size)/10)
	for size in range(1, 15):
		search_space.append(2+float(size)/5)
	for size in range(1, 20):
		search_space.append(5+float(size)/4)
	# 2->5 small, 5->

	sizeSelector = SizeSelector(network, up_heats, down_heats, down_heats_diffused, diffuser, PERMUTE)	
	sizeSelector.generateBackground(search_space)
	graph = sizeSelector.getGraph()

	min = (1,0)
	for size in search_space:
		score, pval = graph[size]
		if pval < min[0]:
			min = (pval, size)

	search_space = []
	start = int(min[1])-1
	end = int(min[1])+1
	for size in range(start*10+1, end*10):
		search_space.append(float(size)/10)

	sizeSelector.generateBackground(search_space)
	graph = sizeSelector.getGraph()
	for size in search_space:
		score, pval = graph[size]
		if pval < min:
			min = (pval, size)
	
	sys.exit(0)

#
# Extract a subnetwork solution from the diffused heats
#
subnet_soln, subnet_soln_nodes, alpha_score, linker_scores = extractSubnetwork(up_heats, down_heats, up_heats_diffused, down_heats_diffused, size_control, opts.alpha)
#print "Real RI Score:\t"+str(alpha_score)

# 
# Generate linker stats and output
# 
out_degrees = getOutDegrees(subnet_soln)
print "Writing network node stats to "+output_folder+"/node.stats"
out_file = output_folder+"/node.stats"
out = open(out_file, 'w')
out.write("NODE\tCONNECTING\tMIN_HEAT\tOUT_DEGREE\n")
for node in subnet_soln_nodes:
	out_deg = out_degrees[node]
	linker_heat = linker_scores[node]
	connecting = "0" 
	if node not in up_heats:
		if down_heats is not None and node not in down_heats:
			connecting = "1"	
	out.write("\t".join([node, connecting, str(linker_heat), str(out_deg)])+"\n")
out.close()

print "Writing "+output_folder+"/tiedie.sif result "
writeNetwork(subnet_soln, args[0])

# 
# Find logically consistent paths
#
TP, FP, validated = findConsistentPaths(up_signs, down_signs, subnet_soln, output_folder, True)

# Open the report file
report_file = args[1]
report_fh = open(report_file, 'w')

# Score Sets According to Compactness Score
print "Writing Report to "+report_file+" :compactness analysis\n"
S = set(up_heats.keys())
T = set(down_heats.keys())
C = subnet_soln_nodes.difference(S).difference(T)
U = S.union(T)
Sr = S.intersection(subnet_soln_nodes)
Tr = T.intersection(subnet_soln_nodes)
report_fh.write(str(float(len(Sr))/len(S))+"\t"+"of source nodes"+str(len(Sr))+" out of "+str(len(S))+"\n")
report_fh.write(str(float(len(Tr))/len(T))+"\t"+"of target nodes"+str(len(Tr))+" out of "+str(len(T))+"\n")
report_fh.write("And "+str(len(C))+" connecting nodes\n")
penalty = float(len(C))/len(U)
score = float(len(Sr))/(len(S)*2) + float(len(Tr))/(len(T)*2) - SCORE_MU*penalty
report_fh.write("Compactness Score:"+str(score)+"\n")

print "Running permutation tests... (could take several minutes for inputs of hundreds of genes @1000 permutations)"
perObj = NetBalancedPermuter(network, up_heats)
permutedHeats = perObj.permute(PERMUTE)
permuted_scores = []
for heats in permutedHeats:
	# diffuse
	diffused_heats =  diffuser.diffuse(heats)
	cutoff, score = findLinkerCutoff(heats, down_heats, diffused_heats, down_heats_diffused, size_control)
	#print score
	permuted_scores.append(score)

# write out the distribution for significance plot
#sig_fh = open(output_folder+"/score.txt", 'w')
#sig_fh.write(str(alpha_score)+"\n")
#sig_fh.close()
#sig_fh = open(output_folder+"/permuted_scores.txt", 'w')
#for val in sorted(permuted_scores, reverse=True):
#	sig_fh.write(str(val)+"\n")
#sig_fh.close()

no_gte = 0.0
for val in sorted(permuted_scores, reverse=True):
	if val >= alpha_score:
		no_gte += 1
	else:
		break

# Davison & Hinkley (1997): true empirical p-value is (r+1)/(n+1)
print "Writing Report to "+report_file+" :empirical p-value..." 
pval = (no_gte+1)/(PERMUTE+1)
report_fh.write("P-value: "+str(pval)+" (with "+str(PERMUTE)+" random permutations)\n")
report_fh.close()
