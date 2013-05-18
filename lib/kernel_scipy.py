#!/usr/bin/env  python2.7

from numpy import genfromtxt, dot
import sys
import math
from array import array
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm

class SciPYKernel:

	def __init__(self, network_file):
		""" 
			Input:
				kernel_file - a tab-delimited matrix file with both a header
				and first-row labels, in the same order. 
		"""

		# might have multiple kernels here
		self.labels = {}
		self.ncols = {}
		self.nrows = {}

		# parse the network, build the graph laplacian
		edges, nodes, node_out_degrees = self.parseNet(network_file)
		num_nodes = len(nodes)
		node_order = list(nodes)
		index2node = {}
		node2index = {}
		for i in range(0, num_nodes):
			index2node[i] = node_order[i]	
			node2index[node_order[i]] = i
	
		# construct the diagonals
		row = array('i')
		col = array('i')
		data = array('f')
		for i in range(0, num_nodes):
			# diag entries: out degree
			degree = 0 
			if index2node[i] in node_out_degrees:
				degree = node_out_degrees[index2node[i]]	
			# append to the end
			data.insert(len(data), degree)	
			row.insert(len(row), i)	
			col.insert(len(col), i)	

		# add edges
		for i in range(0, num_nodes):
			for j in range(0, num_nodes):
				if i == j:
					continue
				if (index2node[i], index2node[j]) not in edges:
					continue
				# append index to i-th row, j-th column
				row.insert(len(row), i)
				col.insert(len(col), j)
				# -1 for laplacian edges
				data.insert(len(data), -1)

		# graph laplacian
		L = coo_matrix((data,(row, col)), shape=(num_nodes,num_nodes)).tocsc()
		time_T = -0.1
		self.laplacian = L
		self.index2node = index2node
		self.kernel = expm(time_T*L)
		self.labels = node_order
	
		#self.printLaplacian()

	def printLaplacian(self):

		cx = self.laplacian.tocoo()
		for i,j,v in zip(cx.row, cx.col, cx.data):
			a = self.index2node[i]
			b = self.index2node[j]
			print "\t".join([a,b,str(v)])

	def parseNet(self, network):

		edges = set()
		nodes = set()	
		degrees = {}
		for line in open(network, 'r'):

			parts = line.rstrip().split("\t")
			source = parts[0]
			target = parts[2]

			# if inputting a multi-graph, skip this
			if (source, target) in edges:
				continue

			edges.add((source, target))
			edges.add((target, source))
			nodes.add(source)
			nodes.add(target)

			if source not in degrees:
				degrees[source] = 0
			if target not in degrees:
				degrees[target] = 0

			degrees[source] += 1
			degrees[target] += 1

		return (edges, nodes, degrees)


	def kernelMultiplyOne(self, vector):
		"""
			Input:
				vector: A hash mapping gene labels to floating point values 
		"""
		array = []
		# loop over gene names in the network kernel: add the starting value if 
		# it's present in the supplied input vector
		for label in self.labels:
			if label in vector:
				array.append(vector[label])
			else:
				array.append(0)

		# take the dot product
		value = self.kernel*array

		return_vec = {}
		idx = 0
		for label in self.labels:
			return_vec[label] = float(value[idx])
			idx += 1

		return return_vec

	@staticmethod
	def getAngle(v1, v2):

		arry1 = []
		arry2 = []
		for key in v1:
			arry1.append(float(v1[key]))
			arry2.append(float(v2[key]))

		mag_1 = math.sqrt(dot(arry1,arry1))
		mag_2 = math.sqrt(dot(arry2,arry2))

		cos_theta = dot(arry1,arry2)/(mag_1*mag_2)

		return math.acos(cos_theta)

	def diffuse(self, vector, reverse=False):

		# reverse is not used: heat diffusion is undirected
		diffused_vector = self.kernelMultiplyOne(vector)

		return diffused_vector


