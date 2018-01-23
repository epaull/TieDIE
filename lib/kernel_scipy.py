from __future__ import print_function
from __future__ import division
#from __future__ import unicode_literals
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

            network_file - a tab-delimited file in .sif network format:
            <source> <interaction> <target>

        Returns:

            Kernel object.

        """

        self.labels = {}
        # The number of rows and columns for each kernel
        self.ncols = {}
        self.nrows = {}

        # parse the network, build indexes
        edges, nodes, node_out_degrees = self.parseNet(network_file)
        num_nodes = len(nodes)
        node_order = list(nodes)
        index2node = {}
        node2index = {}
        for i in range(0, num_nodes):
            index2node[i] = node_order[i]
            node2index[node_order[i]] = i

        # construct the diagonals
        # SCIPY uses row and column indexes to build the matrix
        # row and columns are just indexes: the data column stores
        # the actual entries of the matrix
        row = array('i')
        col = array('i')
        data = array('f')
        # build the diagonals, including the out-degree
        for i in range(0, num_nodes):
            # diag entries: out degree
            degree = 0
            if index2node[i] in node_out_degrees:
                degree = node_out_degrees[index2node[i]]
            # append to the end
            # array object: first argument is the index, the second is the data value
            # append the out-degree to the data array
            data.insert(len(data), degree)
            # build the diagonals
            row.insert(len(row), i)
            col.insert(len(col), i)

        # add off-diagonal edges
        for i in range(0, num_nodes):
            for j in range(0, num_nodes):
                if i == j:
                    continue
                if (index2node[i], index2node[j]) not in edges:
                    continue
                # append index to i-th row, j-th column
                row.insert(len(row), i)
                col.insert(len(col), j)
                # -1 for laplacian: i.e. the negative of the adjacency matrix
                data.insert(len(data), -1)

        # Build the graph laplacian: the CSC matrix provides a sparse matrix format
        # that can be exponentiated efficiently
        L = coo_matrix((data,(row, col)), shape=(num_nodes,num_nodes)).tocsc()
        time_T = -0.1
        self.laplacian = L
        self.index2node = index2node
        # this is the matrix exponentiation calculation.
        # Uses the Pade approximiation for accurate approximation. Computationally expensive.
        # O(n^2), n= # of features, in memory as well.
        self.kernel = expm(time_T*L)
        self.labels = node_order

        #self.printLaplacian()

    def getLabels(self):
        """
            Return the set of all node/gene labels used by this kernel object
        """
        #all_labels = set()
        #for label in self.labels:
        all_labels = set(self.labels)

        return all_labels


    def printLaplacian(self):
        """
        Debug function
        """
        cx = self.laplacian.tocoo()
        for i,j,v in zip(cx.row, cx.col, cx.data):
            a = self.index2node[i]
            b = self.index2node[j]
            print("\t".join([a,b,str(v)]))

    def parseNet(self, network):
        """
        Parse .sif network, using just the first and third columns
        to build an undirected graph. Store the node out-degrees
        in an index while we're at it.
        """
        edges = set()
        nodes = set()
        degrees = {}
        for line in open(network, 'r'):

            parts = line.rstrip().split("\t")
            source = parts[0]
            target = parts[2]

            # if inputing a multi-graph, skip this
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
            Multiply the specified kernel by the supplied input heat vector.

            Input:
                vector: A hash mapping gene labels to floating point values
                kernel: a single index for a specific kernel

            Returns:
                A hash of diffused heats, indexed by the same names as the
                input vector
        """
        # Have to convert to ordered array format for the input vector
        array = []
        for label in self.labels:
            # Input heats may not actually be in the network.
            # Check and initialize to zero if not
            if label in vector:
                array.append(vector[label])
            else:
                array.append(0)

        # take the dot product
        value = self.kernel*array

        # Convert back to a hash and return diffused heats
        return_vec = {}
        idx = 0
        for label in self.labels:
            return_vec[label] = float(value[idx])
            idx += 1

        return return_vec

    def diffuse(self, vector, reverse=False):
        """
        Diffuse input heats over the set of kernels, add to this object

        Input:
            {'gene1': float(heat1)
             'gene2' : float(heat2)
              ...
            }

        Returns:
            Diffused heat vector
        """

        diffused_vector = self.kernelMultiplyOne(vector)

        return diffused_vector
