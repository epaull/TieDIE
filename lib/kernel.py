from __future__ import print_function
from __future__ import division
#from __future__ import unicode_literals
from numpy import genfromtxt, dot
import numpy as np
import sys
import math

class Kernel:

    def __init__(self, kernel_files):
        """
            Input:

                kernel_file - a tab-delimited matrix file with both a header
                and first-row labels, in the same order.

            Returns:

                Kernel object.
        """

        # Multiple kernels are supported. Linearity of the diffusion kernel allows
        # this feature.

        # store kernel object
        self.kernels = {}
        # Store the header line for each kernel: will be used for lookup
        self.labels = {}
        # The number of rows and columns for each kernel
        self.ncols = {}
        self.nrows = {}
        # parse each kernel file
        for kernel in kernel_files.split(":"):
            # numpy's genfromtxt format. Relatively memory-intensive
            # FIXME: add option for matlab .mat compressed file format
            self.kernels[kernel] = genfromtxt(kernel,delimiter="\t")[1:,1:]
            self.labels[kernel] = None
            fh = open(kernel,'r')
            # get the header line
            for line in fh:
                self.labels[kernel] = line.rstrip().split("\t")[1:]
                break
            fh.close()

            self.ncols[kernel] = self.kernels[kernel].shape[1]-1
            self.nrows[kernel] = self.kernels[kernel].shape[0]-1

    def getLabels(self):
        """
            Return the set of all node/gene labels used by this kernel object
        """
        all_labels = set()
        for label in self.labels:
            all_labels = all_labels.union(set(self.labels[label]))

        return all_labels

    def kernelMultiplyOne(self, kernel, vector):
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
        for label in self.labels[kernel]:
            # Input heats may not actually be in the network.
            # Check and initialize to zero if not
            if label in vector:
                array.append(vector[label])
            else:
                array.append(0)

        # Matrix mulitply op
        value = dot(self.kernels[kernel], array)

        # Convert back to a hash and return diffused heats
        return_vec = {}
        idx = 0
        for label in self.labels[kernel]:
            return_vec[label] = float(value[idx])
            idx += 1

        return return_vec

    def addVectors(self, vector_list):
        """
        Sum vectors: Add hash / float-valued vectors
        """
        sum = {}

        for vec in vector_list:
            for key in vec:
                val = vec[key]
                if key not in sum:
                    sum[key] = val
                else:
                    sum[key] += val

        return sum

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
        # diffuse separately on each kernel (if more than one), and store.
        return_vectors = []
        for kernel in self.kernels:
            # run diffusion on each constituent kernel
            diffused_vector = self.kernelMultiplyOne(kernel, vector)
            return_vectors.append(diffused_vector)

        return self.addVectors(return_vectors)
