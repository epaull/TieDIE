#!/usr/bin/env  python2.7

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

    @staticmethod
    def getAngle(v1, v2):
        """
        Inactive Module: Get the angle between two vectors in n-space. 
        Could be used for additional null model test & distance function.
        """
        arry1 = []
        arry2 = []
        for key in v1:
            arry1.append(float(v1[key]))
            arry2.append(float(v2[key]))

        norm_arry1 = [ a/sum(arry1) for a in arry1 ]
        norm_arry2 = [ a/sum(arry2) for a in arry2 ]

        mag_1 = math.sqrt(dot(norm_arry1,norm_arry1))
        mag_2 = math.sqrt(dot(norm_arry2,norm_arry2))

        cos_theta = dot(norm_arry1,norm_arry2)/(mag_1*mag_2)

        return math.acos(cos_theta)

    @staticmethod
    def getSYMKLDiv(v1, v2):
        """
        <Development module> Get the symmetric Kullback-Leibler divergence between input vectors.
        """
        return (Kernel.getKLDiv(v1, v2) + Kernel.getKLDiv(v2, v1) )/ 2

    @staticmethod
    def getKLDiv(v1, v2):
        """
        <Development module> Get the Kullback-Leibler divergence between input vectors, 
        which typically represent diffused heat values. 
        Input:
            2 diffused heat vectors
        Output:
            (float) The KL-divergence metric between vectors
        """

        # convert values to ordered array
        arry1 = []
        arry2 = []
        for key in v1:
            arry1.append(float(v1[key]))
            arry2.append(float(v2[key]))

        # normalize both vectors, adding an epsilon (pseudo-count) value for zero-values
        EPSILON = 0.00001
        norm_arry1 = [ a+EPSILON/(sum(arry1)+len(arry1)*EPSILON) for a in arry1 ]
        norm_arry2 = [ a+EPSILON/(sum(arry2)+len(arry2)*EPSILON) for a in arry2 ]

        # calculate the KL-Div
        div_sum = 0 
        for i in range(0, len(norm_arry1)):
            div_sum += np.log( (norm_arry1[i]/norm_arry2[i]) )*norm_arry1[i]

        return div_sum


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
        return_vectors = []
        for kernel in self.kernels:
            diffused_vector = self.kernelMultiplyOne(kernel, vector)
            return_vectors.append(diffused_vector)

        return self.addVectors(return_vectors)


