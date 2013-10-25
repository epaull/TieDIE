#!/usr/bin/env python2.7

###
### MakeKernel: Generate a heat diffusion kernel from .sif file
###
###	Authors: 
###
###		Evan Paull (epaull@soe.ucsc.edu)
###
###	Requirements:
###
### 	python 2.7.X
###		numpy 1.7+ 
###		scipy 0.12+ 
###


import os, sys
from collections import defaultdict
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-n","--network",dest="network",action="store",default=None,help="Network in .sif file format")
parser.add_option("-o","--output",dest="output",action="store",type="string",default=None,help="Output file to write")
(opts, args) = parser.parse_args()

# local imports assume the directory structure from github . 
sys.path.append(os.path.dirname(sys.argv[0])+'/../lib')
from kernel_scipy import SciPYKernel

# build the kernel
kernel_diffuser = SciPYKernel(opts.network)
# write to file
kernel_diffuser.writeKernel(opts.output)

