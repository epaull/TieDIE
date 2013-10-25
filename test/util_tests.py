#!/usr/bin/env	python

import sys, os
sys.path.append(os.path.dirname(sys.argv[0])+'/../lib')

import doctest
import tiedie_util 
import linkers 

doctest.testmod(tiedie_util)
doctest.testmod(linkers)
