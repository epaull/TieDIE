TieDIE: Tied Diffusion for Subnetwork Discovery. 
========

Current Version 
--------

1.0
	
Authors
--------

Evan O. Paull, Daniel Carlin and Joshua M. Stuart.

Additional Contributors
--------

Srikanth Bezawada (TieDIE Cytoscape Plugin)
Josh L. Espinoza (Quick kernel loading feature)
Dana Silverbush (MATLAB kernel generation code updates to newer versions)

Requirements
--------

Python 2.7 and the python numpy module are required to run the tiedie 
executable, when using a pre-computed diffusion kernel. 

Either MATLAB or the python scipy module, version 0.12 or later, is 
required for diffusion kernel computation: the latter is free, though
not as computationally efficient as the MATLAB implementation.

* [python 2.7](http://www.python.org/): all modules.
   * [scipy](http://www.scipy.org/): >= 0.12.0 kernel generation
   * [numpy](http://numpy.scipy.org/)
   * [networkx](http://networkx.github.io/)
* [MATLAB](http://www.mathworks.com/products/matlab/): kernel generation

Installation
-------

- Install dependencies
- Download the TieDIE repository to the desired location
- (Recommended) Pre-Generate kernel file with MATLAB (bin/makeKernel.sh)
- Run TieDIE/bin/tiedie

Examples
-------
- **GBM.test** An example signaling network from Glioblastoma (TCGA Network, 2012) is provided, along with input heats 
for an upstream set of genes (mutated genes) and a downstream set of nodes (transcriptional responses). To run:

	cd examples/GBM.test
	make

A tutorial for this simple example is provided in doc/Tutorial.pdf. 

Programs
-------

- **tiedie** Python executable to run the TieDIE algorithm. 
- **makeKernel.sh** Shell script executable that calls MATLAB for diffusion kernel file generation.
- (Auxillary) **span.R** An R-implementation of the Prize Collecting Steiner Tree network formulation, that calls 
the BioNet package. 

Folders
------
* bin : executables and matlab source files
* lib : python code libraries for the tiedie executable
* test : doctest unit tests, functional tests and regression tests
* examples : GBM and BRCA inputs for demonstration purposes
* galaxy : Galaxy web-server wrapper for tiedie to run through the web interface. (https://main.g2.bx.psu.edu/)
* pathways : the "superpathway" described in the TieDIE paper, used with the TCGA BRCA dataset

In Press
------
TieDIE was first featured in the 2013 Nature paper "Comprehensive molecular characterization of clear cell renal cell carcinoma". In this TCGA (The Cancer Genome Atlas) network publication, a TieDIE analysis was used to connnect frequently mutated genes involving the SWI/SNF chromatin remodelling complex to a diverse set of gene expression changes characteristic of tumor development and progression. The TieDIE manuscript was not yet published at the time of the Nature publication and so is cited by name and author only. The TieDIE network solution is shown in figure 4 of the main text, which can be found at this link: http://www.nature.com/nature/journal/v499/n7456/full/nature12222.html . 



Contact
------
Feature requests, comments and requests for clarification should all be sent to the author at <epaull@soe.ucsc.edu>. 
I will try to respond quickly to all requests, so feel free to email me!
