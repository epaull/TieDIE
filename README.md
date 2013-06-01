TieDIE: Tied Diffusion for Subnetwork Discovery. 
========

Current Version 
--------

1.0
	
Authors
--------

Evan O. Paull, Daniel Carlin


Requirements
--------

Python 2.7 and the python numpy module are required to run the tiedie 
executable, when using a pre-computed diffusion kernel. 

Either MATLAB or the python scipy module, version 0.12 or later, is 
required for diffusion kernel computation: the latter is free, though
not as computationally efficient as the MATLAB implementation.

* [python 2.7](http://www.python.org/): all modules.
   * [scipy](http://www.scipy.org/): kernel generation
   * [numpy](http://numpy.scipy.org/)
* [MATLAB](http://www.mathworks.com/products/matlab/): kernel generation

Installation
-------

- Install dependencies
- Download the TieDIE repository to the desired location
- (Optional) Pre-Generate kernel file with MATLAB (bin/makeKernel.sh)
- Run TieDIE/bin/tiedie
	
 
Modules
-------

- **tiedie** Python executable to run the TieDIE algorithm. 
- **makeKernel.sh** Shell script executable that calls MATLAB for diffusion kernel file generation.
- (Auxillary) **span.R** An R-implementation of the Prize Collecting Steiner Tree network formulation, that calls 
the BioNet package. 
