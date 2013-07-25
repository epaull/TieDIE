#!/bin/sh

py_version=`python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))'`
if [ ! $py_version = 2.7 ]; then
	echo "Error: Python 2.7 not installed, or not set as the default python installation"
	exit 1
fi
	
echo "Running Kernel Tests..."
./kernel_tests.py
echo "Running General Utility Functional Tests..."
./util_tests.py

echo "Running Regression Tests..."
./regression.py
