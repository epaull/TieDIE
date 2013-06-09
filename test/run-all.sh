#!/bin/sh

echo "Running Kernel Tests..."
./kernel_tests.py
echo "Running General Utility Functional Tests..."
./util_test.py

echo "Running Regression Tests..."
./regression.py
