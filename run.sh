#!/bin/bash

# Run the application
mpirun -np 4 ./build/dynamic_scc < ./tests/input/i5.txt