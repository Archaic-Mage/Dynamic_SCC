#!/bin/bash

# Run the application
mpirun -np 4 ./build/dynamic_scc < ./tests/input/i2.txt