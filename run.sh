#!/bin/bash

# Run the application
mpirun -np 3 ./build/dynamic_scc < ./tests/input/i2.txt