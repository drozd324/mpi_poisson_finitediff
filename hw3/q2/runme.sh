#!/bin/bash
cd "${0%/*}" || exit

make

echo
echo "Running mpi exercise"
mpirun -n 4 ./main

echo
echo "Running python test"
python3 test.py 

echo
make clean

