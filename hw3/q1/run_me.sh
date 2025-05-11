#!/bin/bash
cd "${0%/*}" || exit

SOLVE_SIZE=15
MAXN=$((SOLVE_SIZE + 2))
NUM_PROCS=4

# For each message passing type
for ((i=0; i<3; i++)) do
    # Replace second-last character in line 29 with $i
	sed -i "29s/^\(.*\)\(.\)\(.\)$/\1${i}\3/" main.c

    make
    mpirun -n ${NUM_PROCS} ./poiss2d ${SOLVE_SIZE} > /dev/null 2>&1
done

make clean

python3 plotting.py
