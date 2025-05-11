#!/bin/bash

SOLVE_SIZE=15
MAXN=$(SOLVE_SIZE)+2

NUM_PROCS=4

# for each message passing type
for ((i=1; i<4; i++)) do
	sed -i "29s/3/$(i)/g" main.c 
	
	make 
	mpirun -n $(NUM_PROCS) ./poiss2d $(SOLVE_SIZE)
done

make clean
