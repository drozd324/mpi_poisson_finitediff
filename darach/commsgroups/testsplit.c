#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>

#include "mpi.h"

int main(int argc, char *argv[])
{
  int world_size, myid, nrank, nsize;
  int colour;
  MPI_Comm world, newcomm;


  MPI_Init(&argc, &argv);
  world = MPI_COMM_WORLD;
  MPI_Comm_rank(world, &myid);
  MPI_Comm_size(world, &world_size);

  colour = 0;
  if(myid%2==1){
    colour = 1;
  }

  MPI_Comm_split(world, colour, 0, &newcomm);
  MPI_Comm_rank(newcomm, &nrank);
  MPI_Comm_size(newcomm, &nsize);

  printf("(world rank, size:  %d, %d): newcomm size: %d; newcomm rank: %d\n",myid, world_size, nsize, nrank);

  MPI_Barrier(world);
  printf("We're all done\n");
  
  MPI_Finalize();

  return 0;
}


