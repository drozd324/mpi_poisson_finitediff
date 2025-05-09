#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define MAXDIMS 5

int main(int argc, char *argv[])
{
     int rank, size;
     int ndims;
     int dims[1];
     int periods[1];
     int reorder;
     MPI_Comm cartcomm;
     int coords[MAXDIMS];
     int i;

     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &size);

     ndims      = 1;
     dims[0]    = size;
     periods[0] = 0;
     reorder    = 0;

     MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);

     MPI_Cart_coords(cartcomm, rank, MAXDIMS, coords);
     
     for(i=0; i < MAXDIMS; i++){
     	  printf("(rank %d) coords[%d]=%d\n",rank,i, coords[i]);
     }

     MPI_Finalize();
     
     return 0;
}
