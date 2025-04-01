#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define MAXDIMS 5

int main(int argc, char *argv[])
{
     int rank, size;
     int ndims;
     int dims[2] = {size} ;
     int periods[2] = {0, 0};
     int reorder;
     MPI_Comm cartcomm;
     int source, dest;

     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &size);

     ndims      = 1;
     dims[0]    = size;
     reorder    = 0;

     // makes a new communicator to which the cartesian topology info can be attached
     MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);

     /*
     int MPI_Cart_shift(MPI_Comm comm, int direction, int disp,
			int *rank_source, int *rank_dest)
     */
     MPI_Cart_shift(cartcomm, 0, 1, &source, &dest);

     printf("UPSHIFT(right) (rank %d) source: %d; dest: %d\n",rank, source, dest);

     MPI_Cart_shift(cartcomm, 0, -1, &source, &dest);
     printf("DOWNSHIFT(left) (rank %d) source: %d; dest: %d\n",rank, source, dest);

     MPI_Finalize();
     
     return 0;
}
