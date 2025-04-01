#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define MAXDIMS 5

int main(int argc, char *argv[])
{
     int rank, size;
     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &size);

     int ndims = 1;
     int dims[1] = {size};
     int periods[1] = {1};
     int reorder = 0;
     MPI_Comm cartcomm;
     int nbrleft, nbrright;

     // makes a new communicator to which the cartesian topology info can be attached
     MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);

     //int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_nbrleft, int *rank_nbrright)
     MPI_Cart_shift(cartcomm, 0, 1, &nbrleft, &nbrright);
     printf("UPSHIFT(right) (rank %d) nbrleft: %d; nbrright: %d\n",rank, nbrleft, nbrright);

     MPI_Cart_shift(cartcomm, 0, -1, &nbrleft, &nbrright);
     printf("DOWNSHIFT(left) (rank %d) nbrleft: %d; nbrright: %d\n",rank, nbrleft, nbrright);

     MPI_Finalize();
     
     return 0;
}
