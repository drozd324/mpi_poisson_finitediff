#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
     int rank, size;
     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &size);

     int ndims = 2;
     int dims[2] = {0, 0}; 
     int periods[2] = {0, 0};
     int reorder = 0;
     MPI_Comm cartcomm;
     int nbrleft, nbrright, nbrabove, nbrbelow;

     // makes a new communicator to which the cartesian topology info can be attached
     MPI_Dims_create(size, ndims, dims);
     MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);

     /*
     int MPI_Cart_shift(MPI_Comm comm, int direction, int disp,
			int *rank_nbrleft, int *rank_nnbrright)
     */
     MPI_Cart_shift(cartcomm, 0, 1, &nbrleft, &nbrright);
     printf("HORIZONTAL(right) (rank %d) nbrleft: %d; nbrright: %d\n", rank, nbrleft , nbrright);

     MPI_Cart_shift(cartcomm, 1, 1, &nbrbelow, &nbrabove);
     printf("VERTICAL(right) (rank %d) nbrbelow: %d; nbrabove: %d\n",rank, nbrbelow, nbrabove);

	 MPI_Comm_free(&cartcomm);
     MPI_Finalize();
     return 0;
}
