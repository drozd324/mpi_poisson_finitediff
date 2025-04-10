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

	MPI_Dims_create(size, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);

	/*
	int MPI_Cart_shift(MPI_Comm comm, int direction, int disp,
			int *rank_nbrleft, int *rank_nnbrright)
	*/

	//x
	MPI_Cart_shift(cartcomm, 1, 1, &nbrleft, &nbrright);
	printf("HORIZONTAL (rank %d) nbrleft: %d; nbrright: %d\n", rank, nbrleft , nbrright);
	//y
	MPI_Cart_shift(cartcomm, 0, 1, &nbrbelow, &nbrabove);
	printf("VERTICAL (rank %d) nbrbelow: %d; nbrabove: %d\n",rank, nbrabove, nbrbelow);

	int coords[2];
	MPI_Cart_coords(cartcomm, rank, 2, coords);  // For current rank	
	//int dims[2], periods[2], coords[2];
	//MPI_Cart_get(cart_comm, 2, dims, periods, coords);

	printf("PROPERTIES (rank %d): coords = (%d, %d): dims = %d, %d \n\n",rank, coords[0], coords[1], dims[0], dims[1]);
		

	MPI_Comm_free(&cartcomm);
	MPI_Finalize();
	return 0;
}
