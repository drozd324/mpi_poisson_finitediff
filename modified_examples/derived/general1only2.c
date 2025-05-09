/*
  Example where only _one_ instance of the MPI type is being used;
  __not__ an array of the types
 */

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

struct {
  int a[4];
  char b;
} mystruct;

int main(int argc, char *argv[])
{
  int blockcounts[2] = {4, 1};
  MPI_Datatype types[2];
  MPI_Aint     displs[2];
  MPI_Datatype structtype;
  int rank, nprocs;

  int typesize;
  MPI_Aint typelb, typeextent;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MPI_Get_address(&mystruct.a, &displs[0]);
  MPI_Get_address(&mystruct.b, &displs[1]);

  types[0] = MPI_INT;
  types[1] = MPI_CHAR;

  displs[1] -= displs[0];
  displs[0] = 0;

  MPI_Type_create_struct(2, blockcounts, displs, types, &structtype);
  MPI_Type_commit(&structtype);

  MPI_Type_size(structtype, &typesize);
  printf("\n_Size_ of structtype: %d\n",typesize);

  MPI_Type_get_extent(structtype, &typelb, &typeextent);
  /* MPI_Aint has type long int here */
  printf("structtype: lb: %ld; extent: %ld\n\n",typelb, typeextent);
  printf("C sizeof structtype: %ld\n\n",sizeof(mystruct));

  MPI_Type_free(&structtype);

  MPI_Finalize();

   return 0;
}
