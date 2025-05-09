/*
  Example where only _one_ instance of the MPI type is being used;
  __not__ an array of the types
 */

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

struct {
  int a;
  char b;
} mystruct;

int main(int argc, char *argv[])
{
  int blockcounts[2] = {1, 1};
  MPI_Datatype types[2];
  MPI_Aint     displs[2];
  MPI_Datatype structtype;
  int rank, nprocs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if( rank == 0 ){
    if (nprocs < 2 ){
      fprintf(stderr,"\n=======> Error: this example requires at least two processors\n\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  types[0] = MPI_INT;
  types[1] = MPI_CHAR;

  MPI_Get_address(&mystruct.a, &displs[0]);
  MPI_Get_address(&mystruct.b, &displs[1]);

  printf("displs[1]: %ld; displ0: %ld\n", displs[1], displs[0]);

  /* UNLESS using MPI_BOTTOM, make displacements offsets from the
     beginning of the structure */
  displs[1] -= displs[0];
  displs[0] -= displs[0];
  /* displs[1] = MPI_Aint_diff(displs[1], displs[0]); */
  /* displs[0] = MPI_Aint_diff(displs[0], displs[0]); */

  printf("displs[1]: %ld; displ0: %ld\n", displs[1], displs[0]);

  MPI_Type_create_struct(2, blockcounts, displs, types, &structtype);
  MPI_Type_commit(&structtype);

  mystruct.a = 0;
  mystruct.b = 'a';

  if( rank != 0 ){
    printf("BEFORE BCAST: (myid: %d): mystruct.a: %d; mystruct.b: %c\n",rank, mystruct.a, mystruct.b);
  }
  if( rank == 0 ){
      mystruct.a = 11;
      mystruct.b = 'z';
  }

  MPI_Bcast(&mystruct, 1, structtype, 0, MPI_COMM_WORLD);
  /* Using MPI_BOTTOM, we don't have to make displacements relative to
     beginning of mystruct */
  /* MPI_Bcast(MPI_BOTTOM, 1, structtype, 0, MPI_COMM_WORLD); */

  printf("AFTER BCAST (myid: %d): mystruct.a: %d; mystruct.b: %c\n",rank, mystruct.a, mystruct.b);

  MPI_Type_free(&structtype);

  MPI_Finalize();

   return 0;
}
