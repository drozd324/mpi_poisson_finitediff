#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

struct {
  int a;
  char b;
/* } __attribute__ ((__packed__)) mystruct, struct_array[10]; */
} mystruct, struct_array[10];

/* } __attribute__ ((__packed__)) mystruct, struct_array[10]; */


int main(int argc, char *argv[])
{
  int blockcounts[2] = {1, 1};
  MPI_Datatype types[2];
  MPI_Aint     displs[2];
  MPI_Datatype structtype;
  MPI_Aint typelb, typextent;
  int rank, nprocs;
  int typesize;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MPI_Get_address(&struct_array[0].a, &displs[0]);
  MPI_Get_address(&struct_array[0].b, &displs[1]);

  types[0] = MPI_INT;
  types[1] = MPI_CHAR;

  displs[1] -= displs[0];
  displs[0] = 0;

  MPI_Type_create_struct(2, blockcounts, displs, types, &structtype);
  /* check that struct has correct extent */
  MPI_Type_get_extent(structtype, &typelb, &typextent);
  MPI_Type_size(structtype, &typesize);
  printf("\n_Size_ of old structtype: %d (sizeof(struct_array[0]): %ld; mpi extent: %ld\n",
	 typesize, sizeof(struct_array[0]), typextent);
  if( typextent != sizeof(struct_array[0]) ){
    printf("!! sizeof(struct_array[0]): %ld; mpi typextent: %ld\nCreate new datatype...",
	   sizeof(struct_array[0]), typextent);
    MPI_Datatype sold = structtype;

    /* int MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb,
       MPI_Aint extent, MPI_Datatype *newtype); */
    MPI_Type_create_resized(sold, 0, sizeof(mystruct), &structtype);
    MPI_Type_free(&sold);
  }
  MPI_Type_commit(&structtype);

  MPI_Type_size(structtype, &typesize);
  printf("\nmpi _size_ of new structtype: %d; sizeof(struct_array[0]): %ld\n",
	 typesize, sizeof(struct_array[0]));
  MPI_Type_get_extent(structtype, &typelb, &typextent);
  /* MPI_Aint has type long int here */
  printf("new structtype: lb: %ld; extent: %ld\n\n",typelb, typextent);

  MPI_Type_free(&structtype);

  MPI_Finalize();

   return 0;
}
