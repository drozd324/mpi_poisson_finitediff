#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "mpi.h"

#define M 5
#define N 10

void print_2d_array(int arr[][N], char *mesg);
void print_in_order(int arr[][N], char *mesg);

int main(int argc, char *argv[])
{
  int arr[M][N];
  int recv1[M];
  int recv2[N];
  MPI_Datatype stridedcol;
  int rank, size;
  MPI_Status stat;
  int count;
  int typesize;
  MPI_Aint typelb, typeextent;
  int i, j;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if( rank == 0 ){
    if (size !=2 ){
      fprintf(stderr,"\n=======> Error: this example requires exactly two processors\n\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  if( rank == 0 ){
    count = 0;
    for(i=0; i < M; i++){
      for(j=0; j < N; j++){
	arr[i][j] = count++;
      }
    }
  }
  print_in_order(arr, "");

  for(i=0; i < M; i++)
    recv1[i] = -1.0;
  for(i=0; i < N; i++)
    recv2[i] = -1.0;

  if( rank == 0 ){
    MPI_Send(&arr[1][0], N, MPI_INT, 1, 0, MPI_COMM_WORLD);
  } else if ( rank == 1 ){
    MPI_Recv(recv2, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
    printf("(rank %d) ",rank);
    for(i = 0;i < N; i++){
      printf("%d ",recv2[i]);
    }
    printf("\n");
  } else {
    fprintf(stderr,"ERROR\n");
    exit(1);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if( rank == 0 )
    printf("\nStrided Data\n");

  MPI_Type_vector(M, 1, N, MPI_INT, &stridedcol);
  MPI_Type_commit(&stridedcol);

  MPI_Type_size(stridedcol, &typesize);
  printf("\n_Size_ of stridedcol: %d\n",typesize);
  /* Extent for this example: sizeof(int) + (M-1)*N*sizeof(int) */
  MPI_Type_get_extent(stridedcol, &typelb, &typeextent);
  /* MPI_Aint has type long int here */
  printf("stridedcol: lb: %ld; extent: %ld\n\n",typelb, typeextent);

  if( rank == 0 ){
    MPI_Send(&arr[0][6], 1, stridedcol, 1, 0, MPI_COMM_WORLD);
  } else if ( rank == 1 ){
    MPI_Recv(recv1, M, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
    printf("(rank %d (strided)) ",rank);
    for(i = 0;i < M; i++){
      printf("%d ",recv1[i]);
    }
    printf("\n");
  } else {
    fprintf(stderr,"ERROR\n");
    exit(1);
  }
  
  MPI_Type_free(&stridedcol);

  MPI_Finalize();

  return 0;
}

void print_2d_array(int arr[][N], char *mesg)
{
  int i, j;

  printf("============= %s =============\n",mesg);
  for(i=0; i < M; i++){
    for(j=0; j < N; j++){
      printf("%2d ",arr[i][j]);
    }
    printf("\n");
  }

}

void print_in_order(int arr[][N], char *mesg)
{
  int rank, size;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  fflush(stdout);usleep(500);
  MPI_Barrier(MPI_COMM_WORLD);
  for(i=0; i < size; i++){
    if( rank == i ){
      printf("---------> rank %d <---------\n",rank);
      print_2d_array(arr, mesg);
      printf("-----------------------------\n");
    }
    fflush(stdout);usleep(500);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
