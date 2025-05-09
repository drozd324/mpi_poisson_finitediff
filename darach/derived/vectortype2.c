#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "mpi.h"

#define M 5
#define N 10

void print_2d_array(int arr[][N], char *mesg);
void print_in_order(int arr[][N], char *mesg);
void init_arr(int arr[][N]);

int main(int argc, char *argv[])
{
  int arr[M][N];
  MPI_Datatype stridedcol;
  int rank, size;
  MPI_Status stat;
  int count;
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

  init_arr(arr);

  if( rank == 0 ){
    count = 0;
    for(i=0; i < M; i++){
      for(j=0; j < N; j++){
	arr[i][j] = count++;
      }
    }
  }
  print_in_order(arr, "row(BEFORE)");

  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* MPI_Finalize(); */
  /* return 1; */

  if( rank == 0 ){
    MPI_Send(&arr[1][0], N, MPI_INT, 1, 0, MPI_COMM_WORLD);
  } else if ( rank == 1 ){
    MPI_Recv(arr[1], N, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
    /* printf("(rank %d) ",rank); */
    /* for(i = 0;i < N; i++){ */
    /*   printf("%d ",recv2[i]); */
    /* } */
    /* printf("\n"); */
  } else {
    fprintf(stderr,"ERROR\n");
    exit(1);
  }
  print_in_order(arr, "row(AFTER)");
  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* MPI_Finalize(); */
  /* return 1; */

  if( rank == 0 )
    printf("\nStrided Data\n");

  MPI_Type_vector(M, 1, N, MPI_INT, &stridedcol);
  MPI_Type_commit(&stridedcol);

  if( rank == 0 ){
    MPI_Send(&arr[0][1], 1, stridedcol, 1, 0, MPI_COMM_WORLD);
  } else if ( rank == 1 ){
    MPI_Recv(&arr[0][5], 1, stridedcol, 0, 0, MPI_COMM_WORLD, &stat);
  } else {
    fprintf(stderr,"ERROR\n");
    exit(1);
  }

  print_in_order(arr, "col(AFTER)");
  
  MPI_Type_free(&stridedcol);

  MPI_Finalize();

  return 0;
}

void init_arr(int arr[][N])
{
  int i, j;

  for(i = 0;i < M; i++){
    for(j = 0;j < N; j++){
      arr[i][j] = -1;
    }
  }
  
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
