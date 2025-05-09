#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

typedef struct {
  double x, y;
} contigtype;

int main(int argc, char *argv[])
{
  MPI_Datatype ctype2;
  int rank, nprocs;
  MPI_Status status;
  contigtype a[4];
  contigtype b[4];
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if( rank==0 && nprocs != 2 ){
    fprintf(stderr, "\n=======> Error: this examples runs on two processors only\n\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  MPI_Type_contiguous(2, MPI_DOUBLE, &ctype2);
  MPI_Type_commit( &ctype2 );

  if(rank == 0){
    for(i=0; i<4; i++){
      a[i].x = i;
      a[i].y = -i;
    }
    MPI_Send(a, 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    MPI_Send(a, 8, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
  } else {
    int count = 0;
    MPI_Recv(b, 4, ctype2, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, ctype2, &count);
    printf("First Recv: (rank: %d) MPI_Get_count: %d\n",rank, count);
    count = -1;
    MPI_Get_elements(&status, ctype2, &count);
    printf("First Recv: (rank: %d) MPI_Get_elements: %d\n",rank, count);
    
    MPI_Recv(b, 4, ctype2, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, ctype2, &count);
    printf("Second Recv: (rank: %d) MPI_Get_count: %d\n",rank, count);
    count = -1;
    MPI_Get_elements(&status, ctype2, &count);
    printf("Second Recv: (rank: %d) MPI_Get_elements: %d\n",rank, count);
  }

  MPI_Finalize();
  return 0;
}

