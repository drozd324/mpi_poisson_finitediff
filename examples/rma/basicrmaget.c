#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

int main(int argc, char **argv)
{
  MPI_Comm comm;
  MPI_Win win;
  double sendbuf[] = {0.0,0.0,0.0,0.0,0.0};
  double recvbuf[] = {0.0,0.0,0.0,0.0,0.0};
  int count = 5;
  int myid, nprocs;
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if( myid == 0 ){
    if (nprocs < 2 ){
      fprintf(stderr,"\n=======> Error: this example requires at least two processors\n\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_split(MPI_COMM_WORLD,  myid <= 1, myid, &comm);

  /* init sendbuf, recvbuf */
  if( myid == 0 ){
    for(i=0; i < count; i++)
      sendbuf[i] = (double)i + 1;
  } else if ( myid == 1 ){
    for(i=0; i < count; i++)
      recvbuf[i] = -10000.0;
  }

  if(myid > 1) return 0;
  
  if (myid == 0){
    MPI_Win_create(sendbuf, count*sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);

  } else if (myid == 1){
    /* "getting" process does not need to set aside memory */
    MPI_Win_create(MPI_BOTTOM, 0, sizeof(double), MPI_INFO_NULL, comm, &win);
  }
  
  /* start epoch */
  MPI_Win_fence(0, win);

  if( myid == 1 ){
    MPI_Get(recvbuf, count, MPI_DOUBLE, 0, 0, count, MPI_DOUBLE, win);
  }
  /* do other computations here */

  MPI_Win_fence(0, win);
  /* end epoch */

  MPI_Win_free(&win);

  printf("Result of send and receive: \n");
  for(i=0; i < count; i++){
    printf("(myid: %d) r[%d] = %lf; ",myid, i,recvbuf[i]);
  }
  printf("\n");
  printf("DONE\n");

  MPI_Comm_free(&comm);

  MPI_Finalize();
  return 0;
}
