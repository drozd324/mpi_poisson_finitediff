#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

int main(int argc, char **argv)
{
  MPI_Status status;
  double sendbuf[] = {0.0,0.0,0.0,0.0,0.0};
  double recvbuf[] = {0.0,0.0,0.0,0.0,0.0};
  MPI_Request request;
  int count = 5;
  int nprocs, myid;
  int tag=10;
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if( myid == 0 ){
    if (nprocs !=2 ){
      fprintf(stderr,"\n=======> Error: this example requires exactly two processors\n\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  /* init sendbuf, recvbuf */
  if( myid == 0 ){
    for(i=0; i < count; i++)
      sendbuf[i] = (double)i;
  } else if ( myid == 1 ){
    for(i=0; i < count; i++)
      recvbuf[i] = -10000.0;
  }
  
  if (myid == 0){

    MPI_Isend(sendbuf, count, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &request);

    /**** do some computation to mask latency ****/

    MPI_Wait(&request, &status);

  } else if (myid == 1){
    
    MPI_Irecv(recvbuf, count, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &request);

    /**** do some computation to mask latency ****/

    MPI_Wait(&request, &status);
    /* MPI_Wait deallocates request: */
    if( request == MPI_REQUEST_NULL ){
      printf("YES; recv request deallocated after MPI_Wait\n");
    }
    printf("RECV STATUS: status.MPI_SOURCE: %d; status.MPI_TAG: %d; status.MPI_ERROR: %d\n",status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    printf("Result of send and receive: \n");
  }
  
  for(i=0; i < count; i++){
    printf("(myid: %d) r[%d] = %lf; ",myid, i,recvbuf[i]);
  }
  printf("\n");
  printf("DONE\n");

  MPI_Finalize();
  return 0;
}
