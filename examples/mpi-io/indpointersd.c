#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

int main(int argc, char* argv[])
{
  MPI_File fh;
  MPI_Status status;
  char fname[100];
  int fname_len;
  double *data;
  int rank, nprocs, n_per_proc;
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if( argc != 3 ){
    fprintf(stderr, "Usage: %s <number of doubles per processor> <binary file>\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 10);
    exit(1);
  }

  if(rank == 0){
    n_per_proc = atoi(argv[1]);
    sprintf(fname, "%s", argv[2]);
    fname_len=strlen(fname);
  }

  MPI_Bcast(&n_per_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&fname_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(fname, (fname_len+1), MPI_CHAR, 0, MPI_COMM_WORLD);
  printf("rank: %d): fname: %s\n",rank,fname);
  
  data = (double *)malloc(n_per_proc*sizeof(double));

  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

  /* int MPI_File_seek(MPI_File fh, MPI_Offset offset, */
  /*      int whence) */
  MPI_File_seek(fh, rank*n_per_proc*sizeof(double), MPI_SEEK_SET);

       /* int MPI_File_read(MPI_File fh, void *buf, */
       /*      int count, MPI_Datatype datatype, MPI_Status *status) */
  MPI_File_read(fh, data, n_per_proc, MPI_DOUBLE, &status);

  MPI_File_close(&fh);

  for(i=0; i < n_per_proc; i++){
    printf("(rank: %d): %lf\n",rank, data[i]);
  }
  
  free(data);
  MPI_Finalize();

  return 0;

}
