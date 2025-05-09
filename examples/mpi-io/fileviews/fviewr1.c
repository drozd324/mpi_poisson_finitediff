#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

int main(int argc, char* argv[])
{
  MPI_File fh;
  MPI_Status status;
  MPI_Datatype filetype;
  char fname[100];
  int fname_len;
  int *data;
  int rank, nprocs, n_per_blk, nrep;
  int nitems_read;
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if( argc != 4 ){
    fprintf(stderr, "Usage: %s <number of ints per blk> <number of repeats> <binary file>\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 10);
    exit(1);
  }

  if(rank == 0){
    n_per_blk = atoi(argv[1]);
    nrep = atoi(argv[2]);
    sprintf(fname, "%s", argv[3]);
    fname_len=strlen(fname);
  }

  MPI_Bcast(&n_per_blk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nrep, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&fname_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(fname, (fname_len+1), MPI_CHAR, 0, MPI_COMM_WORLD);
  if(rank==0){
    printf("rank: %d): fname: %s\n",rank,fname);
    printf("rank: %d: nperblk: %d; nrep: %d\n",rank,n_per_blk,nrep);
  }

  data = (int *)malloc(nrep*n_per_blk*sizeof(int));

  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

  MPI_Type_vector(nrep, n_per_blk, n_per_blk*nprocs, MPI_INT, &filetype);
  MPI_Type_commit(&filetype);

  MPI_File_set_view(fh, n_per_blk*sizeof(int)*rank, MPI_INT, filetype, "native", MPI_INFO_NULL);

  MPI_File_read_all(fh, data, nrep*n_per_blk, MPI_INT, &status);

  MPI_File_close(&fh);

  MPI_Type_free(&filetype);

  MPI_Get_count(&status, MPI_INT, &nitems_read);
  printf("(rank: %d): items read: %d; \n",rank, nitems_read);

  for(i=0;i<nrep*n_per_blk; i++){
    printf("(rank: %d): data[%d]: %d\n",rank,i, data[i]);
  }
  
  free(data);
  MPI_Finalize();

  return 0;
}
