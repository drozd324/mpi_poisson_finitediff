#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#define MY_BUFSIZE 10

int main(int argc, char *argv[])
{
  MPI_Aint lb, extent;
  MPI_Datatype etype, filetype, contig;
  MPI_Offset disp;
  MPI_File fh;
  int rank, nprocs;
  char fname[200];
  int fname_len;
  int buf[MY_BUFSIZE];
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(rank == 0){
    if( argc != 2 ){
      fprintf(stderr,"Error: usage: %s <output filename>\n",argv[0]);
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    sprintf(fname, "%s", argv[1]);
    fname_len=strlen(fname);
  }
  if(nprocs > 1){
    fprintf(stderr,"This example should only be run on one process\n");
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  MPI_Bcast(&fname_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(fname, (fname_len+1), MPI_CHAR, 0, MPI_COMM_WORLD);
  printf("rank: %d): fname: %s\n",rank,fname);

  for(i=0; i < MY_BUFSIZE; i++){
    buf[i] = i+1;
  }

  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR,
		MPI_INFO_NULL, &fh);

  MPI_Type_contiguous(2, MPI_INT, &contig);
  lb = 0;
  extent = 6*sizeof(int);
  MPI_Type_create_resized(contig, lb, extent, &filetype);
  MPI_Type_commit(&filetype);

  //  disp = 5*sizeof(int);
  disp=0;

  etype = MPI_INT;

  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);

  MPI_File_write(fh, buf, MY_BUFSIZE, MPI_INT, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
  
  MPI_Finalize();
  return 0;
}
