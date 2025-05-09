#include <stdlib.h>
#include <stdio.h>

/* for the sleep command below */
#include <unistd.h>

#include <mpi.h>

int main(int argc, char *argv[])
{
  int MPCW_rank, size, count;
  int send_buf, send_buf2;
  int recv_buf, recv_buf2;

  MPI_Group group_world, grprem;
  MPI_Comm smallcomm;
  int smallcommrank;
  int excluderank[] = {0};

  MPI_Init(&argc, &argv);
  
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  MPI_Comm_rank(MPI_COMM_WORLD, &MPCW_rank); /* local */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* local */

  if( size < 3 ){
    if( MPCW_rank == 0 ){
      fprintf(stderr,"\n=====> Error: example won't work properly on less than 3 procs\n\n");
    }
    MPI_Abort(MPI_COMM_WORLD, 2);
  }

  if( MPCW_rank == 0 ){
    printf("\n====> removing rank %d from MPCW group to make new group\n\n",excluderank[0]);
  }
  MPI_Group_excl(group_world, 1, excluderank, &grprem); /* local */
  
  MPI_Comm_create(MPI_COMM_WORLD, grprem, &smallcomm); /* non local */
  
  send_buf = send_buf2 = MPCW_rank + 10;
  recv_buf = recv_buf2 = -1;
  count = 1;

  printf("(MPCW: %d): send_buf: %d, %d\n",MPCW_rank, send_buf, send_buf2);
    
  if(MPCW_rank != excluderank[0]){

    MPI_Comm_rank(smallcomm, &smallcommrank);
  /* compute on slave */
    /* sleep(2); */
    /* printf("(smallcomm rank: %d (MPCW rank: %d)): finished sleeping\n",smallcommrank, MPCW_rank); */

    MPI_Reduce(&send_buf, &recv_buf, count, MPI_INT, MPI_SUM, 1, smallcomm);

    printf("(smallcomm rank: %d (MPCW rank: %d): recv_buf: %d\n",smallcommrank, MPCW_rank, recv_buf);
    
    MPI_Comm_free(&smallcomm);
  }
  /* zero falls through immediately to this reduce, others do
     later... */
  MPI_Reduce(&send_buf2, &recv_buf2, count, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  printf("(MPCW: %d): recv_buf2: %d\n",MPCW_rank, recv_buf2);
  
  MPI_Group_free(&group_world);
  MPI_Group_free(&grprem);
  
  MPI_Finalize();
  return 0;
}
