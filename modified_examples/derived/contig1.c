#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

typedef struct {
  double x, y, z;
  double mass;
} Particle;

void InitParticles(Particle *, const int );

#define MAX_PARTICLES 400
#define MAX_P         128

int main(int argc, char *argv[])
{
  MPI_Datatype particletype;
  Particle particles[MAX_PARTICLES];
  int npart;
  int rank, nprocs;
  int typesize;
  MPI_Aint typeextent, typelb;
  MPI_Request req=MPI_REQUEST_NULL;
  MPI_Status status;
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if( rank==0 && nprocs != 2 ){
    fprintf(stderr, "\n=======> Error: this examples runs on two processors only\n\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  
  if (rank == 0 && argc < 2) { 
    fprintf( stderr, "\n===========> Usage: %s n\n\n", argv[0]);
    MPI_Abort( MPI_COMM_WORLD, 2 );
  }
  MPI_Barrier(MPI_COMM_WORLD);
  npart = atoi(argv[1]);

  MPI_Bcast(&npart, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Type_contiguous( 4, MPI_DOUBLE, &particletype );
  MPI_Type_commit( &particletype );

  MPI_Type_size(particletype, &typesize);
  printf("\n_Size_ of structtype: %d\n",typesize);
  MPI_Type_get_extent(particletype, &typelb, &typeextent);
  /* MPI_Aint has type long int here */
  printf("particletype: lb: %ld; extent: %ld\n\n",typelb, typeextent);
  
  InitParticles(particles, npart);
  
  if( rank == 1 ){
    printf("=============== BEFORE Send/Recv: particles on rank 1 ===============\n");
    for(i=0; i < npart; i++){
      printf("particle: %d (%lf, %lf, %lf, m: %lf)\n",
	     i,particles[i].x, particles[i].y, particles[i].z, particles[i].mass);
    }
    printf("==================== BEFORE list complete ==========================\n");
  }

  if( rank == 0 ){
    MPI_Isend(particles, npart, particletype, 1, 2, MPI_COMM_WORLD, &req);
  } else if( rank == 1 ){
    MPI_Irecv(particles, npart, particletype, 0, 2, MPI_COMM_WORLD, &req);
  } else {
    fprintf(stderr, "ERROR: shouldn't be here\n");
    MPI_Abort( MPI_COMM_WORLD, 4 );
  }

  MPI_Wait(&req, &status);
  if( req != MPI_REQUEST_NULL ){
    fprintf(stderr, "ERROR: req should be MPI_REQUEST_NULL\n");
    MPI_Abort( MPI_COMM_WORLD,  5);
  }

  if( rank == 1 ){
    printf("(myid: %d): STATUS: status.MPI_SOURCE: %d; status.MPI_TAG: %d; status.MPI_ERROR: %d\n",
	   rank,status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);
  }

  if( rank == 1 ){
    printf("=============== particles on rank 1 ===============\n");
    for(i=0; i < npart; i++){
      printf("particle: %d (%lf, %lf, %lf, m: %lf)\n",
	     i,particles[i].x, particles[i].y, particles[i].z, particles[i].mass);
    }
    printf("=============== list complete       ===============\n");
  }

  MPI_Type_free(&particletype);

  MPI_Finalize();

  return 0;
}

void InitParticles(Particle *particles, const int npart)
{
  int rank;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for(i=0; i < npart; i++){
    if( rank == 0 ){
      particles[i].x	  = drand48();
      particles[i].y	  = drand48();
      particles[i].z	  = drand48();
      particles[i].mass = 1.0;
    } else {
      particles[i].x	  = 0.0;
      particles[i].y	  = 0.0;
      particles[i].z	  = 0.0;
      particles[i].mass = -1.0;
    }
  }

}
