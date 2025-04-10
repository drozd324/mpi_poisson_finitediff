#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "poisson1d.h"
#include "jacobi.h"

void sweep1d(double a[][maxn], double f[][maxn], int nx,
	     int s, int e, double b[][maxn])
{
  double h;
  int i,j;

  h = 1.0/((double)(nx+1));

  for(i=s; i<=e; i++){
    for(j=1; j<nx+1; j++){
      b[i][j] = 0.25 * ( a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1]  - h*h*f[i][j] );
    }
  }
}

void sweep2d(double a[][maxn], double f[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y, double b[][maxn]){
	double h;

	h = 1.0/((double)(nx+1));

	for(int i=s_x; i<=e_x; i++){
		for(int j=s_y; j<=e_y; j++){
			b[i][j] = 0.25 * ( a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1]  - h*h*f[i][j] );
		}
	}
}

/* sendrecv */
//void exchang3(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
//	      int nbrleft, int nbrright)
//{
//
//  MPI_Sendrecv(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, &x[s-1][1], nx, MPI_DOUBLE, nbrleft,
//	       0, comm, MPI_STATUS_IGNORE);
//  MPI_Sendrecv(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, &x[e+1][1], nx, MPI_DOUBLE, nbrright,
//	       1, comm, MPI_STATUS_IGNORE);
//
//}

void exchang3_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){
	MPI_Datatype vect;
	
	int row_lenght = e_x - s_x + 1;
	int col_lenght = e_y - s_y + 1;

	MPI_Sendrecv(&x[s_x][s_y]  , row_lenght, MPI_DOUBLE, nbrleft , 0, 
		     &x[s_x-1][s_y], row_lenght, MPI_DOUBLE, nbrright, 0, 
		     comm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(&x[e_x][s_y]  , row_lenght, MPI_DOUBLE, nbrright, 1,
		     &x[s_x+1][s_y], row_lenght, MPI_DOUBLE, nbrleft , 1,
     		     comm, MPI_STATUS_IGNORE);

	MPI_Type_vector(row_lenght, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);

	MPI_Sendrecv(&x[s_x][s_y]  , 1, vect, nbrup  , 2, 
		     &x[s_x][s_y-1], 1, vect, nbrdown, 2, 
		     comm, MPI_STATUS_IGNORE);

	MPI_Sendrecv(&x[e_x][s_y]  , 1, vect, nbrdown, 3,
		     &x[s_x][s_y+1], 1, vect, nbrup  , 3,
     		     comm, MPI_STATUS_IGNORE);
}

//void exchangi1(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
//	       int nbrleft, int nbrright){
//  MPI_Request reqs[4];
//
//  MPI_Irecv(&x[s-1][1], nx, MPI_DOUBLE, nbrleft , 0, comm, &reqs[0]);
//  MPI_Isend(&x[e][1],   nx, MPI_DOUBLE, nbrright, 0, comm, &reqs[2]);
//  MPI_Irecv(&x[e+1][1], nx, MPI_DOUBLE, nbrright, 0, comm, &reqs[1]);
//  MPI_Isend(&x[s][1],   nx, MPI_DOUBLE, nbrleft , 0, comm, &reqs[3]);
//  /* not doing anything useful here */
//
//  MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
//}

void exchangi2(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){
	MPI_Request reqs[8];
	MPI_Datatype vect;
	
	int row_lenght = e_x - s_x + 1;
	int col_lenght = e_y - s_y + 1;
	
	// x-direction
	MPI_Irecv(&x[s_x-1][s_y], row_lenght, MPI_DOUBLE, nbrleft , 0, comm, &reqs[0]);
	MPI_Isend(&x[e_x][s_y],   row_lenght, MPI_DOUBLE, nbrright, 0, comm, &reqs[2]);
	MPI_Irecv(&x[s_x+1][s_y], row_lenght, MPI_DOUBLE, nbrright, 0, comm, &reqs[1]);
	MPI_Isend(&x[s_x][s_y],   row_lenght, MPI_DOUBLE, nbrleft , 0, comm, &reqs[3]);
	
	// y-direction
	MPI_Type_vector(row_lenght, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);

	MPI_Irecv(&x[s_x][s_y-1], 1, vect, nbrdown, 0, comm, &reqs[4]);
	MPI_Isend(&x[s_x][e_y],   1, vect, nbrup,   0, comm, &reqs[6]);
	MPI_Irecv(&x[s_x][e_y+1], 1, vect, nbrup,   0, comm, &reqs[5]);
	MPI_Isend(&x[s_x][s_y],   1, vect, nbrdown, 0, comm, &reqs[7]);

	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	MPI_Type_free(&vect);
}

double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e)
{
  double sum;
  double tmp;
  int i, j;

  sum = 0.0;

  for(i=s; i<=e; i++){
    for(j=1;j<nx+1;j++){
      tmp = (a[i][j] - b[i][j]);
      sum = sum + tmp*tmp;
    }
  }

  return sum;

}


double griddiff_2d(double a[][maxn], double b[][maxn], int s_x, int e_x, int s_y, int e_y){
  double sum;
  double tmp;
  int i, j;

  sum = 0.0;

  for(i=s_x; i<=e_x; i++){
    for(j=s_y; j<=e_y; j++){
      tmp = (a[i][j] - b[i][j]);
      sum = sum + tmp*tmp;
    }
  }

  return sum;
}
