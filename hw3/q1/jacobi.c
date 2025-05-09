#include "poisson.h"
#include "jacobi.h"

void sweep2d(double a[][maxn], double f[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y, double b[][maxn]){

	double h;
	h = 1.0/((double)maxn);

	for(int i=s_x; i<=e_x; i++){
		for(int j=s_y; j<=e_y; j++){
			b[i][j] = 0.25 * (a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1]  - h*h*f[i][j]);
		}
	}
}

void Isend_Irecv_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){
	MPI_Request reqs[8];
	MPI_Datatype vect;
	
	int row_lenght = e_x+1 - s_x;
	int col_lenght = e_y+1 - s_y;

	MPI_Irecv(&x[s_x-1][s_y], col_lenght, MPI_DOUBLE, nbrleft , 0, comm, &reqs[0]);
	MPI_Irecv(&x[e_x+1][s_y], col_lenght, MPI_DOUBLE, nbrright, 0, comm, &reqs[1]);
	MPI_Isend(&x[e_x][s_y],   col_lenght, MPI_DOUBLE, nbrright, 0, comm, &reqs[2]);
	MPI_Isend(&x[s_x][s_y],   col_lenght, MPI_DOUBLE, nbrleft , 0, comm, &reqs[3]);
	

	MPI_Type_vector(row_lenght, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);

	MPI_Irecv(&x[s_x][s_y-1], 1, vect, nbrdown, 0, comm, &reqs[4]);
	MPI_Irecv(&x[s_x][e_y+1], 1, vect, nbrup,   0, comm, &reqs[5]);
	MPI_Isend(&x[s_x][e_y],   1, vect, nbrup,   0, comm, &reqs[6]);
	MPI_Isend(&x[s_x][s_y],   1, vect, nbrdown, 0, comm, &reqs[7]);

	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	MPI_Type_free(&vect);
}

double griddiff_2d(double a[][maxn], double b[][maxn], int s_x, int e_x, int s_y, int e_y){
	double sum;
	double tmp;
	int i, j;

	sum = 0.0;
	for(i=s_x; i<=e_x; i++){
		for(j=s_y; j<=e_y; j++){
			tmp = a[i][j] - b[i][j];
			sum = sum + tmp*tmp;
		}
	}

	return sum;
}
