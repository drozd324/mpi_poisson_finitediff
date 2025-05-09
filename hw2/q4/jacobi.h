#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

void sweep2d(double a[][maxn], double f[][maxn], int nx,
		int s_x, int e_x, int s_y, int e_y, double b[][maxn]);

void Sendrecv_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm);


void Isend_Irecv_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm);

double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e);
double griddiff_2d(double a[][maxn], double b[][maxn], int s_x, int e_x, int s_y, int e_y);
