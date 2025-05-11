#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

void sweep2d(double a[][maxn], double f[][maxn], int nx,
		int s_x, int e_x, int s_y, int e_y, double b[][maxn]);

void Isend_Irecv_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm);

void PutGet_winfence_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm);

void PutGet_activeSync_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm);
