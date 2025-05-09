#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include "poisson.h" 

void init_basic_bv_2d(double a[][maxn], double b[][maxn], double f[][maxn],
             int nx, int ny, int s_x, int e_x, int s_y, int e_y);

int MPE_Decomp1d(int n, int size, int rank, int *s, int *e);

void init_full_grid(double g[][maxn]);
void init_full_grids(double a[][maxn], double b[][maxn] ,double f[][maxn]);


void print_full_grid(double x[][maxn]);
void print_in_order(double x[][maxn], MPI_Comm comm);
void print_grid_to_file(char *fname, double x[][maxn], int nx, int ny);

double analytic_sol(double x, double y);
void analytic_sol_matrix(double u[][maxn], int s, int e, int nx, int id, int size);
double compute_mse(double a[][maxn], double b[][maxn], int s, int e, int nx);

void write_grid(double grid[][maxn], char* file_name, int nx);
void GatherGrid(double grid[][maxn], double proc_grid[][maxn], int nx, MPI_Comm comm);
void decomp1d(int n, int p, int myid, int* s, int* e);

void GatherGrid2d(double grid[][maxn], double proc_grid[][maxn], int nx, int ny, int* dims, int* coords, MPI_Comm comm);
void analytic_sol_matrix_2d(double u[][maxn], int nx, int ny, int* coords, int* dims);
double compute_mse_2d(double a[][maxn], double b[][maxn], int nx, int ny, int* coords, int* dims);

