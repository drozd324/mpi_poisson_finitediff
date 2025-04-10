#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <mpi.h>

#include "poisson1d.h"
#include "jacobi.h"

#define maxit 2000

#include "decomp1d.h"

void init_full_grid(double g[][maxn]);
void init_full_grids(double a[][maxn], double b[][maxn] ,double f[][maxn]);

void init_basic_bv(double a[][maxn], double b[][maxn], double f[][maxn],
			 int nx, int ny, int s, int e);

void print_full_grid(double x[][maxn]);
void print_in_order(double x[][maxn], MPI_Comm comm);
void print_grid_to_file(char *fname, double x[][maxn], int nx, int ny);

//q2
double analytic_sol(double x, double y);
void analytic_sol_matrix(double u[][maxn], int s, int e, int nx, int id, int size);
double compute_mse(double a[][maxn], double b[][maxn], int s, int e, int nx);

//q3
void write_grid(double grid[][maxn], char* file_name, int nx);
void GatherGrid(double grid[][maxn], double proc_grid[][maxn], int nx, MPI_Comm comm);
void decomp1d(int n, int p, int myid, int* s, int* e);

//q4
void GatherGrid2d(double grid[][maxn], double proc_grid[][maxn], int nx, int ny, int* dims, int* coords, MPI_Comm comm);
void init_basic_bv_2d(double a[][maxn], double b[][maxn], double f[][maxn], int nx, int ny, int s_x, int e_x, int s_y, int e_y);
void analytic_sol_matrix_2d(double u[][maxn], int nx, int ny, int* coords, int* dims);
double compute_mse_2d(double a[][maxn], double b[][maxn], int nx, int ny, int* coords, int* dims);

int main(int argc, char **argv)
{
	double a[maxn][maxn], b[maxn][maxn], f[maxn][maxn];
	double u[maxn][maxn];
	int nx, ny;
	int myid, nprocs;
	/* MPI_Status status; */
	int nbrleft, nbrright, nbrup, nbrdown, s_x, e_x, s_y, e_y, it;
	double glob_diff;
	double ldiff;
	double t1, t2;
	double tol=1.0E-11;
	char name[1024];
	int namelen;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	MPI_Get_processor_name(name, &namelen);
	printf("(myid %d): running on node: %s\n",myid, name);

	if( myid == 0 ){
		/* set the size of the problem */
		if(argc > 2){
			fprintf(stderr,"---->Usage: mpirun -np <nproc> %s <nx>\n",argv[0]);
			fprintf(stderr,"---->(for this code nx=ny)\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if( argc == 2 ){
			nx = atoi(argv[1]);
		}

		if( nx > maxn-2 ){
			fprintf(stderr,"grid size too large\n");
			exit(1);
		}
	}

	MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	printf("(myid: %d) nx = %d\n",myid,nx);
	ny = nx;

	init_full_grids(a, b, f);
	
	//////////////////////////////////////// q4 ///////////////////////////////////////////////////////////
	int ndims = 2;
        int dims[2] = {0, 0};
	int periods[2] = {0, 0};
	int reorder = 0;
	MPI_Comm cartcomm;
	int coords[2];

	MPI_Dims_create(nprocs, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);
	MPI_Cart_coords(cartcomm, myid, 2, coords);	

	//x
        MPI_Cart_shift(cartcomm, 1, 1, &nbrleft, &nbrright);
	//y                                
	MPI_Cart_shift(cartcomm, 0, 1, &nbrdown, &nbrup);

	for (int i=0; i<2; i++){
		if (coords[i] == -1){
			coords[i] = MPI_PROC_NULL;
		}
	}

	//////////////////////////////////////// q4 ///////////////////////////////////////////////////////////
	MPE_Decomp1d(nx, dims[1], coords[1], &s_x, &e_x);
	MPE_Decomp1d(ny, dims[0], coords[0], &s_y, &e_y);

	printf("(myid: %d) nx=ny: %d; s_x: %d; e_x: %d; s_y: %d; e_y: %d; nbrleft: %d; nbrright: %d; nbrup: %d; nbrdown: %d\n",
		myid, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown);

	init_basic_bv_2d(a, b, f, nx, ny, s_x, e_x, s_y, e_y);

	print_in_order(a, MPI_COMM_WORLD);

	/* MPI_Barrier(MPI_COMM_WORLD); */
	/* MPI_Abort(MPI_COMM_WORLD, 1); */

	t1 = MPI_Wtime();

	/////////////////////////////////////////////////////////////////////////////////////////////////
	glob_diff = 1000;
	for(it=0; it<maxit; it++){

		exchangi2(a, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
		sweep2d(a, f, nx, s_x, e_x, s_y, e_y, b);

		exchangi2(b, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
		sweep2d(b, f, nx, s_x, e_x, s_y, e_y, a);

		ldiff = griddiff_2d(a, b, s_x, e_x, s_y, e_y);
		MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if(myid==0 && it%10==0){
			printf("(myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff, glob_diff);
		}
		/* if(it%5==0){ */
		/*	 print_in_order(a, MPI_COMM_WORLD); */
		/* } */
		if( glob_diff < tol ){
			if(myid==0){
		printf("iterative solve converged\n");
			}
			break;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	t2=MPI_Wtime();
	
	printf("DONE! (it: %d)\n",it);

	if( myid == 0 ){
		if( it == maxit ){
			fprintf(stderr,"Failed to converge\n");
		}
		printf("Run took %lf s\n",t2-t1);
	}

	print_in_order(a, MPI_COMM_WORLD);
	if( nprocs == 1	){
		print_grid_to_file("grid", a,	nx, ny);
		print_full_grid(a);
	}
	
	// mean squared error calcluation for part 2 updated to part 4	
	analytic_sol_matrix_2d(u, nx, ny, coords, dims);
	double whole_u[maxn][maxn];
	GatherGrid2d(whole_u, u, nx, ny, dims, coords, MPI_COMM_WORLD);
	if (myid == 0){
		write_grid(whole_u, "./grids/analytic_grid.txt", nx);
	}
	//print_in_order(u, MPI_COMM_WORLD);
	double local_mse = compute_mse_2d(a, u, nx, ny, coords, dims);
	double global_mse = 0;
	MPI_Reduce(&local_mse, &global_mse, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	global_mse = global_mse / 3;
	if (myid == 0){
		printf("\nMSE compared with analytic solution = %.15lf\n\n", global_mse);
	}

	MPI_Barrier(MPI_COMM_WORLD);
		
	// q3
	double whole_grid[maxn][maxn];
	GatherGrid2d(whole_grid, a, nx, ny, dims, coords, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 0){
		printf("Priting whole grid on rank = 0\n");
		print_full_grid(whole_grid);
		write_grid(whole_grid, "./grids/grid.txt", nx);
	}

	MPI_Finalize();
	return 0;
}


void decomp1d(int n, int p, int myid, int* s, int* e){
	int remainder = n % p;
	int base = n / p;

	if (myid < remainder) {
		*s = myid * (base + 1);
		*e = *s + (base + 1);
	} else {
		*s = remainder * (base + 1) + (myid - remainder) * base;
		*e = *s + base;
	}
}

// writes grid to file
void write_grid(double grid[][maxn], char* file_name, int nx){
	FILE *fp = fopen(file_name, "w");

	if (fp == NULL) {
		printf("Error opening file!\n");
	}
	
	for (int i=0; i<nx+2; i++){
		for (int j=0; j<nx+2; j++){
			fprintf(fp, "%lf ", grid[i][j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void GatherGrid(double grid[][maxn], double proc_grid[][maxn], int nx, MPI_Comm comm){
	int myid, size;
	MPI_Comm_rank(comm, &myid);
	MPI_Comm_size(comm, &size);
	MPI_Barrier(comm);
	MPI_Status status;
	
	if (myid != 0){	
		MPI_Send(proc_grid, (nx+2)*(nx+2), MPI_DOUBLE, 0, 100 + myid, comm);
	}

	if (myid == 0){
		// combinting on proc 0
		int s1;
		int e1;
		decomp1d(nx+2, size, myid, &s1, &e1);
		for (int i=s1; i<e1; i++){
			for (int j=0; j<nx+2; j++){
				grid[i][j] = proc_grid[i][j];
			}
		}	

		// collecting from other procs
		double temp_grid[maxn][maxn];
		for (int k=1; k<size; k++){
			init_full_grid(temp_grid);
			MPI_Recv(temp_grid, (nx+2)*(nx+2), MPI_DOUBLE, k, 100 + k, comm, &status);
			decomp1d(nx+2, size, k, &s1, &e1);
			for (int i=s1; i<=e1; i++){
				for (int j=0; j<nx+2; j++){
					grid[i][j] = temp_grid[i][j];
				}
			}	
		}
	}
}

void GatherGrid2d(double grid[][maxn], double proc_grid[][maxn], int nx, int ny, int* dims, int* coords, MPI_Comm comm){
	int myid, size;
	MPI_Comm_rank(comm, &myid);
	MPI_Comm_size(comm, &size);
	MPI_Barrier(comm);
	MPI_Status status;
	
	if (myid != 0){	
		MPI_Send(proc_grid, (nx+2)*(nx+2), MPI_DOUBLE, 0, 100 + myid, comm);
	}

	if (myid == 0){
		int s1_x, s1_y, e1_x, e1_y;

		// combinting on proc 0
		decomp1d(nx+2, dims[1], coords[1], &s1_x, &e1_x);
		decomp1d(ny+2, dims[0], coords[0], &s1_y, &e1_y);
		for (int i=s1_x; i<e1_x; i++){
			for (int j=s1_y; j<e1_y; j++){
				grid[i][j] = proc_grid[i][j];
			}
		}	

		// collecting from other procs
		double temp_grid[maxn][maxn];
		for (int k=1; k<size; k++){
			init_full_grid(temp_grid);
			MPI_Recv(temp_grid, (nx+2)*(nx+2), MPI_DOUBLE, k, 100 + k, comm, &status);
			decomp1d(nx+2, dims[1], coords[1], &s1_x, &e1_x);
			decomp1d(ny+2, dims[0], coords[0], &s1_y, &e1_y);
			for (int i=s1_x; i<e1_x; i++){
				for (int j=s1_y; j<e1_y; j++){
					grid[i][j] = temp_grid[i][j];
				}
			}	
		}
	}
}

double analytic_sol(double x, double y){
	return y / ( (1+x)*(1+x) + y*y);
} 

void analytic_sol_matrix_2d(double u[][maxn], int nx, int ny, int* coords, int* dims){
	int s1_x, e1_x;
	decomp1d(nx+2, dims[1], coords[1], &s1_x, &e1_x);
	int s1_y, e1_y;
	decomp1d(ny+2, dims[0], coords[0], &s1_y, &e1_y);
	
	for (int i=s1_x; i<e1_x; i++){
		for (int j=s1_y; j<e1_y; j++){
			u[i][j] = analytic_sol((double)i/(double)nx, (double)j/(double)ny);
		}
	}
}

double compute_mse_2d(double a[][maxn], double b[][maxn], int nx, int ny, int* coords, int* dims){
	double sum = 0;
	int s1_x, e1_x;
	decomp1d(nx+2, dims[1], coords[1], &s1_x, &e1_x);
	int s1_y, e1_y;
	decomp1d(ny+2, dims[0], coords[0], &s1_y, &e1_y);

	for (int i=s1_x; i<e1_x; i++){
		for (int j=s1_y; j<e1_y; j++){
			sum += (a[i][j] - b[i][j]) * (a[i][j] - b[i][j]);
		}
	}
	return sum / (((double)e1_y - (double)s1_y + 1) * ((double)e1_x - (double)s1_x + 1));
}


void init_basic_bv(double a[][maxn], double b[][maxn], double f[][maxn],
			 int nx, int ny, int s, int e)
{
	int i,j;
	double left, bottom, right, top;
	double u;

	/* set everything to 0 first */
	for(i=s-1; i<=e+1; i++){
		for(j=0; j <= nx+1; j++){
			a[i][j] = 0.0;
			b[i][j] = 0.0;
			f[i][j] = 0.0;
		}
	}

	/* deal with boundaries */
	for(i=s; i<=e; i++){
		a[i][0] = 0; //bottom = 0
		b[i][0] = 0;
		u = 1 / ( (1 + ((double)i/(double)nx))*(1 + ((double)i/(double)nx)) + 1); // 1/((1+x)*(1+x) + 1)
		a[i][nx+1] = u;
		b[i][nx+1] = u;
	}

	/* this is true for proc 0 */
	if( s == 1 ){
		for(j=1; j<nx+1; j++){
			u = ((double)j/(double)nx)/(1 + ((double)j/(double)nx)*((double)j/(double)nx));  // y/(1 + y*y)
			a[0][j] = u;	
			b[0][j] = u;
		}
	}
 
	/* this is true for proc size-1 */
	if( e == nx ){
		for(j=1; j<nx+1; j++){
			u = ((double)j/(double)nx)/(4 + ((double)j/(double)nx)*((double)j/(double)nx));  // y/(4 + y*y)
			a[nx+1][j] = u;
			b[nx+1][j] = u;
		}

	}

}

void init_basic_bv_2d(double a[][maxn], double b[][maxn], double f[][maxn],
			 int nx, int ny, int s_x, int e_x, int s_y, int e_y)
{
	int i,j;
	double left, bottom, right, top;
	double u, w;

	/* set everything to 0 first */
	for(i=s_x-1; i<=e_x+1; i++){
		for(j=s_y-1; j<=e_y+1; j++){
			a[i][j] = 0.0;
			b[i][j] = 0.0;
			f[i][j] = 0.0;
		}
	}

	/* deal with boundaries */
	for(i=s_x; i<=e_x; i++){
		a[i][0] = 0; //bottom = 0
		b[i][0] = 0;

		w = 1 / ( (1 + ((double)i/(double)nx))*(1 + ((double)i/(double)nx)) + 1); // 1/((1+x)*(1+x) + 1)
		a[i][nx+1] = w;
		b[i][nx+1] = w;
	}

	for(j=s_y; j<=e_y; j++){
		u = ((double)j/(double)nx)/(1 + ((double)j/(double)nx)*((double)j/(double)nx));  // y/(1 + y*y)
		a[0][j] = u;	
		b[0][j] = u;

		w = ((double)j/(double)nx)/(4 + ((double)j/(double)nx)*((double)j/(double)nx));  // y/(4 + y*y)
		a[nx+1][j] = w;
		b[nx+1][j] = w;
	}
}

void init_full_grid(double g[][maxn])

{
	int i,j;
	const double junkval = -5;

	for(i=0; i < maxn; i++){
		for(j=0; j<maxn; j++){
			g[i][j] = junkval;
		}
	}
}

/* set global a,b,f to initial arbitrarily chosen junk value */
void init_full_grids(double a[][maxn], double b[][maxn] ,double f[][maxn])
{
	int i,j;
	const double junkval = -5;

	for(i=0; i < maxn; i++){
		for(j=0; j<maxn; j++){
			a[i][j] = junkval;
			b[i][j] = junkval;
			f[i][j] = junkval;
		}
	}

}

/* prints to stdout in GRID view */
void print_full_grid(double x[][maxn])
{
	int i,j;
	for(j=maxn-1; j>=0; j--){
		for(i=0; i<maxn; i++){
			if(x[i][j] < 10000.0){
	printf("|%2.6lf| ",x[i][j]);
			} else {
	printf("%9.2lf ",x[i][j]);
			}
		}
		printf("\n");
	}

}


void print_in_order(double x[][maxn], MPI_Comm comm)
{
	int myid, size;
	int i;

	MPI_Comm_rank(comm, &myid);
	MPI_Comm_size(comm, &size);
	MPI_Barrier(comm);
	printf("Attempting to print in order\n");
	sleep(1);
	MPI_Barrier(comm);

	for(i=0; i<size; i++){
		if( i == myid ){
			printf("proc %d\n",myid);
			print_full_grid(x);
		}
		fflush(stdout);
		usleep(500);	
		MPI_Barrier(comm);
	}
}

void print_grid_to_file(char *fname, double x[][maxn], int nx, int ny)
{
	FILE *fp;
	int i,j;

	fp = fopen(fname, "w");
	if( !fp ){
		fprintf(stderr, "Error: can't open file %s\n",fname);
		exit(4);
	}

	for(j=ny+1; j>=0; j--){
		for(i=0; i<nx+2; i++){
			fprintf(fp, "%lf ",x[i][j]);
			}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

