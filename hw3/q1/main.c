#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include "poisson.h"
#include "jacobi.h"
#include "tools.h"

#define maxit 2000

int main(int argc, char **argv){
	double a[maxn][maxn], b[maxn][maxn], f[maxn][maxn];
	double u[maxn][maxn];
	int nx, ny;
	int myid, nprocs;
	/* MPI_Status status; */
	int nbrleft, nbrright, nbrup, nbrdown;
	int s_x, e_x, s_y, e_y;
	int it;
	double glob_diff;
	double ldiff;
	double t1, t2;
	double tol = 1e-11;
	char name[1024];
	int namelen;

	// picks out which message passing functions to use
	int mode = 2;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	MPI_Get_processor_name(name, &namelen);
	printf("(myid %d) running on node: %s\n in mode: %d",myid, name, mode);

	if (myid == 0){
		/* set the size of the problem */
		if (argc > 2){
			fprintf(stderr,"---->Usage: mpirun -np <nproc> %s <nx>\n",argv[0]);
			fprintf(stderr,"---->(for this code nx=ny)\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (argc == 2){
			nx = atoi(argv[1]);
		}

		if (nx > maxn-2){
			fprintf(stderr,"grid size too large\n");
			exit(1);
		}
	}
	
	MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	printf("(myid: %d) nx = %d\n",myid,nx);
	ny = nx;

	// inserts junk values -5
	init_full_grids(a, b, f);
	
	// init 2d grid with cart create
	int ndims = 2;
    int dims[2] = {0, 0};
	int periods[2] = {0, 0};
	int reorder = 0;
	MPI_Comm cartcomm;
	int coords[2];

	MPI_Dims_create(nprocs, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);
	MPI_Cart_coords(cartcomm, myid, 2, coords);	

	MPI_Cart_shift(cartcomm, 0, 1, &nbrdown, &nbrup);
	MPI_Cart_shift(cartcomm, 1, 1, &nbrleft, &nbrright);

	MPE_Decomp1d(nx, dims[0], coords[0], &s_x, &e_x);
	MPE_Decomp1d(ny, dims[1], coords[1], &s_y, &e_y);

	printf("(myid: %d) nx: %d; s_x: %d; e_x: %d; s_y: %d; e_y: %d; nbrleft: %d; nbrright: %d; nbrup: %d; nbrdown: %d\n",
		myid, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown);

	init_basic_bv_2d(a, b, f, nx, ny, s_x, e_x, s_y, e_y);

	t1 = MPI_Wtime();
	glob_diff = 1000;
	
	switch (mode) {
		case 0:
			for (it=0; it<maxit; it++){
				Isend_Irecv_2d(a, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
				sweep2d(a, f, nx, s_x, e_x, s_y, e_y, b);
			
				Isend_Irecv_2d(b, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
				sweep2d(b, f, nx, s_x, e_x, s_y, e_y, a);

				ldiff = griddiff_2d(a, b, s_x, e_x, s_y, e_y);
				MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				if(myid==0 && it%10==0){
					printf("(myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff, glob_diff);
				}
				if(glob_diff < tol ){
					if(myid==0){
						printf("iterative solve converged\n");
					}
					break;
				}
			}
			break;
		case 1:
			for (it=0; it<maxit; it++){
				PutGet_winfence_2d(a, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
				sweep2d(a, f, nx, s_x, e_x, s_y, e_y, b);
			
				PutGet_winfence_2d(b, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
				sweep2d(b, f, nx, s_x, e_x, s_y, e_y, a);

				//if (check_covergence(a, b, s_x, e_x, s_y, e_y, ldiff, glob_diff, tol, it)) break;
				ldiff = griddiff_2d(a, b, s_x, e_x, s_y, e_y);
				MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				if(myid==0 && it%10==0){
					printf("(myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff, glob_diff);
				}
				if(glob_diff < tol ){
					if(myid==0){
						printf("iterative solve converged\n");
					}
					break;
				}
			}
			break;
		case 2:
			for (it=0; it<maxit; it++){
				PutGet_activeSync_2d(a, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
				sweep2d(a, f, nx, s_x, e_x, s_y, e_y, b);
			
				PutGet_activeSync_2d(b, nx, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, MPI_COMM_WORLD);
				sweep2d(b, f, nx, s_x, e_x, s_y, e_y, a);

				ldiff = griddiff_2d(a, b, s_x, e_x, s_y, e_y);
				MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				if(myid==0 && it%10==0){
					printf("(myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff, glob_diff);
				}
				if(glob_diff < tol ){
					if(myid==0){
						printf("iterative solve converged\n");
					}
					break;
				}
			}
			break;
	}

	t2=MPI_Wtime();
    printf("DONE! (it: %d)\n",it);
    if (myid == 0){
        if (it == maxit){
            fprintf(stderr,"Failed to converge\n");
        }
        printf("Run took %lf s\n",t2-t1);
    }
    print_in_order(a, MPI_COMM_WORLD);
    if (nprocs == 1){
        write_grid(a, "grid", nx);
        print_full_grid(a);
    }

    double whole_grid[maxn][maxn];
    GatherGrid2d(whole_grid, a, nx, ny, dims, coords, cartcomm);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0){
        printf("\nPriting whole grid on rank = 0\n");
        print_full_grid(whole_grid);

		char save_filename[32];
		snprintf(save_filename, sizeof(save_filename), "./grids/grid%d.txt", mode);
        write_grid(whole_grid, save_filename, nx);
    }

    // analytic solution
    if (myid == 0){
    	analytic_sol_matrix_2d(u, nx, ny);
        write_grid(u, "./grids/analytic_grid.txt", nx);
        printf("\nPriting analytic solution on rank = 0\n");
        print_full_grid(u);

		double max_abs_diff = max_diff(whole_grid, u, nx, ny, coords, dims);
		printf("MAX DIFF: %lf\n", max_abs_diff);
    }


	MPI_Finalize();
	return 0;
}


