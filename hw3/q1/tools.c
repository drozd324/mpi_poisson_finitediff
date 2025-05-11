#include "tools.h"

void init_full_grids(double a[][maxn], double b[][maxn] ,double f[][maxn]){
    int i, j;
    double junkval = 9.999999;

    for(i=0; i < maxn; i++){
        for(j=0; j<maxn; j++){
            a[i][j] = junkval;
            b[i][j] = junkval;
            f[i][j] = junkval;
        }
    }
}

void init_full_grid(double g[][maxn]){
    int i,j;
    double junkval = 9.999999;

    for(i=0; i<maxn; i++){
        for(j=0; j<maxn; j++){
            g[i][j] = junkval;
        }
    }
}

void init_basic_bv_2d(double a[][maxn], double b[][maxn], double f[][maxn],
             int nx, int ny, int s_x, int e_x, int s_y, int e_y){
    double u, w;
	double x, y;
	int i, j;

    /* set everything to 0 first */
    for (i=s_x; i<=e_x; i++){
        for (j=s_y; j<=e_y; j++){
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            f[i][j] = 0.0;
        }
    }

	// top bottom
    //for(i=s_x; i<=e_x; i++){
    for (i=0; i<nx+2; i++){
		j = 0;
			x = ((double)i) / ((double)(nx+2));
			y = ((double)j) / ((double)(ny+2));
		
			a[i][0] = analytic_sol(x, y); //bottom = 0
			b[i][0] = analytic_sol(x, y);

		j = nx+1;
			x = ((double)i) / ((double)(nx+2));
			y = ((double)j) / ((double)(ny+2));

			a[i][nx+1] = analytic_sol(x, y);;
			b[i][nx+1] = analytic_sol(x, y);;
		}

	// left right
    //for(j=s_y; j<=e_y; j++){
    for (j=0; j<nx+2; j++){
		i = 0;
			x = ((double)i) / ((double)(nx+2));
			y = ((double)j) / ((double)(ny+2));
	
			a[0][j] = analytic_sol(x, y);
			b[0][j] = analytic_sol(x, y);

		i = nx+1;
			x = ((double)i) / ((double)(nx+2));
			y = ((double)j) / ((double)(ny+2));

			a[nx+1][j] = analytic_sol(x, y);
			b[nx+1][j] = analytic_sol(x, y);
    }
}

int MPE_Decomp1d(int n, int size, int rank, int *s, int *e){
    int nlocal, deficit;

    nlocal = n / size;
    *s = rank * nlocal + 1;

    deficit = n % size;
    *s = *s + ((rank < deficit) ? rank : deficit);

    if (rank < deficit) nlocal++;
    	*e = *s + nlocal - 1;
    if (*e > n || rank == size-1) *e = n;
    	return MPI_SUCCESS;
}

void GatherGrid2d(double grid[][maxn], double proc_grid[][maxn], int nx, int ny, int* dims, int* coords, MPI_Comm comm){
    int myid, size;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &size);
    MPI_Barrier(comm);
    MPI_Status status;

    if (myid != 0){
        MPI_Send(proc_grid, (nx+2)*(ny+2), MPI_DOUBLE, 0, 100 + myid, comm);
		//fprintf(stderr ,"SENDING PART OF GRID FROM RANK=%d\n", myid);
    }

    if (myid == 0){
        int s_x, s_y, e_x, e_y;

        // combinting on proc 0
        MPE_Decomp1d(nx+2, dims[0], coords[0], &s_x, &e_x);
        MPE_Decomp1d(ny+2, dims[1], coords[1], &s_y, &e_y);
        for (int i=s_x-1; i<=e_x; i++){
            for (int j=s_y-1; j<=e_y; j++){
                grid[i][j] = proc_grid[i][j];
            }
        }

        // collecting from other procs
        double temp_grid[maxn][maxn];
        for (int k=1; k<size; k++){
            init_full_grid(temp_grid);
            MPI_Recv(temp_grid, (nx+2)*(ny+2), MPI_DOUBLE, k, 100 + k, comm, &status);

            int k_coords[2];
            MPI_Cart_coords(comm, k, 2, k_coords);

            MPE_Decomp1d(nx+2, dims[0], k_coords[0], &s_x, &e_x);
            MPE_Decomp1d(ny+2, dims[1], k_coords[1], &s_y, &e_y);

            for (int i=s_x-1; i<=e_x; i++){
                for (int j=s_y-1; j<=e_y; j++){
                    grid[i][j] = temp_grid[i][j];
                }
            }
        }
    }
}

double analytic_sol(double x, double y){
    return y / ( (1+x)*(1+x) + y*y);
}


void analytic_sol_matrix_2d(double u[][maxn], int nx, int ny){
	double x, y;
    for (int i=0; i<nx+2; i++){
        for (int j=0; j<ny+2; j++){
			x = ((double)i) / ((double)(nx+2));
			y = ((double)j) / ((double)(ny+2));
	
			u[i][j] = analytic_sol(x, y);
        }
    }
}

double max_diff(double a[][maxn], double b[][maxn], int nx, int ny, int* coords, int* dims){
    double temp_max_diff = 0;

    for (int i=0; i<nx+2; i++){
        for (int j=0; j<ny+2; j++){
			double diff = (a[i][j] - b[i][j]) * (a[i][j] - b[i][j]);
			if (diff > temp_max_diff){
				temp_max_diff = diff;
			}
        }
    }
    return temp_max_diff;
}

double griddiff_2d(double a[][maxn], double b[][maxn], int s_x, int e_x, int s_y, int e_y){
	double sum = 0;
	double tmp;

	for(int i=s_x; i<=e_x; i++){
		for(int j=s_y; j<=e_y; j++){
			tmp = a[i][j] - b[i][j];
			sum = sum + tmp*tmp;
		}
	}

	return sum;
}

int check_covergence(double a[][maxn], double b[][maxn], int s_x, int e_x, int s_y, int e_y, double* ldiff, double* glob_diff, int tol, int it){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	*ldiff = griddiff_2d(a, b, s_x, e_x, s_y, e_y);
	MPI_Allreduce(ldiff, glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if (rank==0 && it%10==0){
		printf("(rank %d) locdiff: %lf; glob_diff: %lf\n", rank, *ldiff, *glob_diff);
	}
	if (*glob_diff < tol){
		if (rank==0){ 
			printf("iterative solve converged\n");
		}
		return 1;
	}
	return 0;
}

/* prints to stdout in GRID view */
void print_full_grid(double x[][maxn]){
    for(int i=0; i<maxn; i++){
    	for(int j=0; j<maxn; j++){
            if(x[i][j] < 10000.0){
				printf("|%2.6lf| ", x[i][j]);
			} else {
				printf("%9.2lf "  , x[i][j]);
			}
				}
        printf("\n");
    }
}


void print_in_order(double x[][maxn], MPI_Comm comm){
    int myid, size;

    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &size);
    MPI_Barrier(comm);
    sleep(1);
    MPI_Barrier(comm);

    for (int i=0; i<size; i++){
        if (i == myid ){
            printf("proc %d\n",myid);
            print_full_grid(x);
        }
        fflush(stdout);
        usleep(500);
        MPI_Barrier(comm);
    }
}

// writes grid to file
void write_grid(double grid[][maxn], char* file_name, int nx){
    FILE *fp = fopen(file_name, "w");

    if (fp == NULL) {
        printf("Error opening file!\n");
    }

	for (int j=nx+2-1; j>=0; j--){
		for (int i=0; i<nx+2; i++){
            fprintf(fp, "%lf ", grid[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

