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
    int i, j;
    double u, w;

    /* set everything to 0 first */
    for(i=s_x; i<=e_x; i++){
        for(j=s_y; j<=e_y; j++){
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            f[i][j] = 0.0;
        }
    }

  /* deal with boundaries */
  // top bottom
    //for(i=s_x; i<=e_x; i++){
    for(i=0; i<maxn; i++){
        a[i][0] = 0; //bottom = 0
        b[i][0] = 0;

        w = 1 / ( (1 + ((double)i/(double)maxn))*(1 + ((double)i/(double)maxn)) + 1); // 1/((1+x)*(1+x) + 1)
        a[i][maxn-1] = w;
        b[i][maxn-1] = w;
    }

	// left right
    //for(j=s_y; j<=e_y; j++){
    for(j=0; j<maxn; j++){
        u = ((double)j/(double)maxn)/(1 + ((double)j/(double)maxn)*((double)j/(double)maxn));  // y/(1 + y*y)
        a[0][j] = u;
        b[0][j] = u;

        w = ((double)j/(double)maxn)/(4 + ((double)j/(double)maxn)*((double)j/(double)maxn));  // y/(4 + y*y)
        a[maxn-1][j] = w;
        b[maxn-1][j] = w;
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
        MPI_Send(proc_grid, (nx+2)*(nx+2), MPI_DOUBLE, 0, 100 + myid, comm);
    }

    if (myid == 0){
        int s1_x, s1_y, e1_x, e1_y;

        // combinting on proc 0
        MPE_Decomp1d(nx+2, dims[1], coords[1], &s1_x, &e1_x);
        MPE_Decomp1d(ny+2, dims[0], coords[0], &s1_y, &e1_y);
        for (int i=s1_x; i<=e1_x; i++){
            for (int j=s1_y; j<=e1_y; j++){
                grid[i][j] = proc_grid[i][j];
            }
        }

        // collecting from other procs
        double temp_grid[maxn][maxn];
        for (int k=1; k<size; k++){
            init_full_grid(temp_grid);
            MPI_Recv(temp_grid, (nx+2)*(nx+2), MPI_DOUBLE, k, 100 + k, comm, &status);

            int k_coords[2];
            MPI_Cart_coords(comm, k, 2, k_coords);

            MPE_Decomp1d(nx+2, dims[1], k_coords[1], &s1_x, &e1_x);
            MPE_Decomp1d(ny+2, dims[0], k_coords[0], &s1_y, &e1_y);

            for (int i=s1_x; i<=e1_x; i++){
                for (int j=s1_y; j<=e1_y; j++){
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
    MPE_Decomp1d(nx+2, dims[1], coords[1], &s1_x, &e1_x);
    int s1_y, e1_y;
    MPE_Decomp1d(ny+2, dims[0], coords[0], &s1_y, &e1_y);

    for (int i=s1_x; i<e1_x; i++){
        for (int j=s1_y; j<e1_y; j++){
            u[i][j] = analytic_sol((double)i/(double)nx, (double)j/(double)ny);
        }
    }
}

double compute_mse_2d(double a[][maxn], double b[][maxn], int nx, int ny, int* coords, int* dims){
    double sum = 0;
    int s1_x, e1_x;
    MPE_Decomp1d(nx+2, dims[1], coords[1], &s1_x, &e1_x);
    int s1_y, e1_y;
    MPE_Decomp1d(ny+2, dims[0], coords[0], &s1_y, &e1_y);

    for (int i=s1_x; i<e1_x; i++){
        for (int j=s1_y; j<e1_y; j++){
            sum += (a[i][j] - b[i][j]) * (a[i][j] - b[i][j]);
        }
    }
    return sum / (((double)e1_y - (double)s1_y + 1) * ((double)e1_x - (double)s1_x + 1));
}


/* prints to stdout in GRID view */
void print_full_grid(double x[][maxn]){
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


void print_in_order(double x[][maxn], MPI_Comm comm){
    int myid, size;
    int i;

    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &size);
    MPI_Barrier(comm);
    sleep(1);
    MPI_Barrier(comm);

    for(i=size-1; i>=0; i--){
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

