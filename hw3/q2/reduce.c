#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

void mat_vect_prod(double* A, double* x, double* b, int n, int m){
	memset(b, 0, n*sizeof(double));

    for (int i=0; i<n; i++){
		for (int j=0; j<m; j++){
			b[i] += A[i*m + j] * x[j];
		}
	}
}

int main(int argc, char* argv[]){
	int rank, size;
	int root = 1;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double* blocks = malloc(n*N * sizeof(double));
	double* sub_vec = malloc(n * sizeof(double));

	double* b = calloc(N, sizeof(double));
	double* b_end = calloc(N, sizeof(double));

	mat_vect_prod(blocks, sub_vec, b, N, n);
	MPI_Reduce(b, b_end, N, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

	if (rank == root){
		for (int i=0; i<N; i++){
			printf("%lf ", b_end[i]); 
		}
		printf("\n"); 
		printf("\n"); 
	}

	
	free(b);
	free(b_end);
	free(sub_vec);
	free(blocks);
	MPI_Finalize();
	return 0;
}
