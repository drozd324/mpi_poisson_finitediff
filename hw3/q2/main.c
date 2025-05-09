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
	int root = 2;
	char mat_fname[100] = "mat-d20-b5-p4.bin";
	char vec_fname[100] = "x-d20.txt.bin";

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == root)
		printf("ROOT=%d\n\n", root);
	if (size != 4){
		printf("THIS CASE SHOULD BE RAN WITH 4 PROCESSES\n");
		MPI_Abort(MPI_COMM_WORLD, 3);	
	}


	//==================== reading matrix ========================
	MPI_File fh_mat;
	int N;
	int n = 5;

	MPI_File_open(MPI_COMM_WORLD, mat_fname, MPI_MODE_RDWR, MPI_INFO_NULL, &fh_mat);
	MPI_File_read(fh_mat, &N, 1, MPI_INT, MPI_STATUS_IGNORE);
	
	double* blocks = malloc(size*n*n * sizeof(double));

	MPI_Datatype filetype_mat;
	MPI_Type_vector(size*n*n, size*n*n, 0, MPI_DOUBLE, &filetype_mat);	
	MPI_Type_commit(&filetype_mat);

	printf("rank=%d, offset=%d\n", rank, rank*(size*n*n)*sizeof(double));
	
	MPI_File_set_view(fh_mat, sizeof(int) + rank*(size*n*n)*sizeof(double), MPI_DOUBLE, filetype_mat, "native", MPI_INFO_NULL);
	MPI_File_read_all(fh_mat, blocks, size*n*n, MPI_DOUBLE, MPI_STATUS_IGNORE);
	
	if (rank == root){
		for (int k=0; k<size; k++){	
			for (int i=0; i<n; i++){	
				for (int j=0; j<n; j++){	
					printf("%lf ", blocks[k*n*n + i*n + j]);
				}
				printf("\n");
			}
			printf("\n");
		}	
	}
	

	MPI_File_close(&fh_mat);
	MPI_Type_free(&filetype_mat);

	//==================== reading vector ========================
	MPI_File fh_vec;
	int vec_len;

	MPI_File_open(MPI_COMM_WORLD, vec_fname, MPI_MODE_RDWR, MPI_INFO_NULL, &fh_vec);
	MPI_File_read(fh_vec, &vec_len, 1, MPI_INT, MPI_STATUS_IGNORE);
	
	int sub_vec_len = vec_len/size;
	double* sub_vec = malloc(sub_vec_len * sizeof(double));

	MPI_Datatype filetype_vec;
	MPI_Type_vector(sub_vec_len, sub_vec_len, 0, MPI_DOUBLE, &filetype_vec);	
	MPI_Type_commit(&filetype_vec);

	MPI_File_set_view(fh_vec, rank*sub_vec_len*sizeof(double), MPI_DOUBLE, filetype_vec, "native", MPI_INFO_NULL);
	MPI_File_read_all(fh_vec, sub_vec, sub_vec_len, MPI_DOUBLE, MPI_STATUS_IGNORE);
	
	if (rank == root){
		for (int i=0; i<n; i++){	
			printf("%lf ", sub_vec[i]);
		}
		printf("\n");
		printf("\n");
	}
	

	MPI_File_close(&fh_vec);
	MPI_Type_free(&filetype_vec);

	//==================== distributed matrix vector product ========================
	
	double* b = calloc(N, sizeof(double));
	double* b_end = calloc(N, sizeof(double));

	mat_vect_prod(blocks, sub_vec, b, N, sub_vec_len);
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
