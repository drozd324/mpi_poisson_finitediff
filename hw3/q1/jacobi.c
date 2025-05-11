#include "poisson.h"
#include "jacobi.h"

void sweep2d(double a[][maxn], double f[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y, double b[][maxn]){

	double h;
	h = 1.0/((double)(nx+2));

	for(int i=s_x; i<=e_x; i++){
		for(int j=s_y; j<=e_y; j++){
			b[i][j] = 0.25 * (a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1]  - h*h*f[i][j]);
		}
	}
}


void PutGet_winfence_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){
	
	int rank;
	MPI_Comm_rank(comm, &rank);
	
	MPI_Datatype vect;
	int col_len = e_x+1 - s_x;
	int row_len = e_y+1 - s_y;

	MPI_Type_vector(col_len, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);

	MPI_Win win;
	MPI_Win_create(x, (nx+2)*(nx+2) * sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);
	
	MPI_Win_fence(0, win);

	MPI_Put(&x[s_x][s_y], row_len, MPI_DOUBLE, nbrdown , (s_x  )*(nx+2) + s_y  , row_len, MPI_DOUBLE, win);
	MPI_Put(&x[s_x][e_y], 1      , vect      , nbrright, (s_x  )*(nx+2) + e_y  , 1      , vect      , win);
	MPI_Put(&x[s_x][s_y], 1      , vect      , nbrleft , (s_x  )*(nx+2) + s_y  , 1      , vect      , win);
	MPI_Put(&x[e_x][s_y], row_len, MPI_DOUBLE, nbrup   , (e_x  )*(nx+2) + s_y  , row_len, MPI_DOUBLE, win);

	MPI_Win_fence(0, win);

	MPI_Type_free(&vect);
	MPI_Win_free(&win);
}

void PutGet_activeSync_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){
	
	MPI_Datatype vect;

	int col_len = e_x+1 - s_x;
	int row_len = e_y+1 - s_y;

	MPI_Type_vector(col_len, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);
		
	int nbrs[4] = {nbrleft, nbrright, nbrup, nbrdown};
	int group_size = 4;
	int valid_nbrs[4];
	int num_valid_nbrs = 0;
	for (int i=0; i<4; i++){
		if (nbrs[i] != MPI_PROC_NULL){
			valid_nbrs[num_valid_nbrs] = nbrs[i];	
			num_valid_nbrs++;
		}
	}

	MPI_Group nbr_group;
	MPI_Group world_group;
	MPI_Comm_group(comm, &world_group);
	MPI_Group_incl(world_group, num_valid_nbrs, valid_nbrs, &nbr_group);

	MPI_Win win;
	MPI_Win_create(x, (nx+2)*(nx+2) * sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);

	MPI_Win_post(nbr_group, 0, win);
	MPI_Win_start(nbr_group, 0, win);

	MPI_Put(&x[s_x][s_y], row_len, MPI_DOUBLE, nbrdown , (s_x  )*(nx+2) + s_y  , row_len, MPI_DOUBLE, win);
	MPI_Put(&x[s_x][e_y], 1      , vect      , nbrright, (s_x  )*(nx+2) + e_y  , 1      , vect      , win);
	MPI_Put(&x[s_x][s_y], 1      , vect      , nbrleft , (s_x  )*(nx+2) + s_y  , 1      , vect      , win);
	MPI_Put(&x[e_x][s_y], row_len, MPI_DOUBLE, nbrup   , (e_x  )*(nx+2) + s_y  , row_len, MPI_DOUBLE, win);
	
	MPI_Win_complete(win);
	MPI_Win_wait(win);

	MPI_Group_free(&nbr_group);
	MPI_Group_free(&world_group);
	MPI_Type_free(&vect);
	MPI_Win_free(&win);
}

void Isend_Irecv_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){

	int col_len = e_x+1 - s_x;
	int row_len = e_y+1 - s_y;

	MPI_Request reqs[8];
	MPI_Datatype vect;
	MPI_Type_vector(row_len, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);

	MPI_Isend(&x[s_x][e_y], 1      , vect      , nbrup   , 0, comm, &reqs[7]);
	MPI_Isend(&x[e_x][s_y], col_len, MPI_DOUBLE, nbrright, 1, comm, &reqs[6]);
	MPI_Isend(&x[s_x][s_y], col_len, MPI_DOUBLE, nbrleft , 2, comm, &reqs[5]);
	MPI_Isend(&x[s_x][s_y], 1      , vect      , nbrdown , 3, comm, &reqs[4]);

	MPI_Irecv(&x[s_x][s_y-1], 1      , vect      , nbrdown , 0, comm, &reqs[0]);
	MPI_Irecv(&x[s_x-1][s_y], col_len, MPI_DOUBLE, nbrleft, 1, comm, &reqs[1]);
	MPI_Irecv(&x[e_x+1][s_y], col_len, MPI_DOUBLE, nbrright , 2, comm, &reqs[2]);
	MPI_Irecv(&x[s_x][e_y+1], 1      , vect      , nbrup   , 3, comm, &reqs[3]);

	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	MPI_Type_free(&vect);
}
