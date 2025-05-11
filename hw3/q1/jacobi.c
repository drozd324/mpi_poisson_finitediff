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
	
	int rank;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm newcomm;
	MPI_Datatype vect;

	int col_len = e_x+1 - s_x;
	int row_len = e_y+1 - s_y;

	MPI_Type_vector(col_len, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);
		
	MPI_Group local_group;
	MPI_Comm_group(comm, &local_group);

	//int group_ranks[5] = {rank, nbrup, nbrdown, nbrleft, nbrright};
	//MPI_Group_incl(local_group, 5, group_ranks, local_group);

//	MPI_Comm_split(comm, colour, 0, &newcomm);
//	MPI_Comm_rank(newcomm, &nrank);
//  	MPI_Comm_size(newcomm, &nsize);

	MPI_Win win;
	MPI_Win_create(x, (nx+2)*(nx+2) * sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);

	MPI_Win_post(local_group, 0, win);
	MPI_Win_start(local_group, 0, win);
	
	MPI_Put(&x[s_x][s_y], row_len, MPI_DOUBLE, nbrdown , (s_x  )*(nx+2) + s_y  , row_len, MPI_DOUBLE, win);
	MPI_Put(&x[s_x][e_y], 1      , vect      , nbrright, (s_x  )*(nx+2) + e_y  , 1      , vect      , win);
	MPI_Put(&x[s_x][s_y], 1      , vect      , nbrleft , (s_x  )*(nx+2) + s_y  , 1      , vect      , win);
	MPI_Put(&x[e_x][s_y], row_len, MPI_DOUBLE, nbrup   , (e_x  )*(nx+2) + s_y  , row_len, MPI_DOUBLE, win);
	
	MPI_Win_complete(win);
	MPI_Win_wait(win);

	MPI_Type_free(&vect);
	MPI_Win_free(&win);
}

void Isend_Irecv_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){

	MPI_Request reqs[8];
	MPI_Datatype vect;
	
	int col_len = e_x+1 - s_x;
	int row_len = e_y+1 - s_y;
	MPI_Type_vector(col_len, 1, nx+2, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);

	MPI_Irecv(&x[s_x][s_y-1], row_len, MPI_DOUBLE, nbrdown , 0, comm, &reqs[0]);
	MPI_Irecv(&x[s_x-1][s_y], 1      , vect      , nbrleft , 1, comm, &reqs[1]);
	MPI_Irecv(&x[e_x+1][s_y], 1      , vect      , nbrright, 2, comm, &reqs[2]);
	MPI_Irecv(&x[s_x][e_y+1], row_len, MPI_DOUBLE, nbrup   , 3, comm, &reqs[3]);

	MPI_Isend(&x[s_x][e_y], row_len, MPI_DOUBLE, nbrup   , 0, comm, &reqs[7]);
	MPI_Isend(&x[e_x][s_y], 1      , vect      , nbrright, 1, comm, &reqs[6]);
	MPI_Isend(&x[s_x][s_y], 1      , vect      , nbrleft , 2, comm, &reqs[5]);
	MPI_Isend(&x[s_x][s_y], row_len, MPI_DOUBLE, nbrdown , 3, comm, &reqs[4]);


	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	MPI_Type_free(&vect);
}



void Sendrecv_2d(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm){
	
	int rank;
	MPI_Comm_rank(comm, &rank);	

	int row_len = e_x + 1 - s_x;
	int col_len = e_y + 1 - s_y;

	MPI_Datatype vect;
	MPI_Type_vector(row_len, 1, maxn, MPI_DOUBLE, &vect);
	MPI_Type_commit(&vect);
	
	
	/*
	int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 int dest, int sendtag,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int source, int recvtag,
				 MPI_Comm comm, MPI_Status * status)
	*/

	// to left
	MPI_Sendrecv(&x[s_x  ][s_y], 1, vect,
				 nbrleft, 0, 
		         &x[s_x-1][s_y], 1, vect,
				 rank, 0, 
		         comm, MPI_STATUS_IGNORE);
	// to right
	MPI_Sendrecv(&x[e_x  ][s_y], 1, vect,
				 nbrright , 1,
		         &x[s_x+1][s_y], 1, vect,
				 rank, 1,
     		     comm, MPI_STATUS_IGNORE);

	MPI_Type_free(&vect);
	
	// to top
	MPI_Sendrecv(&x[s_x][s_y  ], row_len, MPI_DOUBLE,
				 nbrup  , 2, 
		         &x[s_x][s_y-1], row_len, MPI_DOUBLE,
				 rank, 2,
		         comm, MPI_STATUS_IGNORE);
	// send to bottom
	MPI_Sendrecv(&x[e_x][s_y  ], row_len, MPI_DOUBLE,
				 nbrdown, 3,
		         &x[s_x][s_y+1], row_len, MPI_DOUBLE,
				 rank, 3,
     		     comm, MPI_STATUS_IGNORE);
}

