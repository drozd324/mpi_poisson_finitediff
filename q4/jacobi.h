void sweep1d(double a[][maxn], double f[][maxn], int nx,
	     int s, int e, double b[][maxn]);

void sweep2d(double a[][maxn], double f[][maxn], int nx,
		int s_x, int e_x, int s_y, int e_y, double b[][maxn]);

void exchang1(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
	      int nbrbottom, int nbrtop);

void exchang2(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
	      int nbrleft, int nbrright);

void exchang3(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
	      int nbrleft, int nbrright);

void exchang3_2d(double x[][maxn], int nx,
                int s_x, int e_x, int s_y, int e_y ,
                int nbrleft, int nbrright, int nbrup, int nbrdown,
		MPI_Comm comm);


void exchangi1(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
	       int nbrleft, int nbrright);

void exchangi2(double x[][maxn], int nx, 
		int s_x, int e_x, int s_y, int e_y ,
		int nbrleft, int nbrright, int nbrup, int nbrdown, 
		MPI_Comm comm);

void nbxchange_and_sweep(double u[][maxn], double f[][maxn], int nx, int ny,
			 int s, int e, double unew[][maxn], MPI_Comm comm,
			 int nbrleft, int nbrright);

double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e);

double griddiff_2d(double a[][maxn], double b[][maxn], int s_x, int e_x, int s_y, int e_y);
