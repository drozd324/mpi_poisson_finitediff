#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "poisson.h"
#include "jacobiseq.h"

#define maxit 1000

void init_full_grids(double a[][maxn], double b[][maxn] ,double f[][maxn]);
void print_full_grid(double x[][maxn]);
void init_basic_bv(double a[][maxn], double b[][maxn], double f[][maxn], int n);


int main(int argc, char **argv)
{
  double a[maxn][maxn], b[maxn][maxn], f[maxn][maxn];
  int n;

  int it;
  double ldiff;
  /* double t1, t2; */
  double tol=1.0E-11;

  /* set the size of the problem */
  if(argc > 2){
    fprintf(stderr,"---->Usage: %s <n>\n",argv[0]);
    exit(1);
  }
  if( argc == 2 ){
    n = atoi(argv[1]);
  }

  if( n > maxn-2 ){
    fprintf(stderr,"grid size too large\n");
    exit(1);
  }

  printf("n = %d\n",n);

  init_full_grids(a, b, f);
  /* print_full_grid(a); */

  init_basic_bv(a, b, f, n);
  print_full_grid(a);

  ldiff = 1000;
  for(it=0; it<maxit; it++){

    sweep(a, f, n, b);
    sweep(b, f, n, a);

    ldiff = griddiffseq(a, b, n);

    if(it%10==0){
      printf("diff beteween grids: %lf; \n",ldiff);
    }
    if( ldiff < tol ){
  	printf("iterative solve converged\n");
      break;
    }

  }
    
  printf("DONE! (it: %d)\n",it);

  if( it == maxit ){
    fprintf(stderr,"Failed to converge\n");
  }
  /* printf("Run took %lf s\n",t2-t1); */

    /* print_grid_to_file("grid", a,  n, ny); */
    print_full_grid(a);

  return 0;
}

void init_basic_bv(double a[][maxn], double b[][maxn], double f[][maxn], int n)
{
  int i,j;
  double left, bottom, right, top;

  left   = -1.0;
  bottom = 1.0;
  right  = 2.0;
  top    = 3.0;  

  /* set everything to 0 first; not really necessary */
  for(i=0; i<=n+1; i++){
    for(j=0; j <= n+1; j++){
      a[i][j] = 0.0;
      b[i][j] = 0.0;
      f[i][j] = 0.0;
    }
  }

  /* deal with boundaries */
  for(i=1; i<=n; i++){
    a[i][0] = bottom;
    b[i][0] = bottom;
    a[i][n+1] = top;
    b[i][n+1] = top;
  }

  for(j=1; j<n+1; j++){
      a[0][j] = left;
      b[0][j] = left;
      a[n+1][j] = right;
      b[n+1][j] = right;
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
