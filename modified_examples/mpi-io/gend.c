#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[])
{
  FILE *fp;
  char fname[100];
  int nprocs, items_per_proc;
  int totalitems;
  double *data;
  int icnt;
  int i, j;
  int ret=0;

  if( argc != 3 ){
    fprintf(stderr, "Usage: %s <number of processors> <number of doubles per processor>\n", argv[0]);
    exit(1);
  }
  nprocs = atoi(argv[1]);
  items_per_proc = atoi(argv[2]);
  totalitems = nprocs * items_per_proc;
  printf("nprocs: %d; items_per_proc: %d; So: total:  %d\n",
	 nprocs,items_per_proc,totalitems);
  sprintf(fname, "doublempiotfile_%d_%d",nprocs,items_per_proc);
  printf("Output file name: %s\n",fname);
  fp = fopen(fname,"w");
  if(fp==NULL){
    fprintf(stderr,"Error: failure to open file\n");
    exit(4);
  }
  
  data = (double *)malloc(totalitems * sizeof(double));
  
  icnt = 0;
  for(i=0; i < nprocs; i++){
    for(j=0; j<items_per_proc; j++){
      data[icnt] = (double)i;
      icnt++;
    }
  }
  if(icnt != totalitems){
    fprintf(stderr,"counter(%d) and totalitems(%d) mismatch\n",icnt, totalitems);
    exit(2);
  }

  ret = fwrite(data, sizeof(double), totalitems, fp);

  if(ret != totalitems){
    fprintf(stderr,"return items not equal to total items\n");
    exit(3);
  }

  ret=fclose(fp);
  if(ret != 0){
    fprintf(stderr,"Failure to close file\n");
  }

  return 0;
}
