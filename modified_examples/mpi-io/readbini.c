#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[])
{
  FILE *fp;
  char *fname;
  int tmp, cnt;
  
  if(argc != 2){
    fprintf(stderr,"Usage: %s <filename>\n",argv[0]);
    exit(10);
  }

  fname = argv[1];

  fp = fopen(fname, "rb");
  
  cnt=0;
  while(!feof(fp)){
    fread(&tmp, sizeof(int), 1, fp);
    printf("int %d: %d\n",cnt,tmp);
    cnt++;
  }
  printf("DONE\n");
  return 0;
}
