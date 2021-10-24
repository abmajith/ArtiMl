#include <stdio.h>
#include <stdlib.h>
#define Maxlen 1000



/* test this c program by 
 * gcc -Wall -o a reading_space_speratredfile.c 
 * ./a testwriting.txt   */
int getarray(char *str, double *arr, int k){
  double d;
  int i = 0;
  char *ptr = NULL;
  ptr = (char *) malloc(sizeof(char) * Maxlen);
  char str1[Maxlen];
  int j = 0;
  while( *str != '\0'){
    if ((str1[i] = *str) != ' '){
      i++;
    } else{
      d = strtod(str1, &ptr);
      if (j < k){
        *(arr + j) = d;
        j++;
      }
      i = 0;
    }
    str++;
  }
  return 0;
}

int main(int argc, char *argv[]){ /* give the space separated file name here, it will access data file, and extract data vector one by one */
  FILE *fp;
  int c;
  fp = fopen(argv[1], "r");
  if (fp == NULL){
    fprintf(stderr, "%s can't open %s\n", argv[0],argv[1] );
    exit(1);
  }
  int i;
  char str[Maxlen];
  int k = 3;
  double *array = NULL;
  array = (double*) malloc(sizeof(double) * k);
  while ( (c = getc(fp)) != EOF){
    if (c != '\n'){
      str[i] = c;
      i++;
    }else {
      str[i] = ' ';
      str[i+1] = '\0';
      i = 0;
      getarray(str, array, k);/* here we will get the k double values in this array */
      for (int l = 0; l < k; l++)
        printf("%lf\t",*(array + l));
      printf("\n");
    }
  }
  if (ferror(stdout)){
    fprintf(stderr, "%s: error wirting stdout\n", argv[1]);
  }
  return 0;
}
