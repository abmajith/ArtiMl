#include <stdio.h>
#include <stdlib.h>
#include "knn.h"

int main(){
  int vec_dim = 2;
  int No_data_pts = 3;
  double *std;
  double *mean;
  std = (double*) malloc(sizeof(double) * vec_dim);
  mean = (double*) malloc(sizeof(double) * vec_dim);
  double *data[No_data_pts];
  for (int i = 0; i < No_data_pts; i++){
    data[i] = (double*) malloc(sizeof(double) * vec_dim);
  }
  for (int i = 0; i < No_data_pts; i++){
    for (int j = 0; j < vec_dim; j++){
      *(data[i] + j) = (i + 1) * (j + 2);
    }
  }
  Normalize_KNN_vec_data_Db_prec(vec_dim, No_data_pts, data, mean, std);
  for (int i = 0; i < No_data_pts; i++){
    for (int j = 0; j < vec_dim; j++){
      printf("%f \t", *(data[i] + j) );
    }
    printf("\n");
  }
  for (int i = 0; i < vec_dim; i++){
    printf("%f \t %f \n",mean[i], std[i]);
  }
  int data_label[3];
  data_label[0] = 1;
  data_label[1] = 21;
  data_label[2] = 2;
  int querry_data[2];
  querry_data[0] = 4.1;
  querry_data[1] = 6.1;
  int *querry_result_label = NULL;
  querry_result_label = (int*) malloc(sizeof(int));
  *querry_result_label = 0;
  querry_KNN_maj_voting_exhaustive(vec_dim, No_data_pts, data, data_label, querry_data, 1, querry_result_label);
  printf("The argmax is %d\n", *querry_result_label);
  return 0;
}
