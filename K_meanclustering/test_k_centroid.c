#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "kmean_clustering.h"

int main(void){
  int vec_dim = 2;
  int No_Data_pts = 3;
  double *data[3];
  double rangeMin[vec_dim];
  double rangeMax[vec_dim];
  double absmax[vec_dim];
  for (int i = 0; i < No_Data_pts; i++){
    data[i] = (double *) malloc(sizeof(double) * vec_dim);
    for (int j = 0; j < vec_dim; j++ ){
      *(data[i] + j) = 1 * (i + 1) * (j + 1);
      printf("%f\t",*(data[i] + j));
    }
    printf("\n");
  }

  KMeanCl_Normalize_data_Db(No_Data_pts, vec_dim, data, rangeMin, rangeMax, absmax);
  for (int i = 0; i < No_Data_pts; i++){
    for (int j = 0; j < vec_dim; j++){
      printf("%f\t", *(data[i] + j));
    }
    printf("\n");
  }

  double *Centroid[2];
  for (int l = 0; l < 2; l++){
    Centroid[l] = (double *) malloc(sizeof(double) * vec_dim);
  }
  double *NewCentroid[2];
  for (int l = 0; l < 2; l++){
    NewCentroid[l] = (double *) malloc(sizeof(double) * vec_dim);
  }
  KMeanCl_Initialize_Centroids_Db(2,vec_dim, rangeMin, rangeMax, absmax,Centroid);
  for (int l = 0; l < 2; l++){
    for (int j = 0; j < vec_dim; j++){
      printf("%f\t", *(Centroid[l] + j) );
    }
    printf("\n");
  }
  int argmin[3];
  KMeanCl_minarg_to_centroid_Db(2,vec_dim,No_Data_pts, data,Centroid,argmin);
  for (int i = 0; i < No_Data_pts; i++){
    printf("%d\n", argmin[i]);
  }
  for (int z = 0; z < 5; z++){
    KMeanCl_update_centroids_Db(2, vec_dim, No_Data_pts, data, Centroid, NewCentroid, argmin);
    for (int l = 0; l < 2; l++){
      for (int j = 0; j < vec_dim; j++){
        printf("%f\t", *(NewCentroid[l] + j) );
      }
      printf("\n");
    }
  }
  return 0;
}
