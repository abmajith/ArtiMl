#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>

int KMeanCl_Normalize_data_Db(int No_Data_pts, int vec_dim, double *data[], double Min_range[], double Max_range[], double Abs_max[]){
int d = vec_dim;
  int N = No_Data_pts;
  if (d <= 0){
    printf("The dimension of input datas is a NULL vector, exiting from the program\n");
    exit(0);
  }
  if (N <= 0){
    printf("There are no data points, exiting from the program\n");
    exit(0);
  }
  double max[d];
  double min[d];
  double diff[d];
  double absmax[d];
  /* transfering the data from the dimenstion R^d to [-1,1]^d
   * */
  for (int j = 0; j < d; j++){
    max[j] = *(data[0] + j);
    min[j] = *(data[0] + j);
    diff[j] = 0.0;
    for (int i = 0; i < N; i++){
      if (max[j] < *(data[i] + j))
        max[j] = *(data[i] + j);
      if (min[j] > *(data[i] + j))
        min[j] = *(data[i] + j);
    }
    diff[j] = max[j] - min[j];
    if (min[j] > 0.0)
      absmax[j] = max[j];
    else
      if (max[j] < -1 * min[j]){
        absmax[j] = -1 * min[j];
      } else
        absmax[j] = max[j];
    if (diff[j] == 0.0){
      printf("The %d (st or th) component of all the data points are %f, one can reduce the dimension by eliminating the %d compoents in the datas\n",(j + 1), (diff[j]), (j + 1));
    }
  }
  for (int j = 0; j < d; j++){
    for(int i = 0; i < N; i++)
      *(data[i] + j) /= absmax[j];
  }
  for (int j = 0; j < d; j++){
    Min_range[j] = max[j];
    Max_range[j] = min[j];
    Abs_max[j] = absmax[j];
  }
  printf("Input data are normalized\n");
  return 0;
}


int KMeanCl_Initialize_Centroids_Db(int k, int vec_dim, double Min_range[], double Max_range[], double Abs_max[], double *Centroid[]){
  if (k <= 0){
    printf("There is nothing to intialize for %d cluster cetroids\n program exiting\n", k);
    exit(0);
  }
  if (vec_dim <= 0){
    printf("The centroid vector dimension is NULL, nothing to be done, exiting the program!\n");
    exit(0);
  }
  double rand_num;
  for (int l = 0; l < k; l++){
    for (int j = 0; j < vec_dim; j++){
      if (Min_range[j] == Max_range[j]){
        *(Centroid[l] + j) = Min_range[j];
      } else {
        srand((j + 1) * (l + 1) * time(NULL));
        rand_num = (double) rand() / (double) RAND_MAX;
        *(Centroid[l] + j) = Min_range[j] + rand_num * (Max_range[j] - Min_range[j]);
        *(Centroid[l] + j) /= Abs_max[j];
      }
    }
  }
  printf("%d centroid vectors of dimension %d are intialized using uniform distribution\n",k,vec_dim);
  return 0;
}

int KMeanCl_minarg_to_centroid_Db(int k, int vec_dim, int No_data_pts, double *data[], double *Centroid[], int argmin[]){
  if (vec_dim <= 0){
    printf("The data points are null dimension, exiting the program\n");
    exit(0);
  }
  if (No_data_pts <= 0){
    printf("There are no data points, exiting the program\n");
    exit(0);
  }
  int N= No_data_pts;
  int d = vec_dim;
  int minarg; // to store the current minimum argument
  double lsqdist; // to store the current mininum square distance
  double sqdist; // calcualte the square distance
  for (int i = 0; i < N; i++){
    lsqdist = DBL_MAX;
    minarg = 1;
    for (int l = 0; l < k; l++){
      sqdist = 0.0;
      for (int j = 0; j < d; j++){
        sqdist += ( *(Centroid[l] + j) - *(data[i] + j)) * ( *(Centroid[l] + j) - *(data[i] + j) );
      }
      if (lsqdist > sqdist){
        minarg = l + 1;
        lsqdist = sqdist;
      }
    }
    argmin[i] = minarg;
  }
  return 0;
}

int KMeanCl_update_centroids_Db(int k, int vec_dim, int No_data_pts, double *data[], double *Centroid[], double *NewCentroid[], int argmin[] ){
  if (k <= 0){
    printf("There is nothing to calcuate, zero number of centroids\n");
    return 0;
  }
  if (vec_dim <= 0){
    printf("The input data has null dimension, nothing can be done, exiting the program\n");
    exit(0);
  }
  if (No_data_pts <= 0){
    printf("There are no data points to process\n");
    exit(0);
  }
  int N = No_data_pts;
  int d = vec_dim;
  int cl_near_dat[k]; // to count the number of nearest data points
  double new_cl[k][d]; // to store the new centroids
  for (int l = 0; l < k; l++){
    cl_near_dat[l] = 0;
    for (int j = 0; j < d; j++){
      new_cl[l][j] = 0.0;
    }
  }
  int z;
  for (int i = 0; i < N; i++){
    z = argmin[i] - 1;
    cl_near_dat[z] += 1;
    for (int j = 0; j < d; j++){
      new_cl[z][j] += *(data[i] + j);
    }
  }
  for (int l = 0; l < k; l++){
    if (cl_near_dat[l] != 0){
      for (int j = 0; j < d; j++){
        *(NewCentroid[l] + j) = new_cl[l][j] / cl_near_dat[l];
      }
    } else{
      for (int j1 = 0; j1 < d; j1++){
        *(NewCentroid[l] + j1) = *(Centroid[l] + j1);
      }
    }
  }
  return 0;
}


