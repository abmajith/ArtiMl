#ifndef K_MEAN_CLUSTERING
#define K_MEAN_CLUSTERING

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>

int KMeanCl_Normalize_data_Db(int No_Data_pts, int vec_dim, double *data[], double Min_range[], double Max_range[], double Abs_max[]);

int KMeanCl_Initialize_Centroids_Db(int k, int vec_dim, double Min_range[], double Max_range[], double Abs_max[], double *Centroid[]);

int KMeanCl_minarg_to_centroid_Db(int k, int vec_dim, int No_data_pts, double *data[], double *Centroid[], int argmin[]);

int KMeanCl_update_centroids_Db(int k, int vec_dim, int No_data_pts, double *data[], double *Centroid[], double *NewCentroid[], int argmin[] );

#endif

