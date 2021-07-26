#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* we construct the hessenberg form for a given square matrix by householder reduction */

/*H -> H' = QHQ, U -> U' = UQ*/
/* matrix H is m by m, U is m by m, Q is m by m, constructed in the same way as householder reflection!*/

int Hess_Itr_QRmthd_Db(double **ptrH, double **ptrU, int m, int k){
  if (m <= 2 | k > (m-2) | k < 1)
    return 4; /* dimenstion mismatch, nothing to do*/
  int len = m - k; /* dimenstion of submatrix Q*/
  double init = 0.0;
  double pivot; /*Rkminus1k */
  double *ptr; /* to work with column [Rk-1k Rkk Rk+2k ...Rmk]T */
  ptr = *(ptrH + k - 1) + k;
  pivot = *(*(ptrH + k - 1) + k);
  double sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (double) Sq_Normtwo_Db(ptr, m - k);
  alpha = (double) sqrt(sigma);

  /* if calculated norm is same as the modulus of pivot then divide by zero occurs in the rotation */

  double nf; /*to  normalization factor */
  /* if sigma is zero then that column is zero nothing to do for that column */
  if (sigma == 0.0)
    return 3;

  /* if calculated norm is same as the pivot then divide by zero occurs in the reflection so */
  if (alpha == pivot)
    alpha = 0.0 - alpha;

  nf = sigma - alpha * pivot;
  double Q[len][len]; /* to store the rotated submatrix */
  /* computing the Q matrix first Ilen*len - v outer product v / nf*/
  Q[0][0] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
  for (int i = 1; i < len; i++) /*filling diagonal element */
    Q[i][i] = 1.0 - (  *(*(ptrH + k - 1) + k + i)   *  *(*(ptrH + k - 1) + k + i) ) / nf;
  for (int i = 1; i < len; i++) /* filling 0th column and 0th row*/
    Q[i][0] = Q[0][i] =  0.0 - (  ( pivot - alpha )  *  *(*(ptrH + k - 1) + k + i) ) / nf;
  for (int i = 1; i < len; i++)
    for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
      Q[j][i] = Q[i][j] = 0.0 - (  *(*(ptrH + k - 1) + k + i)  *   *(*(ptrH + k - 1) + k + j) ) / nf;
    } /* outer product of u1 and u1T is computed and stored in matrix Q*/
    /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/





  if (k > 1){
    double **ptrUTR;
    ptrUTR = (double **) malloc(len * sizeof(double*));
    for (int i = 0; i < len; i++)
      *(ptrUTR + i) = (double *) malloc(k * sizeof(double));

    /* check to delete once the job is done to reuse the memeory */
    /* computing top right block */
    for (int i = 0; i < k; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += *(*(ptrU + k + z ) + j )   *  Q[z][j];
        *(*(ptrUTR + j) + i) = init;
        }
    /* copying the computed blocks to the memory of ptrQ */
    for (int i = 0; i < k; i++)
      for (int j = 0; j < len; j++)
        *(*(ptrU + k + j ) + i) = *(*(ptrUTR + j) + i); 

    for (i = 0; i < len; i++)
      free( *(ptrUTR + i));
    free(ptrUTR);

    double **ptrUBR;
    ptrUBR = (double **) malloc (len * sizeof(double*));
    for (int i = 0; i < len; i++)
      *(ptrUBR + i) = (double *) malloc(len * sizeof(double));

  
    /* computing bottom right block */
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init +=  *(*(ptrU + k + z ) + j + k )    * Q[z][j];
        *(*(ptrUBR + j) + i) = init;
      }
    /* copying the computed blocks to the memory of ptrQ */
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++)
        *(*(ptrU + k + j ) + i + k ) = *(*(ptrUBR + j) + i);
    /*Matrix U is updated properly, it is double checked! */
  
    for (int i = 0; i < len; i++)
      free( *(ptrUBR + i));
    free(ptrUBR);
  }
  else
  {
   **ptrU = 1.0;
   for (int i = 0; i < m; i++)
     *(*(ptrU + i)) = *(*ptrU + i) = 0.0;
   for (int i = 1; i < m; i++)
    for (int j = 1; j <= i ; j++)
     *(*(ptrU + j) + i) = *(*(ptrU + i) + j) = Q[i-1][j-1]; 
  }






  double **ptrH2Q;
  ptrH2Q = (double **) malloc(len * sizeof(double*));
  for (int i = 0; i < len; i++)
    *(ptrH2Q + i) = (double *) malloc(k * sizeof(double));

  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len)
        init += *(*(ptrH + k + z) + i)   *  Q[z][j];
      *(*(ptrH2Q + j) + i) = init;
    }
  /* copying to the H matrix*/
  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++)
      *(*(ptrH + k + j) + i) = *(*(ptrH2Q + j) + i);

  for (int i = 0; i < len; i++)
    free( *(ptrH2Q + i) );
  free( ptrH2Q );


  double **ptrQH3;
  ptrQH3 = (double **) malloc(k * sizeof(double*));
  for (int i = 0; i < k; i++)
    *(ptrQH3 + i) = (double *) malloc(len * sizeof(double));

  for (int i = (0 + k - 1); i < len; i++) /*reremember here the column until k - 1 is zero for H3 */
    for (int j = 0; j < k; j++){
      init = 0.0;
      for (int z = 0; z < len)
        init += Q[i][z]  * *(*(ptrH + j) + k + z);
      *(*(ptrQH3 + j) + i) = init;
    }
  for (int i = (0 + k - 1); i < len; i++)
    for (int j = 0; j < k; j++)
      *(*(ptrH + j ) + k + i) = *(*(ptrQH3 + j) + i);

  for (int i = 0; i < k; i++)
    free( *(ptrQH3 + i) );
  free(ptrQH3);


  double QH4Qi[len][len];
  double QH4Q[len][len];

  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len; z++)
        init +=  *(*(ptrH + k + z) + k + i )   *  Q[z][j];
      QH4Qi[i][j] = init;
    }
  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len; z++)
        init += Q[i][z] * QH4Qi[z][j];
      QH4Q[i][j]  = init;
    }

  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++)
     *(*(ptrH + k + j) + k + i)  = QH4Q[i][j];
  return 0;
}/*single rotation to form the intermediate hessenberg matrix is done for double type this code algorithm is correct just check the proper indices are used are mixed up with different indices! */


int Hess_Itr_QRmthd_Fl(float **ptrH, float **ptrU, int m, int k){
  if (m <= 2 | k > (m-2) | k < 1)
    return 4; /* dimenstion mismatch, nothing to do*/
  int len = m - k; /* dimenstion of submatrix Q*/
  float init = 0.0;
  float pivot; /*Rkminus1k */
  float *ptr; /* to work with column [Rk-1k Rkk Rk+2k ...Rmk]T */
  ptr = *(ptrH + k - 1) + k;
  pivot = *(*(ptrH + k - 1) + k);
  float sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (float) Sq_Normtwo_Fl(ptr, m - k);
  alpha = (float) sqrt(sigma);

  /* if calculated norm is same as the modulus of pivot then divide by zero occurs in the rotation */

  float nf; /*to  normalization factor */
  /* if sigma is zero then that column is zero nothing to do for that column */
  if (sigma == 0.0)
    return 3;

  /* if calculated norm is same as the pivot then divide by zero occurs in the reflection so */
  if (alpha == pivot)
    alpha = 0.0 - alpha;

  nf = sigma - alpha * pivot;
  float Q[len][len]; /* to store the rotated submatrix */
  /* computing the Q matrix first Ilen*len - v outer product v / nf*/
  Q[0][0] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
  for (int i = 1; i < len; i++) /*filling diagonal element */
    Q[i][i] = 1.0 - (  *(*(ptrH + k - 1) + k + i)   *  *(*(ptrH + k - 1) + k + i) ) / nf;
  for (int i = 1; i < len; i++) /* filling 0th column and 0th row*/
    Q[i][0] = Q[0][i] =  0.0 - (  ( pivot - alpha )  *  *(*(ptrH + k - 1) + k + i) ) / nf;
  for (int i = 1; i < len; i++)
    for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
      Q[j][i] = Q[i][j] = 0.0 - (  *(*(ptrH + k - 1) + k + i)  *   *(*(ptrH + k - 1) + k + j) ) / nf;
    } /* outer product of u1 and u1T is computed and stored in matrix Q*/
    /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/


  if ( k > 1){
    float **ptrUTR;
    ptrUTR = (float **) malloc(len * sizeof(float*));
    for (int i = 0; i < len; i++)
      *(ptrUTR + i) = (float *) malloc(k * sizeof(float));

    /* check to delete once the job is done to reuse the memeory */
    /* computing top right block */
    for (int i = 0; i < k; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += *(*(ptrU + k + z ) + j )   *  Q[z][j];
        *(*(ptrUTR + j) + i) = init;
        }
    /* copying the computed blocks to the memory of ptrQ */
    for (int i = 0; i < k; i++)
      for (int j = 0; j < len; j++)
        *(*(ptrU + k + j ) + i) = *(*(ptrUTR + j) + i);

    for (i = 0; i < len; i++)
      free( *(ptrUTR + i));
    free(ptrUTR);

    float **ptrUBR;
    ptrUBR = (float **) malloc (len * sizeof(float*));
    for (int i = 0; i < len; i++)
      *(ptrUBR + i) = (float *) malloc(len * sizeof(float));


    /* computing bottom right block */
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init +=  *(*(ptrU + k + z ) + j + k )    * Q[z][j];
        *(*(ptrUBR + j) + i) = init;
      }
    /* copying the computed blocks to the memory of ptrQ */
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++)
        *(*(ptrU + k + j ) + i + k ) = *(*(ptrUBR + j) + i);
    /*Matrix U is updated properly, it is double checked! */

    for (int i = 0; i < len; i++)
      free( *(ptrUBR + i));
    free(ptrUBR);
  }
  else
  {
   **ptrU = 1.0;
   for (int i = 0; i < m; i++)
     *(*(ptrU + i)) = *(*ptrU + i) = 0.0;
   for (int i = 1; i < m; i++)
    for (int j = 1; j <= i ; j++)
     *(*(ptrU + j) + i) = *(*(ptrU + i) + j) = Q[i-1][j-1];
  }


  float **ptrH2Q;
  ptrH2Q = (float **) malloc(len * sizeof(float*));
  for (int i = 0; i < len; i++)
    *(ptrH2Q + i) = (float *) malloc(k * sizeof(float));

  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len)
        init += *(*(ptrH + k + z) + i)   *  Q[z][j];
      *(*(ptrH2Q + j) + i) = init;
    }
  /* copying to the H matrix*/
  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++)
      *(*(ptrH + k + j) + i) = *(*(ptrH2Q + j) + i);

  for (int i = 0; i < len; i++)
    free( *(ptrH2Q + i) );
  free( ptrH2Q );


  float **ptrQH3;
  ptrQH3 = (float **) malloc(k * sizeof(float*));
  for (int i = 0; i < k; i++)
    *(ptrQH3 + i) = (float *) malloc(len * sizeof(float));

  for (int i = (0 + k - 1); i < len; i++)
    for (int j = 0; j < k; j++){
      init = 0.0;
      for (int z = 0; z < len)
        init += Q[i][z]  * *(*(ptrH + j) + k + z);
      *(*(ptrQH3 + j) + i) = init;
    }
  for (int i = (0 + k - 1); i < len; i++)
    for (int j = 0; j < k; j++)
      *(*(ptrH + j ) + k + i) = *(*(ptrQH3 + j) + i);

  for (int i = 0; i < k; i++)
    free( *(ptrQH3 + i) );
  free(ptrQH3);

  float QH4Qi[len][len];
  float QH4Q[len][len];

  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len; z++)
        init +=  *(*(ptrH + k + z) + k + i )   *  Q[z][j];
      QH4Qi[i][j] = init;
    }
  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len; z++)
        init += Q[i][z] * QH4Qi[z][j];
      QH4Q[i][j]  = init;
    }

  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++)
     *(*(ptrH + k + j) + k + i)  = QH4Q[i][j];
  return 0;
}/*single rotation to form the intermediate hessenberg matrix is done for double type this code algorithm is correct just check the proper indices are used are mixed up with different indices! */

/* House Holder reflection to compute intermediate hessenberg matrix H = QR for the given matrix H  and then find the similar matrix to H' = QHQ keep
 * doing this until finding full hessenberg matrix*/

/* once findinding the full hessenberg similar matrix to the original matrix, now perform QR operation such that H = QR, H1 Q1R1 = = QHQ, H2 = Q1QHQQ1
 * so on.....*/

int QR_HessMat_Db(double **ptrH, int m){
  /* generate Q so that H = QR, new H = R and apply RQ and store in H the output H is similar to input H */
  /* Q = GnG2...G1 */
  /* G1 = [ [c,s;-s,c] 0; 0, I]*/
  /* G2 = [1, 0; vec0, [c,s;-s,c] vec0; vec0t I]*/
  /* so on*/
  /* H1=Gn....G2G1H ,  final barH = H1G1*G2*.....Gn* */
  double *ptrcs;
  ptrcs = (double *) malloc( 2 * (m - 1) * sizeof(double) );
  double sqnrm;
  double nrm;
  double a1,a2;
  for (int k = 0; k < m-1; k++){
    sqnrm = (double) Sq_Normtwo_Db( (*(ptrH + k) + k)  , 2);
    if (sqnrm == 0.0)
      return 2; /*indicating the divide by zero */
    nrm = (double) sqrt(sqnrm);
    *(ptrcs + (2 * k) ) = *(*(ptrH + k) + k) / nrm;
    *(ptrcs + (2 * k) + 1 ) = *(*(ptrH + k) + k + 1) / nrm;
    /*to form [c,s; -s,c] */
    for (j = k; j < m; j++){
      a1 = *(*(ptrH + j) + k);
      a2 = *(*(ptrH + j) + k + 1);
      *(*(ptrH + j) + k) = *(ptrcs + (2 * k) )  * a1 + *(ptrcs + (2 * k) + 1 ) * a2;
      *(*(ptrH + j) + k + 1) = *(ptrcs + (2 * k) ) * a2 - *(ptrcs + (2 * k) + 1 ) * a1;
    }
  } /* now ptrH is the R*/
  /* we have to do now R * G1*G2*G3*...Gn* */

  for (int k = 0; k < m - 1; k++){
    for (int j = 0; j < k + 2; j++){
      a1 = *(*(ptrH + k) + j);
      a2 = *(*(ptrH + k + 1) + j);
      *(*(ptrH + k) + j ) = a1  *  *(ptrcs + (2 * k))  +  a2  *  *(ptrcs + (2 * k) + 1);
      *(*(ptrH + k + 1) + j) = a2 *  *(ptrcs + (2 * k)) - a1  *  *(ptrcs + (2 * k) + 1);
    }
  }
  return 0;
}


int QR_HessMat_Fl(float **ptrH, int m){
  /* generate Q so that H = QR, new H = R and apply RQ and store in H the output H is similar to input H */
  /* Q = GnG2...G1 */
  /* G1 = [ [c,s;-s,c] 0; 0, I]*/
  /* G2 = [1, 0; vec0, [c,s;-s,c] vec0; vec0t I]*/
  /* so on*/
  /* H1=Gn....G2G1H ,  final barH = H1G1*G2*.....Gn* */
  float *ptrcs;
  ptrcs = (float *) malloc( 2 * (m - 1) * sizeof(float) );
  float sqnrm;
  float nrm;
  float a1,a2;
  for (int k = 0; k < m-1; k++){
    sqnrm = (float) Sq_Normtwo_Db( (*(ptrH + k) + k)  , 2);
    if (sqnrm == 0.0)
      return 2; /*indicating the divide by zero */
    nrm = (float) sqrt(sqnrm);
    *(ptrcs + (2 * k) ) = *(*(ptrH + k) + k) / nrm;
    *(ptrcs + (2 * k) + 1 ) = *(*(ptrH + k) + k + 1) / nrm;
    /*to form [c,s; -s,c] */
    for (j = k; j < m; j++){
      a1 = *(*(ptrH + j) + k);
      a2 = *(*(ptrH + j) + k + 1);
      *(*(ptrH + j) + k) = *(ptrcs + (2 * k) )  * a1 + *(ptrcs + (2 * k) + 1 ) * a2;
      *(*(ptrH + j) + k + 1) = *(ptrcs + (2 * k) ) * a2 - *(ptrcs + (2 * k) + 1 ) * a1;
    }
  } /* now ptrH is the R*/
  /* we have to do now R * G1*G2*G3*...Gn* */

  for (int k = 0; k < m - 1; k++){
    for (int j = 0; j < k + 2; j++){
      a1 = *(*(ptrH + k) + j);
      a2 = *(*(ptrH + k + 1) + j);
      *(*(ptrH + k) + j ) = a1  *  *(ptrcs + (2 * k))  +  a2  *  *(ptrcs + (2 * k) + 1);
      *(*(ptrH + k + 1) + j) = a2 *  *(ptrcs + (2 * k)) - a1  *  *(ptrcs + (2 * k) + 1);
    }
  }
  return 0;
}


