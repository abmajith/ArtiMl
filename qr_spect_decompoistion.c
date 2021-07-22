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

  double UTR[k][len];
  double UBR[len][len];
  /* check to delete once the job is done to reuse the memeory */
  /* computing top right block */
  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len; z++)
        init += *(*(ptrU + k + z ) + j )   *  Q[z][j];
      UTR[i][j] = init;
      }
  /* computing bottom right block */
  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len; z++)
        init +=  *(*(ptrU + k + z ) + j + k )    * Q[z][j];
      UBR[i][j] = init;
    }
  /* copying the computed blocks to the memory of ptrQ */
  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++)
      *(*(ptrU + k + j ) + i) = UTR[i][j];
  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++)
      *(*(ptrU + k + j ) + i + k ) = UBR[i][j];
  /*Matrix U is updated properly, it is double checked! */

  double H2Q[k][len];
  double QH3[len][k];
  double QH4Qi[len][len];
  double QH4Q[len][len];

  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++){
      init = 0.0;
      for (int z = 0; z < len)
        init += *(*(ptrH + k + z) + i)   *  Q[z][j];
      H2Q[i][j] = init;
    }
  for (int i = 0; i < k; i++)
    for (int j = 0; j < len; j++)
      *(*(ptrH + k + j) + i ) = H2Q[i][j];

  for (int i = 0; i < len; i++)
    for (int j = 0; j < k; j++){
      init = 0.0;
      for (int z = 0; z < len)
        init += Q[i][z]  * *(*(ptrH + j) + k + z);
      QH3[i][j] = init;
    }
  for (int i = 0; i < len; i++)
    for (int j = 0; j < k; j++)
      *(*(ptrH + j ) + k + i) = QH3[i][j];

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
}

/*lot of memeory optimization or reduce the memory conception can be done, so try later, now stop */
