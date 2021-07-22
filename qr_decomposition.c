#include <math.h> /*i am using sqrt function, so its needed */



/* Square of Norm Two*/
/* Finding the square of norm two of given pointer to an array of size m*/
/*norm two: square root of sum of square of elements from the array */
/* note we are returning the square of norm two*/
double Sq_Normtwo_Db(double *col, int m){
  /* col is a pointer to array of size m whose data type is double */
  double x = 0.0;
  for (int i = 0; i < m; i++)
    x += *(col + i)  *  *(col + i);
  return x; /*this function return the square of norm two of given array of size m */
} /* This code was double checked */

float Sq_Normtwo_Fl(float *col, int m){
  float x = 0.0;
  for (int i = 0; i < m; i++)
    x += *(col + i)  *  *(col + i);
  return x;
} /* This code is same as previous one, but float type input*/

int Sq_Normtwo_Int(int *col, int m){
  long x = 0;
  for (int i = 0; i < m; i++)
    x += *(col + i)  *  *(col + i);
  return x;
} /* This code is same as previous one, but integer type*/

long Sq_Normtwo_Ln(long *col, int m){
  long x = 0;
  for (int i = 0; i < m; i++)
    x += *(col + i)  *  *(col + i);
  return x;
} /* This code is same as previous one, but long type*/

/* Note: for real programs, do check that the input and output can hold the same data type otherwise we need to cast the output to higher data type */

/*-----square of norm two of given array function is double checked, remember the note:----------- */

/*-------------------------------------------------------------------------------------------------------- */





/* here the pointer to pointer of a matrix storage method*/
/* n by m matrix*/
/*ptrA1 ptrA1 ptrA2 ..... prAm */
/*a11    a12   a13         a1m */
/*a21    a22   a23         a2m  */
/*....               .......... */
/*an1    an2   an3         anm */



/*House Holder Reflection*/
/* Iteration of House holder reflection on the kth column*/
int HHdRl_Sqmat_Itr_Db(double **ptrR, int m, int k1 ){ /* R is the matrix of m cross m, and relect on the ptrRk1 th array*/
  /*k1 should be non zero positive value */
  /*ptrR is the pointer to pointers of m arrays, each array has size m */
  if ( k1 >= m || k1 < 1 ) /* if k1 is less then 1 or greater than or equal to m then nothing to do */
    return 2; /* indicating dimention mismatch*/

  int k = k1 - 1; /* cprogram convention for the index*/
  int len = m - k; /* dimension of the submatrix and this is correct, checked multiple times*/
  double init = 0.0;
  double pivot = *(*(ptrR + k) + k); /*Rkk */

  double *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */
  ptr = *(ptrR + k) + k; /* address of Rk,k */

  double sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (double) Sq_Normtwo_Db(ptr, len);
  alpha = (double) sqrt(sigma);
  if (sigma == 0.0) /* in the future for stable result have to introduce some epsilon value instead of zero*/
    return 3; /* indicating singular nature of the matrix*/

  /* if calculated norm is same as the pivot then divide by zero occurs in the rotation so */
  if (alpha == pivot)
    alpha = 0.0 - alpha;/*changing the sign to avoid divide by zero case */

  double nf; /* normalization factor to compute Q for the len by len submatrix*/
  nf = sigma - alpha * pivot;
  double Q[len][len], BR[len][len]; /* to store the rotated submatrix of R for making the k+1K+2K..m to zero */

  /*Q = Identity matrix len by  len - (1over nf ) outer  product of u1 and u1T */
  Q[0][0] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
  for (int i = 1; i < len; i++) /*filling diagonal element */
    Q[i][i] = 1.0 - (  *(*(ptrR + k) + k + i)   *  *(*(ptrR + k) + k + i) ) / nf;
  for (int i = 1; i < len; i++) /* filling 0th column and 0th row*/
    Q[i][0] = Q[0][i] =  0.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + i) ) / nf;
  for (int i = 1; i < len; i++) 
    for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
      Q[j][i] = Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)  *   *(*(ptrR + k) + k + j) ) / nf;
    } /* outer product of u1 and u1T is computed and stored in matrix Q*/
  /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/


  /* Rnew = Q'R = [I,zero; zero, Q ] * R ([TLR, TRR; BLR, BRR]) = [TLR, TRR; Q*BLR, Q* BRR] */
  /*Q is len by len size, BRR is len by len size, BLR is len by m-len size, TRR is m-len by len size and TLR is m-len by m-len size */
  /* block matrix multiplication scheme*/
  /* let compute, BR = Q * BRR */
  for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j + k) + z + k);
        BR[i][j] = init;
      } /* bottom right of new R is computed*/

  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++)
      *(*(ptrR + j + k) + i + k) = BR[i][j]; /*copying the new BR of R to existing R (memory) */

    /* now we have to compute the new R Bottom left submatrix*/
    /* if k is zero, we dont have to so*/
  if ( k != 0)
  {
    /*let compute, BL = Q * BLR */
    double BL[len][m - len];
    for (int i = 0; i < len ; i++)
      for (int j = 0; j < (m - len); j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j) + z + k);
        BL[i][j] = init;
      } /* Bottom left of new R is computed*/
    for (int i = 0; i < len; i++)
      for (int j = 0; j < (m - len); j++)
        *(*(ptrR + j) + i + k) = BL[i][j]; /*copying the new BL of R to exisitng R (memory) */
  }
  return 0; /* indicating that succesfully computed!*/

}/* This code is double checked!*/



/*House Holder Reflection*/
/* Iteration of House holder reflection on the kth column*/
int HHdRl_Sqmat_Itr_Fl(float **ptrR, int m, int k1 ){ /* R is the matrix of m cross m, and relect on the ptrRk1 th array*/
  /*k1 should be non zero positive value */
  /*ptrR is the pointer to pointers of m arrays, each array has size m */
  if ( k1 >= m || k1 < 1 )
    return 2; /* indicating dimension mismatch*/

  int k = k1 - 1; /* cprogram convention for the index*/
  int len = m - k; /* dimension of the submatrix and this is correct checked multiple times*/
  float init = 0.0;
  float pivot = *(*(ptrR + k) + k); /*Rkk */

  float *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */
  ptr = *(ptrR + k) + k; /* address of Rk,k */

  float sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (float) Sq_Normtwo_Db(ptr, len);
  alpha = (float) sqrt(sigma);
  
  if (sigma == 0.0) /* in the future for stable result have to introduce some epsilon value instead of zero*/
    return 3; /* indicating singular nature of the matrix*/

  /* if calculated norm is same as the pivot then divide by zero occurs in the rotation so */
  if (alpha == pivot)
    alpha = 0.0 - alpha;

  float nf; /* normalization factor to compute Q for the len by len submatrix*/
  nf = sigma - alpha * pivot;
  float Q[len][len], BR[len][len]; /* to store the rotated submatrix of R for making the k+1K+2K..m to zero */

  /*Q = Identity matrix len by  len - (1over nf ) outer  product of u1 and u1T */
  Q[0][0] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
  for (int i = 1; i < len; i++) /*filling diagonal element */
    Q[i][i] = 1.0 - (  *(*(ptrR + k) + k + i)   *  *(*(ptrR + k) + k + i) ) / nf;
  for (int i = 1; i < len; i++) /* filling 0th column and 0th row*/
    Q[i][0] = Q[0][i] =  0.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + i) ) / nf;
  for (int i = 1; i < len; i++)
    for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
      Q[j][i] = Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)  *   *(*(ptrR + k) + k + j) ) / nf;
    } /* outer product of u1 and u1T is computed and stored in matrix Q*/
  /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/

  /* Rnew = Q'R = [I,zero; zero, Q ] * R ([TLR, TRR; BLR, BRR]) = [TLR, TRR; Q*BLR, Q* BRR] */
  /*Q is len by len size, BRR is len by len size, BLR is len by m-len size, TRR is m-len by len size and TLR is m-len by m-len size */
  /* block matrix multiplication scheme*/
  /* let compute, BR = Q * BRR */
  for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j + k) + z + k);
        BR[i][j] = init;
      } /* bottom right of new R is computed*/
  for (int i = 0; i < len; i++)
    for (int j = 0; j < len; j++)
      *(*(ptrR + j + k) + i + k) = BR[i][j]; /*copying the new BR of R to existing R (memory) */

    /* now we have to compute the new R Bottom left submatrix*/
    /* if k is zero, we dont have to so*/
  if ( k != 0)
  {
    /*let compute, BL = Q * BLR */
    float BL[len][m - len];
    for (int i = 0; i < len ; i++)
      for (int j = 0; j < (m - len); j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j) + z + k);
        BL[i][j] = init;
      } /* Bottom left of new R is computed*/
    for (int i = 0; i < len; i++)
      for (int j = 0; j < (m - len); j++)
        *(*(ptrR + j) + i + k) = BL[i][j]; /*copying the new BL of R to exisitng R (memory) */
  }
  return 0;
}/* This code is double checked!*/
/* the above code is same as before, but for float type!*/


/* House holder reflection on the given square matrix iteration was done */

/*------------------------------------------------------------------------------------- */




/*Given matrix A of size m by n, find Q and R decomposition of A*/
/* size of Q is m by m*/


/* need to review */
int QR_HH_Itr_Db (double **ptrA, double **ptrQ, double **ptrR, int m, int n, int k1){
  /* ptrA is the pointer to pointer of n arrays, each arrays has size m */
  /* ptrQ is the pointer to pointer of m arrays, each arrays has size m */
  /* ptrR is the pointer to pointer of n arrays, each arrays has size m*/

  if (m <= n)
    if (k1 >= m )
      return 2; /*dimension mismatch */
  if (n < m)
    if (k1 > n)
      return 2;
  if (k1 < 1)
    return 2;


  int k = k1 - 1;
  int len = m - k; /* dimenstion of submatrix Q*/
  double init = 0.0;
  double pivot; /*Rkk */

  double *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */

  double sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  /* if calculated norm is same as the modulus of pivot then divide by zero occurs in the rotation */
  double nf; /*to  normalization factor */

  if (k == 0){ /*reflection of *ptrA 0th column */
    pivot = **ptrA;
    ptr = *ptrA;
    sigma = (double) Sq_Normtwo_Db(ptr, len);
    alpha = (double) sqrt(sigma);
    /* if sigma is zero then that column is zero nothing to do for that column */
    if (sigma == 0.0)
      return 3;

   /* if calculated norm is same as the pivot then divide by zero occurs in the reflection so */
    if (alpha == pivot)
      alpha = 0.0 - alpha;

    nf = sigma - alpha * pivot;
    /* computing the Q matrix first Im*m - v outer product v / nf*/
    *(*(ptrQ + 0) + 0)  = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
    for (int i = 1; i < m; i++) /*filling diagonal element */
      *(*(ptrQ + i) + i)  = 1.0 - (  *(*ptrA + i)   *  *(*ptrA + i) ) / nf;
    for (int i = 1; i < m; i++) /* filling 0th column and 0th row*/
      *(*(ptr + i)) =  *(*ptrQ + i)  =  0.0 - (  ( pivot - alpha )  *  *(*ptrA + i) ) / nf;
    for (int i = 1; i < m; i++)
      for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
        *(*(ptrQ + j) + i) = *(*(ptrQ + i) + j) = 0.0 - (  *(*ptrA + i)  *   *(*ptrA + j) ) / nf;
      } /* outer product of u1 and u1T is computed and stored in matrix Q*/
    /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/
    /* computing R matrix = QA*/
    /* Q is m by m and A is m by n*/
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++){
        init = 0.0;
        for (int z = 0; z < m; z++)
          init += *(*(ptrQ + z) + i)  *  *(*(ptrA + j) + z);
        *(*(ptrR + j) + i) = init;
      } /* computed R = QA and stored in the memory of R i.e ptrR*/
  }/* for the reflection of 0th column is done i.e Q and R is computed*/

  if ( k != 0){
    pivot = *(*(ptrR + k ) + k );
    ptr = *(ptrR + k) + k;
    sigma = (double) Sq_Normtwo_Db(ptr, len);
    alpha = (double) sqrt(sigma);
    /* if the sigma is zero, then nothing to do there */

    if (sigma == 0.0)
      return 3;

   /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      alpha = 0.0 - alpha;

    nf = sigma - alpha * pivot;
    double Q[len][len]; /* to store the rotated submatrix of R for making the k+1K+2K..m to zero */
    /* computing the Q matrix first Ilen*len - v outer product v / nf*/
    Q[0][0] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
    for (int i = 1; i < len; i++) /*filling diagonal element */
      Q[i][i] = 1.0 - (  *(*(ptrR + k) + k + i)   *  *(*(ptrR + k) + k + i) ) / nf;
    for (int i = 1; i < len; i++) /* filling 0th column and 0th row*/
      Q[i][0] = Q[0][i] =  0.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + i) ) / nf;
    for (int i = 1; i < len; i++)
      for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
        Q[j][i] = Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)  *   *(*(ptrR + k) + k + j) ) / nf;
      } /* outer product of u1 and u1T is computed and stored in matrix Q*/
      /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/

    /* Q' = [I, zero; zero, Q]*/
    /* here I is m -len by m - len */
    /* after each refelection Q = Qold Q'*/
    /* Qold = [a, b; c, d]*/
    /* where a is m - len by m - len, d is len by len, c is len by m - len and b is m - len by len*/
    /* updating new  Q =  [a, b; c, d] [I,zero; zero, Q] = [a, bQ; c, dQ]*/
    /* bQ is m - len by m (i.e k by m), store bQ in QTR */
    /* dQ is len by len, store dQ in QBR */
    /* block multiplication */
    double QTR[k][len];
    double QBR[len][len];
    /* computing top right block */
    for (int i = 0; i < k; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += *(*(ptrQ + k + z ) + j )   *  Q[z][j];
        QTR[i][j] = init;
        }
      /* computing bottom right block */
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init +=  *(*(ptrQ + k + z ) + j + k )    * Q[z][j];
        QBR[i][j] = init;
      }

      /* copying the computed blocks to the memory of ptrQ */
      for (int i = 0; i < k; i++)
        for (int j = 0; j < len; j++)
          *(*(ptrQ + k + j ) + i) = QTR[i][j];
      for (int i = 0; i < len; i++)
        for (int j = 0; j < len; j++)
          *(*(ptrQ + k + j ) + i + k ) = QBR[i][j];
      /*Matrix Q is updated properly, it is double checked! */

      /* now updating R*/
      /* Easiest approach is ptrR = ptrQ times ptrA, but we can reduce the memory conception and computation time by block matrix way */
      /* So, new R = Q' times Rold, where Rold is m by n matrix*/
      /* Recall Q' = [I, zero; zero, Q] m by m matrix */
      /* R = Q' Rold = [I, zero; zero, Q] [a, b; c, d] = [a,b; Qc, Qd]*/
      /* a is m - len by m - len, c is len by m - len, b is n - m + len by m - len and d is len by n - m + len*/
      /* where Qc is len by m - len, store Qc in BLR*/
      /* Qd is len by n - m + len, store Qd in BRR*/
      double BLR[len][k];
      double BRR[len][n-k];         
      /* calculating the BLR */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < k; j++){
          init = 0.0;
          for (int z = 0; z < len; z++)
            init += Q[i][z] * *(*(ptrR + j) + k + z);
          BLR[i][j] = init;
        }
      /* calculating BRR */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < (n-k) j++){
          init = 0.0;
          for (int z = 0; z < len; z++)
            init += Q[i][z] * *(*(ptrR + k + j) + k + z);
          BRR[i][j] = init;
        }
      /* copying the computed blocks to the memory of ptrR */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < k; j++)
          *(*(ptrR + j) + k + i) = BLR[i][j];
      for (int i = 0; i < len; i++)
        for (int j = 0; j < (n-k); j++)
          *(*(ptrR + k + j) + k + i) = BRR[i][j];
      /* matrix R is updated properly, it is double checked!*/
  }
  return 0;
}

int QR_HH_Itr_Fl (float **ptrA, float **ptrQ, float **ptrR, int m, int n, int k1){
  /* ptrA is the pointer to pointer of n arrays, each arrays has size m */
  /* ptrQ is the pointer to pointer of m arrays, each arrays has size m */
  /* ptrR is the pointer to pointer of n arrays, each arrays has size m*/

  if (m <= n)
    if (k1 >= m )
      return 2;
  if (n < m)
    if (k1 > n)
      return 2;
  if (k1 <= 0)
    return 2;

  int k = k1 - 1;
  int len = m - k; /* dimenstion of submatrix Q*/
  float init = 0.0;
  float pivot; /*Rkk */

  float *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */

  float sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  /* if calculated norm is same as the modulus of pivot then divide by zero occurs in the rotation */
  float nf; /*to  normalization factor */

  if (k == 0){ /*reflection of *ptrA 0th column */
    pivot = **ptrA;
    ptr = *ptrA;
    sigma = (float) Sq_Normtwo_Fl(ptr, len);
    alpha = (float) sqrt(sigma);
    /*If the sigma is zero then nothing to do for that column */
    if (sigma == 0.0)
      return 3;

   /* if calculated norm is same as the pivot then divide by zero occurs in the reflection */
    if (alpha == pivot)
      alpha = 0.0 - alpha;

    nf = sigma - alpha * pivot;
    /* computing the Q matrix first Im*m - v outer product v / nf*/
    *(*(ptrQ + 0) + 0)  = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
    for (int i = 1; i < m; i++) /*filling diagonal element */
      *(*(ptrQ + i) + i)  = 1.0 - (  *(*ptrA + i)   *  *(*ptrA + i) ) / nf;
    for (int i = 1; i < m; i++) /* filling 0th column and 0th row*/
      *(*(ptr + i)) =  *(*ptrQ + i)  =  0.0 - (  ( pivot - alpha )  *  *(*ptrA + i) ) / nf;
    for (int i = 1; i < m; i++)
      for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
        *(*(ptrQ + j) + i) = *(*(ptrQ + i) + j) = 0.0 - (  *(*ptrA + i)  *   *(*ptrA + j) ) / nf;
      } /* outer product of u1 and u1T is computed and stored in matrix Q*/
    /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/
    /* computing R matrix = QA*/
    /* Q is m by m and A is m by n*/
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++){
        init = 0.0;
        for (int z = 0; z < m; z++)
          init += *(*(ptrQ + z) + i)  *  *(*(ptrA + j) + z);
        *(*(ptrR + j) + i) = init;
      } /* computed R = QA and stored in the memory of R i.e ptrR*/
  }/* for the reflection of 0th column is done i.e Q and R is computed*/

  if ( k != 0){
    pivot = *(*(ptrR + k ) + k );
    ptr = *(ptrR + k) + k;
    sigma = (float) Sq_Normtwo_Fl(ptr, len);
    alpha = (float) sqrt(sigma);
    /* if the sigma is zero, then nothing to for that column*/
    if (sigma == 0.0)
      return 3;

   /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      alpha = 0.0 - alpha;

    nf = sigma - alpha * pivot;
    float Q[len][len];
    /* computing the Q matrix first Ilen*len - v outer product v / nf*/
    Q[0][0] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
    for (int i = 1; i < len; i++) /*filling diagonal element */
      Q[i][i] = 1.0 - (  *(*(ptrR + k) + k + i)   *  *(*(ptrR + k) + k + i) ) / nf;
    for (int i = 1; i < len; i++) /* filling 0th column and 0th row*/
      Q[i][0] = Q[0][i] =  0.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + i) ) / nf;
    for (int i = 1; i < len; i++)
      for (int j = 1; j < i; j++){ /* limit is i, because we are dealing with symmetric matrix */
        Q[j][i] = Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)  *   *(*(ptrR + k) + k + j) ) / nf;
      } /* outer product of u1 and u1T is computed and stored in matrix Q*/
      /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/

    /* Q' = [I, zero; zero, Q]*/
    /* here I is m -len by m - len */
    /* after each refelection Q = Qold Q'*/
    /* Qold = [a, b; c, d]*/
    /* where a is m - len by m - len, d is len by len, c is len by m - len and b is m - len by len*/
    /* updating new  Q =  [a, b; c, d] [I,zero; zero, Q] = [a, bQ; c, dQ]*/
    /* bQ is m - len by m (i.e k by m), store bQ in QTR */
    /* dQ is len by len, store dQ in QBR */
    /* block multiplication */
    float QTR[k][len];
    float QBR[len][len];
    /* computing top right block */
    for (int i = 0; i < k; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += *(*(ptrQ + k + z ) + j )   *  Q[z][j];
        QTR[i][j] = init;
        }

      /* computing bottom right block */
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init +=  *(*(ptrQ + k + z ) + j + k )    * Q[z][j];
        QBR[i][j] = init;
      }

      /* copying the computed blocks to the memory of ptrQ */
      for (int i = 0; i < k; i++)
        for (int j = 0; j < len; j++)
          *(*(ptrQ + k + j ) + i) = QTR[i][j];
      for (int i = 0; i < len; i++)
        for (int j = 0; j < len; j++)
          *(*(ptrQ + k + j ) + i + k ) = QBR[i][j];
      /*Matrix Q is updated properly, it is double checked! */


      /* now updating R*/
      /* Easiest approach is ptrR = ptrQ times ptrA, but we can reduce the memory conception and computation time by block matrix way */
      /* So, new R = Q' times Rold, where Rold is m by n matrix*/
      /* Recall Q' = [I, zero; zero, Q] m by m matrix */
      /* R = Q' Rold = [I, zero; zero, Q] [a, b; c, d] = [a,b; Qc, Qd]*/
      /* a is m - len by m - len, c is len by m - len, b is n - m + len by m - len and d is len by n - m + len*/
      /* where Qc is len by m - len, store Qc in BLR*/
      /* Qd is len by n - m + len, store Qd in BRR*/
      float BLR[len][k];
      float BRR[len][n-k];
      /* calculating the BLR */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < k; j++){
          init = 0.0;
          for (int z = 0; z < len; z++)
            init += Q[i][z] * *(*(ptrR + j) + k + z);
          BLR[i][j] = init;
        }
      /* calculating BRR */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < (n-k) j++){
          init = 0.0;
          for (int z = 0; z < len; z++)
            init += Q[i][z] * *(*(ptrR + k + j) + k + z);
          BRR[i][j] = init;
        }
      /* copying the computed blocks to the memory of ptrR */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < k; j++)
          *(*(ptrR + j) + k + i) = BLR[i][j];
      for (int i = 0; i < len; i++)
        for (int j = 0; j < (n-k); j++)
          *(*(ptrR + k + j) + k + i) = BRR[i][j];
      /* matrix R is updated properly, it is double checked!*/
  }
  return 0;
}



/*------We have done the single reflection using house holder transformation  */

/*------------------------------------------------------------------------------------------- */


/*QR Decomposition */
/* of matrix A whose size is m cross n */
int QR_HHdRl_Db (double **ptrA, double **ptrQ, double **ptrR, int m, int n){
  /* ptrA is the pointer of pointers of n array, each array has length m */
  /* ptrQ is the pointer of pointers of m arrays, each array has length m */
  /* ptrR is the pointer of pointers of n arrays, each array has length m */
   /* this QR designed for specific dimension, although modifying them is not a big deal!*/
  int status;
  if (m < n)
    for (int i = 0; i < m - 1; i++ ){
      status = QQ_HH_Itr_Db(**ptrA, **ptrQ, **ptrR, m, n, (i + 1) );
      if (status == 3)
        printf("Non - Singular matrix \n");
      if (status == 2)
        printf("Dimension mismatch \n");
    }
  else
    for (int i = 0; i < n - 1; i++ ){
      status = QQ_HH_Itr_Db(**ptrA, **ptrQ, **ptrR, m, n, (i + 1) );
      if (status == 3)
        printf("Non - Singular matrix \n");
      if (status == 2)
        printf("Dimension mismatch \n");
    }
  return 0;
}

int QR_HHdRl_Fl (float **ptrA, float **ptrQ, float **ptrR, int m, int n){
  /* ptrA is the pointer of pointers of n array, each array has length m */
  /* ptrQ is the pointer of pointers of m arrays, each array has length m */
  /* ptrR is the pointer of pointers of n arrays, each array has length m */
  /* this QR designed for specific dimension, although modifying them is not a big deal!*/
  if (m < n)
    for (int i = 0; i < m - 1; i++ ){
      status = QQ_HH_Itr_Db(**ptrA, **ptrQ, **ptrR, m, n, (i + 1) );
      if (status == 3)
        printf("Non - Singular matrix \n");
      if (status == 2)
        printf("Dimension mismatch \n");
    }
  else
    for (int i = 0; i < n - 1; i++ ){
      status = QQ_HH_Itr_Db(**ptrA, **ptrQ, **ptrR, m, n, (i + 1) );
      if (status == 3)
        printf("Non - Singular matrix \n");
      if (status == 2)
        printf("Dimension mismatch \n");
    }
  return 0;
}



/* -------------------------------------------------------------------------*/
/* Least square solvers */


/* For the overdeterminated equations*/
/* Too many equations for unknown number of variables */
/* A m by n matrix m >= n   A x = b, b is  m by one and x is n by one*/
/* solution A = QR by householder transformation, Q is m by m, R is m by n*/
/* QRx = b => Rx = Q^{T} b */
/* then by backward substituion method find x vector*/

int Over_det_solver_Fl (float **ptrA, float **ptrQ, float **ptrR, int m, int n, float *ptrB, float *ptrX){
  /* ptrA, ptrB are the same size, ptrQ is m by m, ptrB is m by one, ptrX is n by one number of unknowns!*/
  if (m < n)
    return 4; /* not an overdeterminated case */
  QR_HHdRl_Fl(float **ptrA, float **ptrQ, float **ptrR, int m, int n); /* To find A = QR*/
  float B[n]; /* to store Q^{T} b */
  float z;
  for (int i = 0; i < n; i++){
    z = 0.0;
    for (int j = 0; j < m; j++)
      z += *(*(ptrQ + i)  + j)   *  *(ptrB + j);
    B[i] = z;
  }
  /* now need to do backward substitution*/
  if ( *(*(ptrR + n - 1) + n - 1) == 0.0 )
    return 3;
  *(ptrX + n - 1) = B[n-1] / *(*(ptrR + n - 1) + n - 1);
  for (int i = n - 2; i >= 0; i-- ){
    z = 0.0;
    for (int j = n - 1; j < i; j--)
      z +=  *(*(ptrR + j) + i)  *  B[j];
    if ( *(*(ptrR + i) + i) == 0.0 )
      return 3;
    *(ptrX + i) = ( B[i] - z ) / *(*(ptrR + i) + i);
  }
  return 0;
}


int Over_det_solver_Db (double **ptrA, double **ptrQ, double **ptrR, int m, int n, double *ptrB, double *ptrX){
  /* ptrA, ptrB are the same size, ptrQ is m by m, ptrB is m by one, ptrX is n by one number of unknowns!*/
  if (m < n)
    return 4; /* not an overdeterminated case */
  QR_HHdRl_Db(double **ptrA, double **ptrQ, double **ptrR, int m, int n); /* To find A = QR*/
  double B[n]; /* to store Q^{T} b */
  double z;
  for (int i = 0; i < n; i++){
    z = 0.0;
    for (int j = 0; j < m; j++)
      z += *(*(ptrQ + i)  + j)   *  *(ptrB + j);
    B[i] = z;
  }
  /* now need to do backward substitution*/
  if ( *(*(ptrR + n - 1) + n - 1) == 0.0 )
    return 3;
  *(ptrX + n - 1) = B[n-1] / *(*(ptrR + n - 1) + n - 1);
  for (int i = n - 2; i >= 0; i-- ){
    z = 0.0;
    for (int j = n - 1; j < i; j--)
      z +=  *(*(ptrR + j) + i)  *  B[j];
    if ( *(*(ptrR + i) + i) == 0.0 )
      return 3;
    *(ptrX + i) = ( B[i] - z ) / *(*(ptrR + i) + i);
  }
  return 0;
}




/* For the undetermined equations */
/* A n by m matrix n <= m, A x = b n by one, x is m by one */
/* A^{T} now m by n matrix*/ 
/* Compute A^{T} = QR by householder transformation, Q is m by m, R is m by n */
/* A^{T} x = R^{T} Q^{T} x = [R1^{T} , zero] [Q1^{T}; Q2^{T}] x = b*/
/* R1^{T} Q1^{T} x = b*/
/* let Q1^{T}x = y*/
/* use backward substituion method to compute y*/
/* now use the orthogonal nature of the Q1 i.e Q^{Inv} = Q^{T}*/
/* x = Q1 y*/


/* Assume that given matrix is came with transpose itself, our ideaology, we never create any memory in the calling functions, so 
 * the main function has to do how to create and destroy memory for that reason we never create or destroy memory in the functions!*/

int Under_det_solver_Fl (float **ptrAT, float **ptrQ, float **ptrR, int m, int n, float *ptrB, float *ptrX){
  /* ptrAT, ptrR have same size m by n, ptrQ is m by m, ptrB is n by one, ptrX is m by one*/
  if (m < n)
    return 4; /* not an underdeterminated case */
  QR_HHdRl_Fl(float **ptrAT, float **ptrQ, float **ptrR, int m, int n); /* To find AT = QR*/
  float Y[n]; /* to compute Q1y*/
  float z;
  /* now need to do backward substitution*/
  if ( *(*(ptrR + n - 1) + n - 1) == 0.0 )
    return 3;
  Y[n - 1] = *(ptrB + n - 1) / *(*(ptrR + n - 1) + n - 1);
  for (int i = n - 2; i >= 0; i-- ){
    z = 0.0;
    for (int j = n - 1; j < i; j--)
      z +=  *(*(ptrR + i) + j)  *  Y[j];
    if ( *(*(ptrR + i) + i) == 0.0 )
      return 3;
    Y[i] = ( *(ptrB + i) - z ) / *(*(ptrR + i) + i);
  }
  /*computing X = Q times Y */
  for (int i = 0; i < m; i++){
    z = 0.0;
    for (int j = 0; j < n; j++){
      z += Y[j] * *(*(ptrQ + j) + i);
    }
    *(ptrX + i) = z;
  }
  return 0;
}


int Under_det_solver_Db (double **ptrAT, double **ptrQ, double **ptrR, int m, int n, double *ptrB, double *ptrX){
  /* ptrAT, ptrR have same size m by n, ptrQ is m by m, ptrB is n by one, ptrX is m by one*/
  if (m < n)
    return 4; /* not an underdeterminated case */
  QR_HHdRl_Db(double **ptrAT, double **ptrQ, double **ptrR, int m, int n); /* To find A^{T} = QR*/
  double Y[n]; /* to compute Q1y*/
  double z;
  /* now need to do backward substitution*/
  if ( *(*(ptrR + n - 1) + n - 1) == 0.0 )
    return 3;
  Y[n - 1] = *(ptrB + n - 1) / *(*(ptrR + n - 1) + n - 1);
  for (int i = n - 2; i >= 0; i-- ){
    z = 0.0;
    for (int j = n - 1; j < i; j--)
      z +=  *(*(ptrR + i) + j)  *  Y[j];
    if ( *(*(ptrR + i) + i) == 0.0 )
      return 3;
    Y[i] = ( *(ptrB + i) - z ) / *(*(ptrR + i) + i);
  }
  /*computing X = Q times Y */
  for (int i = 0; i < m; i++){
    z = 0.0;
    for (int j = 0; j < n; j++){
      z += Y[j] * *(*(ptrQ + j) + i);
    }
    *(ptrX + i) = z;
  }
  return 0;
}

/*Solving the linear least squares was done! */
/* review latter to check the correctness! */
/*-----------------------------------------------------------------------------------------------------------*/







