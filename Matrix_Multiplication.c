#include <math.h> /*i am using sqrt function, so its needed */


/* *ptrA0   a11 a12 a13 .........a1m  
* *ptrA1    a21 a22 a23 .........a2m
  .............
* *ptrAend.  an1 an2 an3 .........anm */

/*b11 b12 b13 ........b1p 
 *b21 b22 b23 ........b2p
  ..................
  bm1 bm2 ............bmp */

/* c11 c12 c13 .... c1p
 * c21 c22 c23 .... c2p
 * ...............
 * cn1 cn2 cn3 .... cnp */

/* This setting (pointer of pointer to a row array) holds only for the matrix multiplication until matrix transpose the onwards 
 * each pointer from the pointer of arrays refers the column array */




/*Matrix multiplication C = A * b where  */
/*A has n cross m and B is m cross p, so C should hold n cross p */
void Mat_Multi_Db(double  **ptrA, double **ptrB, double **ptrC, int n, int m, int p){
  /* ptrA is the pointer to pointers of n arrays (row matrix), each arrays has the size of m*/
  /* ptrB is the pointer to pointers of m arrays, each arrays has the size of p*/
  /* ptrB is the pointer to pointers of n arrays, each arrays has the size of p*/
  /*cij = sumation over all k  aik * bkj */
  double init;
  int i,j,k;
  for (i = 0; i < n; i++){ /*you can reduce the matrix computation time by doing block matrix compuation*/
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + i) + k )   *  *(*(ptrB + k) + j);
      *(*(ptrC + i) + j) = init;
    }
  }
  return; /* this function stores the matrix multiplication value in the pointer of pointer double type ptrC given as input */
} /* This code was double checked, its correct*/

void Mat_Multi_Fl(float  **ptrA, float **ptrB, float **ptrC, int n, int m, int p){
  float init;
  int i,j,k;
  for (i = 0; i < n; i++){
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + i) + k )   *  *(*(ptrB + k) + j);
      *(*(ptrC + i) + j) = init;
    }
  }
  return;
} /* This code is same as before, just float type instead of double type*/

void Mat_Multi_Int(int  **ptrA, int **ptrB, int **ptrC, int n, int m, int p){
  int init;
  int i,j,k;
  for (i = 0; i < n; i++){
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + i) + k )   *  *(*(ptrB + k) + j);
      *(*(ptrC + i) + j) = init;
    }
  }
  return;
}/*This code is same as before, just integer type instead of float type */

void Mat_Multi_Ln(long  **ptrA, long **ptrB, long **ptrC, int n, int m, int p){
  long init;
  int i,j,k;
  for (i = 0; i < n; i++){
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + i) + k )   *  *(*(ptrB + k) + j);
      *(*(ptrC + i) + j) = init;
    }
  }
  return;
}/*This code is same as before, just long type instead of integer type */



/*------matrix multiplication function code is ok -----------------*/

/*---------------------------------------------------------------------------------------*/





/* Diagnonal multiplication*/
/* multiplication of diagonal elements for the given m cross m matrix*/
int Mat_Daig_Multi_Int(int **ptrA, int m){
  /* ptrA is the pointer to pointers of m arrays, each arrays has the size m */
  int d = 1;
  int i;
  for (i = 0; i < m; i++)
    d  *=  *(*(ptrA + i) + i);
  return d; /* return the integer type of multiplication of all diagnonal element */
} /*This code was double checked */

float Mat_Daig_Multi_Fl(float **ptrA, int m){
  float d = 1.0;
  int i;
  for (i = 0; i < m; i++)
    d *= *(*(ptrA + i) + i);
  return d;
}/*This code is same as previous one, but float type*/

double Mat_Daig_Multi_Db(double **ptrA, int m){
  double d = 1.0;
  int i;
  for (i = 0; i < m; i++)
    d *= *(*(ptrA + i) + i);
  return d;
} /* This code is same as previous one, but double type */


long Mat_Daig_Multi_Ln(long **ptrA, int m){
  /* ptrA is the pointer to pointers of m arrays, each arrays has the size m */
  long d = 1;
  int i;
  for (i = 0; i < m; i++)
    d  *=  *(*(ptrA + i) + i);
  return d; /* return the long type of multiplication of all diagnonal element */
} /*This code is same as previous one, but long type */


/*-----matrix diagonal multiplication function code is ok--------*/

/*--------------------------------------------------------------------------------------- */



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

/* Note: for real programs, do check that the input and output can hold the same data type otherwise we need to cast and lift the data
 * type according to the expected output*/

/*-----square of norm two of given array function is ok----------- */

/*-------------------------------------------------------------------------------------------------------- */



/* Matrix Trasnpose */
/* given a matrix A of size m cross n transpose A^T*/
void Mat_Tran_Db(double **ptrA, double **ptrB, int m, int n){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n */
  /* ptrB is the pointer to pointers of n arrays, each arrays has size m */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      *(*(ptrB + j) + i) =  *(*(ptrA + i) + j);
  return;/* It just transpose and store in the given pointers of array*/
} /* This code was double checked*/


void Mat_Tran_Fl(float **ptrA, float **ptrB, int m, int n){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n */
  /* ptrB is the pointer to pointers of n arrays, each arrays has size m */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      *(*(ptrB + j) + i) =  *(*(ptrA + i) + j);
  return;
} /* Same as previous one, but float type*/

void Mat_Tran_Int(int **ptrA, int **ptrB, int m, int n){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n */
  /* ptrB is the pointer to pointers of n arrays, each arrays has size m */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      *(*(ptrB + j) + i) =  *(*(ptrA + i) + j);
  return;
}/*same as previous one, but integer type*/

void Mat_Tran_Ln(long **ptrA, long **ptrB, int m, int n){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n */
  /* ptrB is the pointer to pointers of n arrays, each arrays has size m */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      *(*(ptrB + j) + i) =  *(*(ptrA + i) + j);
  return;
}/* Same as previous one, but long type*/


/*--------Transpose of given matrix function is ok------------*/

/*--------------------------------------------------------------------------------------- */



/* n by m matrix*/
/*ptrA0 ptrA1 ptrA2 ..... prAend */
/*a11    a12   a13         a1m */
/*a21    a22   a23         a2m  */
/*....               .......... */
/*an1    an2   an3         anm */



/*House Holder Reflection*/
/* Iteration of House holder reflection on the kth column*/
void HHdRl_Itr_Db(double **ptrR, int m, int k1 ){ /* R is the matrix of m cross m, and relect on the k1 th column*/
  /*k1 should be non zero positive value */
  /*ptrA is the pointer to pointers of m arrays (column matrix!), each array has size of m */
  if ( k1 >= m)
    return;
  int k = k1 - 1;
  int len = m - k; /* dimension of the submatrix*/
  double init = 0.0;
  double pivot = *(*(ptrR + k) + k); /*Rkk */

  double *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */
  ptr = *(ptrR + k) + k; /* address of Rk,k */

  double sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (double) Sq_Normtwo_Db(ptr, len);
  alpha = (double) sqrt(sigma);
  /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      return;

  double nf; /* normalization factor to compute Q for the len cross len submatrix*/
  nf = sigma - alpha * pivot;
  double Q[len][len], BR[len][len]; /* to store the submatrix rotation for the submatrix of R to make the k+1K+2K..m to zero */
  /*Q = Ilen cross len - (1over nf ) outer  product of u1T and u1 */
  for (int i = 0; i < len; i++)
    for (int j = 0; j <= i; j++){ /* limit is i, because we are dealing with symmetric matrix */
      if ((i != 0) && (j != 0))
        if (i != j){
          Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)  *   *(*(ptrR + k) + k + j) ) / nf;
          Q[j][i] = Q[i][j]; /* can save space by avoiding the lower corner let see in future*/
        }
        else {
          Q[i][j] = 1.0 - (  *(*(ptrR + k) + k + i)   *  *(*(ptrR + k) + k + j) ) / nf;
        }
      else
      { 
        if  ( i != 0 &&  j == 0)
         if ( i != j ){
          Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)   * (  pivot - alpha ) ) / nf;
          Q[j][i] = Q[i][j];
         }
         else {
            Q[i][j] = 1.0 - (  *(*(ptrR + k) + k + i)   * (  pivot - alpha ) ) / nf;
         }
        else 
          if ( i == 0 &&  j != 0){
            if ( i != j) {
              Q[i][j] = 0.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + j) ) / nf;
              Q[j][i] = Q[i][j];
            }
            else
               Q[i][j] = 1.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + j) ) / nf;
          }
          else
          {
            Q[i][j] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
          }
      }
    }
  /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/

  for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j + k) + z + k);
        BR[i][j] = init;
      }
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++)
        *(*(ptrR + j + k) + i + k) = BR[i][j];

  if ( k != 0)
  {
    double BL[len][m - len];
    for (int i = 0; i < len ; i++)
      for (int j = 0; j < (m - len); j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j) + z + k);
        BL[i][j] = init;
      }
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m - len; j++)
        *(*(ptrR + j) + i + k) = BL[i][j];
  }
  return;
}/* This code is ok!*/


void HHdRl_Itr_Fl(float **ptrR, int m, int k1 ){ /* R is the matrix of m cross m, and relect on the k1 th column*/
  /*k1 should be non zero positive value */
  /*ptrA is the pointer to pointers of m arrays (column matrix!), each array has size of m */
  if ( k1 >= m)
    return;
  int k = k1 - 1;
  int len = m - k; /* dimension of the submatrix*/
  float init = 0.0;
  float pivot = *(*(ptrR + k) + k); /*Rkk */

  float *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */
  ptr = *(ptrR + k) + k; /* address of Rk,k */

  float sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (float) Sq_Normtwo_Fl(ptr, len);
  alpha = (float) sqrt(sigma);
  /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      return;

  float nf; /* normalization factor to compute Q for the len cross len submatrix*/
  nf = sigma - alpha * pivot;
  float Q[len][len], BR[len][len]; /* to store the submatrix rotation for the submatrix of R to make the k+1K+2K..m to zero */
  /*Q = Ilen cross len - (1over nf ) outer  product of u1T and u1 */
  for (int i = 0; i < len; i++)
    for (int j = 0; j <= i; j++){ /* limit is i, because we are dealing with symmetric matrix */
      if ((i != 0) && (j != 0))
        if (i != j){
          Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)  *   *(*(ptrR + k) + k + j) ) / nf;
          Q[j][i] = Q[i][j]; /* can save space by avoiding the lower corner let see in future*/
        }
        else {
          Q[i][j] = 1.0 - (  *(*(ptrR + k) + k + i)   *  *(*(ptrR + k) + k + j) ) / nf;
        }
      else
      {
        if  ( i != 0 &&  j == 0)
         if ( i != j ){
          Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)   * (  pivot - alpha ) ) / nf;
          Q[j][i] = Q[i][j];
         }
         else {
            Q[i][j] = 1.0 - (  *(*(ptrR + k) + k + i)   * (  pivot - alpha ) ) / nf;
         }
        else
          if ( i == 0 &&  j != 0){
            if ( i != j) {
              Q[i][j] = 0.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + j) ) / nf;
              Q[j][i] = Q[i][j];
            }
            else
               Q[i][j] = 1.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + j) ) / nf;
          }
          else
          {
            Q[i][j] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
          }
      }
    }
  /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/

  for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j + k) + z + k);
        BR[i][j] = init;
      }
    for (int i = 0; i < len; i++)
      for (int j = 0; j < len; j++)
        *(*(ptrR + j + k) + i + k) = BR[i][j];

  if ( k != 0)
  {
    float BL[len][m - len];
    for (int i = 0; i < len ; i++)
      for (int j = 0; j < (m - len); j++){
        init = 0.0;
        for (int z = 0; z < len; z++)
          init += Q[i][z] * *(*(ptrR + j) + z + k);
        BL[i][j] = init;
      }
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m - len; j++)
        *(*(ptrR + j) + i + k) = BL[i][j];
  }
  return;
}/* This code is ok!, same as previous one but float type*/

/* need to review */
void QR_HH_Itr_Db (double **ptrA, double **ptrQ, double **ptrR, int m, int n, int k1){
  if (k1 >= m)
    return;
  int k = k1 - 1;
  int len = m - k; /* dimension of the submatrix*/
  double init = 0.0;
  double pivot; /*Rkk */

  double *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */

  double sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  /* if calculated norm is same as the modulus of pivot then divide by zero occurs in the rotation */
  double nf; /* normalization factor */
  if (k == 0){
    pivot = **ptrA;
    ptr = *ptrA;
    sigma = (double) Sq_Normtwo_Db(ptr, len);
    alpha = (double) sqrt(sigma);
   /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      return;
    nf = sigma - alpha * pivot;
    /* computing the Q matrix first Im*m - v outer product v / nf*/
    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i ; j++){
         if ((i != 0) && (j != 0))
            if (i != j){
                *(*(ptrQ + j) + i)  =  *(*(ptrQ + i) + j) = 0.0 - (  *(*ptrA + i)  *   *(*ptrA + j) ) / nf;
              }
             else {
                 *(*(ptrQ + i) + j) = 1.0 - (  *(*ptrR + i)   *  *(*ptrR + j) ) / nf;
              }
         else
          {
            if  ( i != 0 &&  j == 0)
              if ( i != j ){
                *(*(ptrQ + j) + i)  =  *(*(ptrQ + i) + j) = 0.0 - (  *(*ptrR + i)   * (  pivot - alpha ) ) / nf;
               }
              else {
                  *(*(ptrQ + i) + j)   = 1.0 - (  *(*ptrR + i)   * (  pivot - alpha ) ) / nf;
              }
            else
              if ( i == 0 &&  j != 0){
                if ( i != j) {
                    *(*(ptrQ + j) + i)  =  *(*(ptrQ + i) + j)  = 0.0 - (  ( pivot - alpha )  *  *(*ptrR + j) ) / nf;
                  }
                else
                  *(*(ptrQ + j) + i)  = 1.0 - (  ( pivot - alpha )  *  *(*ptrR + j) ) / nf;
              }
              else
              {
                *(*(ptrQ + j) + i)  = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
              }
          }
       }
    /* computing R matrix = QA*/
    /* Q is m by m and A is m by n*/
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++){
        init = 0.0;
        for (int z = 0; z < m; z++)
          init += *(*(ptrQ + z) + i)  *  *(*(ptrA + j) + k);
        *(*(ptrR + j) + i) = init;
      }
  }
  if ( k != 0){
    pivot = *(*(ptrR + k ) + k );
    ptr = *(ptrR + k) + k;
    sigma = (double) Sq_Normtwo_Db(ptr, len);
    alpha = (double) sqrt(sigma);
   /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      return;
    nf = sigma - alpha * pivot;
    /* computing the Q matrix first Ilen*len - v outer product v / nf*/
    double Q[len][len];
      for (int i = 0; i < len; i++)
        for (int j = 0; j <= i; j++){ /* limit is i, because we are dealing with symmetric matrix */
          if ((i != 0) && (j != 0))
            if (i != j){
              Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)  *   *(*(ptrR + k) + k + j) ) / nf;
              Q[j][i] = Q[i][j]; /* can save space by avoiding the lower corner let see in future*/
            }
            else {
              Q[i][j] = 1.0 - (  *(*(ptrR + k) + k + i)   *  *(*(ptrR + k) + k + j) ) / nf;
            }
          else
          {
            if  ( i != 0 &&  j == 0)
            if ( i != j ){
              Q[i][j] = 0.0 - (  *(*(ptrR + k) + k + i)   * (  pivot - alpha ) ) / nf;
              Q[j][i] = Q[i][j];
            }
            else {
                Q[i][j] = 1.0 - (  *(*(ptrR + k) + k + i)   * (  pivot - alpha ) ) / nf;
            }
            else
              if ( i == 0 &&  j != 0){
                if ( i != j) {
                  Q[i][j] = 0.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + j) ) / nf;
                  Q[j][i] = Q[i][j];
                }
                else
                  Q[i][j] = 1.0 - (  ( pivot - alpha )  *  *(*(ptrR + k) + k + j) ) / nf;
              }
              else
              {
                Q[i][j] = 1.0 - (  ( pivot - alpha )  *  ( pivot - alpha ) ) / nf;
              }
          }
        }
  /* Note Q is the symmetric matrix, you can reduce the computation space: try later*/
      /* updating the Q = [[a] [b][Q]; [c] [d][Q]] */
     /* block multiplication */
      double QTR[k][len];
      double QBR[len][len];
      /* calculating top right block */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < k; j++){
          init = 0.0;
          for (int z = 0; z < len; z++)
            init += *(*(ptrQ + k + z ) + j )   *  Q[k][j];
          QTR[i][j] = init;
        }
      /* calculating bottom right block */
      for (int i = 0; i < len; i++)
        for (int j = 0; j < len; j++){
          init = 0.0;
          for (int z = 0; z < len; z++)
            init +=  *(*(ptrQ + k + z ) + j + k )    * Q[k][j];
          QBR[i][j] = init;
        }
      /* copying the computed blocks to the address of given pointers*/
      for (int i = 0; i < k; i++)
        for (int j = 0; j < len; j++)
          *(*(ptrQ + k + j ) + i) = QTR[i][j];
      for (int i = 0; i < len; i++)
        for (int j = 0; j < len; j++)
          *(*(ptrQ + k + j ) + i + k ) = QBR[i][j];
      double BLR[len][k];
      double BRR[][];                                                           /*to be contintued here......................................... */
      /* calculating the */
      /* updating the R*/
  }
  return;
}


/*------We have done the single reflection using house holder transformation  */

/*------------------------------------------------------------------------------------------- */


/*QR Decomposition */
/* of matrix A whose size is m cross n */
void QR_HHdRl_Db (double **ptrA, double **ptrQ, double **ptrR, int m, int n){
  /* ptrA is the pointer of pointers of n array, each array has length m */
  /* ptrQ is the pointer of pointers of m arrays, each array has length m */
  /* ptrR is the pointer of pointers of n arrays, each array has length m */
  for (int i = 0; i < m - 1; i++){
    QQ_HH_Itr_Db (**ptrA, **ptrQ, **ptrR, m, n, (i + 1) );
  }
  return;
}

/* each function instead of returning void, putting some number can refer errors from some created list. Think about it*/
