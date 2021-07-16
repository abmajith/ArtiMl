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


/* For this type of matrix storage, matrix multiplication has to do by the following way */
/* for i0ton  */
/* for j0top  */
/*cij=zero*/
/* for k0tom */
/* cij +=Aik.Bkj*/




/*Matrix multiplication C = A * b where  */
/*A has n cross m and B is m cross p, so C should hold n cross p */
void Mat_Multi_Db(double  **ptrA, double **ptrB, double **ptrC, int n, int m, int p){
  /* ptrA is the pointer to pointers of n arrays, each arrays has the size of m*/
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



/*------matrix multiplication function code is double checked -----------------*/

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


/*-----matrix diagonal multiplication function code is double checked--------*/

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

/* Note: for real programs, do check that the input and output can hold the same data type otherwise we need to cast the output to higher data type */

/*-----square of norm two of given array function is double checked, remember the note:----------- */

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


/*--------Transpose of given matrix function is doble checked------------*/
/* Caution: usage of these type in the main program can create the problem of hiding memory blocks*/
/*order of memory allocation*/

/*--------------------------------------------------------------------------------------- */




/* here after the pointer to pointer of a matrix storage method*/
/* n by m matrix*/
/*ptrA1 ptrA1 ptrA2 ..... prAm */
/*a11    a12   a13         a1m */
/*a21    a22   a23         a2m  */
/*....               .......... */
/*an1    an2   an3         anm */



/*House Holder Reflection*/
/* Iteration of House holder reflection on the kth column*/
void HHdRl_Itr_Db(double **ptrR, int m, int k1 ){ /* R is the matrix of m cross m, and relect on the ptrRk1 th array*/
  /*k1 should be non zero positive value */
  /*ptrR is the pointer to pointers of m arrays, each array has size m */
  if ( k1 >= m || k1 < 1 )
    return; /* if k1 is equal or more than m nothing to do*/
  int k = k1 - 1; /* cprogram convention for the index*/
  int len = m - k; /* dimension of the submatrix and this is correct checked multiple times*/
  double init = 0.0;
  double pivot = *(*(ptrR + k) + k); /*Rkk */

  double *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */
  ptr = *(ptrR + k) + k; /* address of Rk,k */

  double sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (double) Sq_Normtwo_Db(ptr, len);
  alpha = (double) sqrt(sigma);
  /* if calculated norm is same as the pivot then divide by zero occurs in the rotation so */
    if (alpha == pivot)
      return;

  double nf; /* normalization factor to compute Q for the len by len submatrix*/
  nf = sigma - alpha * pivot;
  double Q[len][len], BR[len][len]; /* to store the rotated submatrix of R for making the k+1K+2K..m to zero */
  /*Q = Identity matrix len by  len - (1over nf ) outer  product of u1 and u1T */
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
  return;
}/* This code is double checked!*/



/*House Holder Reflection*/
/* Iteration of House holder reflection on the kth column*/
void HHdRl_Itr_Fl(float **ptrR, int m, int k1 ){ /* R is the matrix of m cross m, and relect on the ptrRk1 th array*/
  /*k1 should be non zero positive value */
  /*ptrR is the pointer to pointers of m arrays, each array has size m */
  if ( k1 >= m || k1 < 1 )
    return; /* if k1 is equal or more than m nothing to do*/
  int k = k1 - 1; /* cprogram convention for the index*/
  int len = m - k; /* dimension of the submatrix and this is correct checked multiple times*/
  float init = 0.0;
  float pivot = *(*(ptrR + k) + k); /*Rkk */

  float *ptr; /* to work with column [Rk,k Rk+1k Rk+2k ...Rmk]T */
  ptr = *(ptrR + k) + k; /* address of Rk,k */

  float sigma, alpha; /* sigma for square of norm two, alpha for norm two*/
  sigma = (float) Sq_Normtwo_Db(ptr, len);
  alpha = (float) sqrt(sigma);
  /* if calculated norm is same as the pivot then divide by zero occurs in the rotation so */
    if (alpha == pivot)
      return;

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
  return;
}/* This code is double checked!*/
/* the above code is same as before, but for float type!*/


/* House holder reflection on the given square matrix iteration was done */

/*------------------------------------------------------------------------------------- */




/*Given matrix A of size m by n, find Q and R decomposition of A, m <= n i.e rank of A is atmost m*/
/* size of Q is m by m*/


/* need to review */
void QR_HH_Itr_Db (double **ptrA, double **ptrQ, double **ptrR, int m, int n, int k1){
  /* ptrA is the pointer to pointer of n arrays, each arrays has size m */
  /* ptrQ is the pointer to pointer of m arrays, each arrays has size m */
  /* ptrR is the pointer to pointer of n arrays, each arrays has size m*/

  if ( k1 >= m || k1 < 0 || m >= n) /*when k1 is equal or more than m no point in doing the reflection */
    return;
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
   /* if calculated norm is same as the pivot then divide by zero occurs in the reflection */
    if (alpha == pivot)
      return;
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
   /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      return;
    nf = sigma - alpha * pivot;
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
  return;
}

void QR_HH_Itr_Fl (float **ptrA, float **ptrQ, float **ptrR, int m, int n, int k1){
  /* ptrA is the pointer to pointer of n arrays, each arrays has size m */
  /* ptrQ is the pointer to pointer of m arrays, each arrays has size m */
  /* ptrR is the pointer to pointer of n arrays, each arrays has size m*/

  if ( k1 >= m || k1 < 0 || m >= n) /*when k1 is equal or more than m no point in doing the reflection */
    return;
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
   /* if calculated norm is same as the pivot then divide by zero occurs in the reflection */
    if (alpha == pivot)
      return;
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
   /* if calculated norm is same as the pivot then divide by zero occurs in the rotation */
    if (alpha == pivot)
      return;
    nf = sigma - alpha * pivot;
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
  if (m > n) /* this QR designed for specific dimension, although modifying them is not a big deal!*/
    return;
  for (int i = 0; i < m - 1; i++){
    QQ_HH_Itr_Db(**ptrA, **ptrQ, **ptrR, m, n, (i + 1) );
  }
  return;
}

void QR_HHdRl_Fl (float **ptrA, float **ptrQ, float **ptrR, int m, int n){
  /* ptrA is the pointer of pointers of n array, each array has length m */
  /* ptrQ is the pointer of pointers of m arrays, each array has length m */
  /* ptrR is the pointer of pointers of n arrays, each array has length m */
  if (m > n) /* this QR designed for specific dimension, although modifying them is not a big deal!*/
    return;
  for (int i = 0; i < m - 1; i++){
    QQ_HH_Itr_Fl(**ptrA, **ptrQ, **ptrR, m, n, (i + 1) );
  }
  return;
}

/* each function instead of returning void, putting some number can transfer error type (from some home made dictonary). Think about it*/
