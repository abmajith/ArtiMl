

/*  *ptrA1 *ptrA2 ........*ptrAm */
/*  a11   a12   a13 .........a1m  
*   a21   a22   a23 .........a2m
     .............
*   an1   an2   an3 .........anm */
/* m column and n row, matrix A is n by m*/

/*b11 b12 b13 ........b1p 
 *b21 b22 b23 ........b2p
  ..................
  bm1 bm2 ............bmp */

/* c11 c12 c13 .... c1p
 * c21 c22 c23 .... c2p
 * ...............
 * cn1 cn2 cn3 .... cnp */

/* This setting is consistent thoughout matrix operations */


/* For this type of matrix storage, matrix multiplication has to do by the following way */
/* for i 0ton  */
/*  for j 0top  */
/*   cij=zero*/
/*   for k0tom */
/*     cij +=Aik.Bkj*/




/*Matrix multiplication C = A * b where  */
/*A has n by m and B is m by p, so C should hold n by p */
void Mat_Mul_Db(double  **ptrA, double **ptrB, double **ptrC, int n, int m, int p){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n*/
  /* ptrB is the pointer to pointers of p arrays, each arrays has size m*/
  /* ptrC is the pointer to pointers of p arrays, each arrays has size n*/
  /*cij = sumation over all k  aik * bkj */
  double init;
  int i,j,k;
  for (i = 0; i < n; i++){   /*you can reduce the matrix computation time by doing block matrix compuation*/
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + k) + i )   *  *(*(ptrB + j) + k);
      *(*(ptrC + j) + i) = init;
    }
  }
  return; /* this function stores the matrix multiplication value in the pointer of pointer double type ptrC given as input */
} /* This code was double checked, its correct*/

void Mat_Mul_Fl(float  **ptrA, float **ptrB, float **ptrC, int n, int m, int p){
  float init;
  int i,j,k;
  for (i = 0; i < n; i++){   /*you can reduce the matrix computation time by doing block matrix compuation*/
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + k) + i )   *  *(*(ptrB + j) + k);
      *(*(ptrC + j) + i) = init;
    }
  }
  return;
} /* This code is same as before, just float type instead of double type*/

void Mat_Mul_Int(int  **ptrA, int **ptrB, int **ptrC, int n, int m, int p){
  int init;
  int i,j,k;
  for (i = 0; i < n; i++){
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + k) + i )   *  *(*(ptrB + j) + k);
      *(*(ptrC + j) + i) = init;
    }
  }
  return;
}/*This code is same as before, just integer type instead of float type */

void Mat_Mul_Ln(long  **ptrA, long **ptrB, long **ptrC, int n, int m, int p){
  long init;
  int i,j,k;
  for (i = 0; i < n; i++){
    for (j = 0; j < p; j++){
      init = 0;
      for (k = 0; k < m; k++)
        init += *(*(ptrA + k) + i )   *  *(*(ptrB + j) + k);
      *(*(ptrC + j) + i) = init;
    }
  }
  return;
}/*This code is same as before, just long type instead of integer type */



/*------matrix multiplication function code is double checked -----------------*/

/*---------------------------------------------------------------------------------------*/



/* Diagnonal multiplication*/
/* multiplication of diagonal elements of the given m by m matrix*/
int Mat_Daig_Mul_Int(int **ptrA, int m){
  /* ptrA is the pointer to pointers of m arrays, each arrays has the size m */
  int d = 1;
  int i;
  for (i = 0; i < m; i++)
    d  *=  *(*(ptrA + i) + i);
  return d; /* return the integer type of multiplication of all diagnonal element */
} /*This code was double checked */

float Mat_Daig_Mul_Fl(float **ptrA, int m){
  float d = 1.0;
  int i;
  for (i = 0; i < m; i++)
    d *= *(*(ptrA + i) + i);
  return d;
}/*This code is same as previous one, but float type*/

double Mat_Daig_Mul_Db(double **ptrA, int m){
  double d = 1.0;
  int i;
  for (i = 0; i < m; i++)
    d *= *(*(ptrA + i) + i);
  return d;
} /* This code is same as previous one, but double type */


long Mat_Daig_Mul_Ln(long **ptrA, int m){
  /* ptrA is the pointer to pointers of m arrays, each arrays has the size m */
  long d = 1;
  int i;
  for (i = 0; i < m; i++)
    d  *=  *(*(ptrA + i) + i);
  return d; /* return the long type of multiplication of all diagnonal element */
} /*This code is same as previous one, but long type */


/*-----matrix diagonal multiplication function code is double checked--------*/

/*--------------------------------------------------------------------------------------- */


/* Matrix Trasnpose */
/* given a matrix A of size m by n transpose A^T*/
void Mat_Trans_Db(double **ptrA, double **ptrB, int m, int n){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n */
  /* ptrB is the pointer to pointers of n arrays, each arrays has size m */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      *(*(ptrB + j) + i) =  *(*(ptrA + i) + j);
  return;/* It just transpose and store in the given pointers of array*/
} /* This code was double checked*/


void Mat_Trans_Fl(float **ptrA, float **ptrB, int m, int n){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n */
  /* ptrB is the pointer to pointers of n arrays, each arrays has size m */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      *(*(ptrB + j) + i) =  *(*(ptrA + i) + j);
  return;
} /* Same as previous one, but float type*/

void Mat_Trans_Int(int **ptrA, int **ptrB, int m, int n){
  /* ptrA is the pointer to pointers of m arrays, each arrays has size n */
  /* ptrB is the pointer to pointers of n arrays, each arrays has size m */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      *(*(ptrB + j) + i) =  *(*(ptrA + i) + j);
  return;
}/*same as previous one, but integer type*/

void Mat_Trans_Ln(long **ptrA, long **ptrB, int m, int n){
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

