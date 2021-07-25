#include <math.h> /*i am using sqrt function, so its needed */

/* Given non singular upper trianglular square matrix m by m, find Inverse */
/* By Backward substitution method */
/* note this one can lead severe numerical approximation inaccuracy, be caution with it!
 * its not because the code is badly written, its the draw back of backward substituion method which involves 
 * sequences of many division operation*/

/*[ a11 a12 a13]     [ ia11 ia12 ia13] = [1 0 0]     
 *[  0  a22 a23]     [  0   ia22 ia23] = [0 1 0] 
  [  0   0  a33]     [  0    0   ia33] = [0 0 1] */



int NSing_UpTriag_Mat_Inv_Backsub_Mthd_Db(double **ptrR, double **ptrIR, int m){
  double det = 1.0;
  for (int i = 0; i < m; i++) /* determinant of the upper triangular matrix is computed */
    det *= *(*(ptrR + i) + i);
  if (det == 0.0)
    return 1; /* indicating that the given matrix is not a upper triangular non-singular matrix*/
  for (int i = 0; i < m; i++) /* inverse of diagonal element of R matrix */
    *(*(ptrIR + i) + i) = 1.0 / *(*(ptrR + i) + i);
  double init;
  /* backward substitution method */
  for (int i = m - 2; i >= 0; i--) /* computing the ith row value */
    for (int j = m - 1; j > 0; j--){ /* computing inverse of ith row and jth column value */
      init = 0.0;
      for (int z = i + 1; z < m; z++)
        init += *(*(ptrR + z ) + i ) * *(*(ptrIR + j) + z);
      *(*(ptrIR + j) + i) = 0.0 - (double) (init / *(*(ptrR + i) + i);
    } /* This block is double checked */
  return 0;
}

int NSing_UpTriag_Mat_Inv_Backsub_Mthd_Fl(float **ptrR, float **ptrIR, int m){
  float det = 1.0;
  for (int i = 0; i < m; i++) /* determinant of the upper triangular matrix is computed */
    det *= *(*(ptrR + i) + i);
  if (det == 0.0)
    return 1; /* indicating that the given matrix is not a upper triangular non-singular matrix*/
  for (int i = 0; i < m; i++) /* diagonal of inverse of R matrix is computed */
    *(*(ptrIR + i) + i) = 1.0 / *(*(ptrR + i) + i);
  float init;
  /* backward substitution method */
  for (int i = m - 2; i >= 0; i--) /* computing the ith row value */
    for (int j = m - 1; j > 0; j--){ /* computing inverse of ith row and jth column value */
      init = 0.0;
      for (int z = i + 1; z < m; z++)
        init += *(*(ptrR + z ) + i )  *  *(*(ptrIR + j) + z);
      *(*(ptrIR + j) + i) = 0.0 - (float) (init / *(*(ptrR + i) + i);
    } /* This block is double checked */
  return 0;
}








/* need to document the algorithm properly, later only check, otherwise you going to waste too much time in this checking*/


/* hear we will write more extensive, stable adjoint method rather than backward substitution method to get better numerical stability!*/
void NSing_UpTriag_Mat_Inv_Adj_Abdul_Mthd_Db(double **ptrR, double **ptrIR, int m){
  double det = 1.0;
  for (int i = 0; i < m; i++)
    det *= *(*(ptrR + i) + i);
  for (int i = 0; i < m; i++) /* diagonal element was computed */
    *(*(ptrIR + i) + i) = 1.0  / *(*(ptrR + i) + i); /* same as before */
  for (int i = 0; i < m - 1; i++) /* for these indices just above the diagonal element  */
    *(*(ptrIR + i + 1) + i) = 0.0 - (   *(*(ptrR + i + 1) + i)   / (  *(*(ptrR + i) + i)  *  *(*(ptrR + i + 1) + i + 1) )  );
        /* for this computation we have to take the negative of the element divided by the two nearest diagnonal entries 
         * check the matrix calculation it is correct*/
    /* rest of the element computed by newly devised algorithm, it has to be peer reviewed, i checked twice this algorithm */
  double tbd, btd;
  double TBD[m+1], BTD[m+1];/* TBD = [1,d1,d1d2,d1d2d3,...,d1d2d3...dm ]*/
  /*BTD = [1,dm,dmdm-1,dmdm-1dm-2,....,dmdm-1...d1]  */
  tbd = 1.0;
  btd = 1.0;
  TBD[0] = 1.0;
  BTD[0] = 1.0;
  double U,D; /* now we are dealing with elements above then second diagonal of the matrix */
  for (int i = 0; i < m; i++){
    tbd *= *(*(ptrR + i) + i);
    btd *= *(*(ptrR + m - 1 - i) + m - 1 - i)
    TBD[i+1] = tbd; /* from top to bottom block diagonal multiplication */
    BTD[i+1] = btd; /* from bottom to top block diagonal multiplication */
  }
/* a11   a12 ....  a1jm1    a1jp1 ...  a1n
 *  0    a22 ....  a2jm1    a2jp1 ... .a2n 
 * .................................
 * aim11 aim12     aim1jm1  aim1jp1   aim1n
 * aip11 aip12     aip1jm1  aip1jp1   aip1n
 * ......................................
 * 0      0            0      0        ann*/


  for (int j = 0 ; j < m-2; j++){
    U = TBD[j]; /* top block matrix diagonal multiplication, [a11 a12 .... a1j-1; 0 a22 a23 .. a2j-1; 0 0 0 aj-1j-1] */
    for (int i = j + 2; i < m; i++){
      D = BTD[m-1-i]; /* bottom diagonal multiplication [ai+1i+1, ai+1i+2.....ai+1m; 0 ai+2i+2...ai+2m; 0 0 0 amm] */

/* the cofactor of the iaij will looks as follows */
      /* U */
      /* a11 a12 ..... a1jm1   |   a1jp1 a1jp2 a1jp3 .....   a1m  
       * 0   a22 ......a2jm1   |   a2jp1 a2jp2 a2jp3 .....   a2m 
       * 0   a33 ....  a2jm1   |
       * 0   0         ajm1jm1 |   ajm1jp1 ajm1jp2 ajm1jp3 .....   ajm1m  
       *
       *
       *
       *0 0 0 00      ajjp1    ajjp2 .........................................ajm
       *0 0 0 0       ajp1jp1  ajp1jp2........................................ajp1m
       *0 0 0 00        0      ajp2jp2 .......................................ajp2m
          0 0                     0   ajp3jp3 ............................... ajp3m
                                     
                                     | Im[4]        
                                     |         alpha * D * b - beta * a * D
       *                             | aim2im2| aim2im1| aim2i  aim2ip1 aim2ip2 .......aim2m   |
       *         0      0   0     0  |        | aim1im1| aim1i   aim1ip1 aim1ip2 ........aim1m | Db
       *         0              0    |             0   |   0  D| aip1ip1 aip1ip2 ........aip1m |  
       *         0              0    |             0   |   0   |    0    aip2ip2 ........aip2m |
       *         0              0    |             0   |   0   |    0      0              amm  | */



      double IM[i-j-1+3]; /* for computing the sucessive determinant of the submatrix*/
      /* [0, D, Db, P, P1, P2,...Pl] */
      IM[0] = 0;
      IM[1] = D;
      /* a = *(*(ptrR + i-1) + i - 1), b = *(*(ptrR + i) + i - 1) */
      /*alpha = *(*(ptrR + i-1) + i - 2), beta = *(*(ptrR + i) + i - 2) */
      IM[2] = D * *(*(ptrR + i) + i - 1); /* d*b */
      /* p = alpha b*d - beta * a * d */
      IM[3] = D  *  ( *(*(ptrR + i-1) + i - 2) * *(*(ptrR + i) + i - 1)  -  *(*(ptrR + i) + i - 2) *  *(*(ptrR + i-1) + i - 1) );
      /*For any such off diagnoal other than imediate off diagonal these four values has to be used to compute the determinant of the resultant
       * sub-matrix! */
      /* there will be more depends on the jth position as well!, so*/
      /* the rest of the elements depends on the difference of i and j index*/
      for (int k = 1; k < (i - j - 1); k++){
        /*complex loop has to be written */
        double gamma = 1.0;
        int l = 3 + k - 1;
        int z = 0;
        IM[3 + k] = 0.0;
        for (; l > 0; l--){
          /* alpha' = *(*(ptrR + i-2 + z) + i - 3) */
          IM[3 + k] += IM[l] * gamma * *(*(ptrR + i-2 + z) + i - 3);
          gamma = 0.0 - gamma * *(*(ptrR + i-2 + z) + i - 2 + z);
          z += 1;
        }
      }
       *(*(ptrIR + i) + j) = (double) pow(-1, (i + j)) * (U * IM[i-j+1]) / det;
    }
  }

}/* need to review this code but this is hell out of complex matrix inversion for the given triangular matrix form to get faster and better numerical stable solutions 
    , not using the conventional matrix  inversion!*/
/* it is reviewed once, but not sure that this will reduce the numerical inaccuracy, since division(by a big integer) and multiplication (by real
 * number) leads the same numerical inaccuracy I guess so */



void NSing_UpTriag_Mat_Inv_Adj_Abdul_Mthd_Fl(float **ptrR, float **ptrIR, int m){
  float det = 1.0;
  for (int i = 0; i < m; i++)
    det *= *(*(ptrR + i) + i);
  for (int i = 0; i < m; i++) /* diagonal element was computed */
    *(*(ptrIR + i) + i) = 1.0  / *(*(ptrR + i) + i); /* same as before */
  for (int i = 0; i < m - 1; i++) /* for these indices just above the diagonal element  */
    *(*(ptrIR + i + 1) + i) = 0.0 - (   *(*(ptrR + i + 1) + i)   / (  *(*(ptrR + i) + i)  *  *(*(ptrR + i + 1) + i + 1) )  );
        /* for this computation we have to take the negative of the element divided by the two nearest diagnonal entries 
         * check the matrix calculation it is correct*/
    /* rest of the element computed by newly devised algorithm, it has to be peer reviewed, i checked twice this algorithm */

  double tbd, btd;
  double TBD[m+1], BTD[m+1];/* TBD = [1,d1,d1d2,d1d2d3,...,d1d2d3...dm ]*/
  /*BTD = [1,dm,dmdm-1,dmdm-1dm-2,....,dmdm-1...d1]  */
  tbd = 1.0;
  btd = 1.0;
  TBD[0] = 1.0;
  BTD[0] = 1.0;
  float U,D; /* now we are dealing with elements above then second diagonal of the matrix */

  for (int i = 0; i < m; i++){
    tbd *= *(*(ptrR + i) + i);
    btd *= *(*(ptrR + m - 1 - i) + m - 1 - i)
    TBD[i+1] = tbd; /* from top to bottom block diagonal multiplication */
    BTD[i+1] = btd; /* from bottom to top block diagonal multiplication */
  }

  for (int j = 0 ; j < m-2; j++){
    U = TBD[j]; /* top block matrix diagonal multiplication, [a11 a12 .... a1j-1; 0 a22 a23 .. a2j-1; 0 0 0 aj-1j-1] */
    for (int i = j + 2; i < m; i++){
      D = BTD[m-1-i]; /* bottom diagonal multiplication [ai+1i+1, ai+1i+2.....ai+1m; 0 ai+2i+2...ai+2m; 0 0 0 amm] */

      float IM[i-j-1+3]; /* for computing the sucessive determinant of the submatrix*/
      /* [0, D, Db, P, P1, P2,...Pl] */
      IM[0] = 0;
      IM[1] = D;
      /* a = *(*(ptrR + i-1) + i - 1), b = *(*(ptrR + i) + i - 1) */
      /*alpha = *(*(ptrR + i-1) + i - 2), beta = *(*(ptrR + i) + i - 2) */
      IM[2] = D * *(*(ptrR + i) + i - 1); /* d*b */
      /* p = alpha b*d - beta * a * d */
      IM[3] = D  *  ( *(*(ptrR + i-1) + i - 2) * *(*(ptrR + i) + i - 1)  -  *(*(ptrR + i) + i - 2) *  *(*(ptrR + i-1) + i - 1) );
      /*For any such off diagnoal other than imediate off diagonal these four values has to be used to compute the determinant of the resultant
       * sub-matrix! */
      /* there will be more depends on the jth position as well!, so*/
      /* the rest of the elements depends on the difference of i and j index*/
      for (int k = 1; k < (i - j - 1); k++){
        /*complex loop has to be written */
        float gamma = 1.0;
        int l = 3 + k - 1;
        int z = 0;
        IM[3 + k] = 0.0;
        for (; l > 0; l--){
          /* alpha' = *(*(ptrR + i-2 + z) + i - 3) */
          IM[3 + k] += IM[l] * gamma * *(*(ptrR + i-2 + z) + i - 3);
          gamma = 0.0 - gamma * *(*(ptrR + i-2 + z) + i - 2 + z);
          z += 1;
        }
      }
       *(*(ptrIR + i) + j) = (float) pow(-1, (i + j)) * (U * IM[i-j+1]) / det;
    }
  }

}/* reviewed this code once,  but this is hell out of complex matrix inversion not using the conventional matrix  inversion! 
not sure it increase the numerical stability*/


