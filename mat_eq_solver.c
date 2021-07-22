#include <math.h> /*i am using sqrt function, so its needed */

/* Given non singular upper trianglular square matrix m by m, find I */
/* By Backward substitution method */
/* note this one can lead severe numerical approximation inaccuracy, be caution with it!
 * its not because the code is badly written, its the draw back of backward substituion method which involves 
 * sequences of many division operation*/
void NSing_UpTriag_Mat_Inv_Backsub_Mthd_Db(double **ptrR, double **ptrIR, int m){
  double det = 1.0;
  for (int i = 0; i < m; i++) /* determinant of the upper triangular matrix is computed */
    det *= *(*(ptrR + i) + i);
  for (int i = 0; i < m; i++) /* diagonal of inverse of R matrix is computed */
    *(*(ptrIR + i) + i) = *(*(ptrR + i) + i)  / det;
  double init;
  /* backward substitution method */
  for (int z = m - 1; z > 0; z --) /* computing the zth row value */
    for (int i = z - 1; i >= 0; i--){ /* computing inverse of zth row and i the column value */
      init = 0.0;
      for (int j = i + 1; j < z; j++)
        init += *(*(ptrR + j ) + i ) * *(*(ptrIR + z) + j);
      *(*(ptrIR + z) + i) = 0.0 - (double) (init / *(*(ptrR + i) + i);
    } /* This block is checked once, check another time later*/
  return;
}

void NSing_UpTriag_Mat_Inv_Backsub_Mthd_Fl(float **ptrR, float **ptrIR, int m){
  float det = 1.0;
  for (int i = 0; i < m; i++) /* determinant of the upper triangular matrix is computed */
    det *= *(*(ptrR + i) + i);
  for (int i = 0; i < m; i++) /* diagonal of inverse of R matrix is computed */
    *(*(ptrIR + i) + i) = *(*(ptrR + i) + i)  / det;
  float init;
  /* backward substitution method */
  for (int z = m - 1; z > 0; z --) /* computing the zth row value */
    for (int i = z - 1; i >= 0; i--){ /* computing inverse of zth row and i the column value */
      init = 0.0;
      for (int j = i + 1; j < z; j++)
        init += *(*(ptrR + j ) + i ) * *(*(ptrIR + z) + j);
      *(*(ptrIR + z) + i) = 0.0 - (float) (init / *(*(ptrR + i) + i);
    } /* This block is checked once, check another time later*/
  return;
}

/* hear we will write more extensive, stable adjoint method rather than backward substitution method to get better numerical stability!*/
void NSing_UpTriag_Mat_Inv_Adj_Abdul_Mthd_Db(double **ptrR, double **ptrIR, int m){
  double det = 1.0;
  for (int i = 0; i < m; i++)
    det *= *(*(ptrR + i) + i);
  for (int i = 0; i < m; i++) /* diagonal element was computed */
    *(*(ptrIR + i) + i) = *(*(ptrR + i) + i)  / det;
  for (int i = 0; i < m - 1; i++) /* for these indices just above the diagonal element  */
    *(*(ptrIR + i + 1) + i) = 0.0 - (   *(*(ptrR + i + 1) + i)   / (  *(*(ptrR + i) + i) * *(*(ptrR + i + 1) + i + 1) )  );
    /* rest of the element computed by newly devised algorithm, it has to be peer reviewed, i checked twice this algorithm */
  double tbd, btd;
  double TBD[m+1], BTD[m+1];/* TBD = [1,d1,d1d2,d1d2d3,...,d1d2d3...dm ]*/
  /*BTD = [1,dm,dmdm-1,dmdm-1dm-2,....,dmdm-1...d1]  */
  tbd = 1.0;
  btd = 1.0;
  TBD[0] = 1.0;
  BTD[0] = 1.0;
  double U,D;
  for (int i = 0; i < m; i++){
    tbd *= *(*(ptrR + i) + i);
    btd *= *(*(ptrR + m - 1 - i) + m - 1 - i)
    Tbd[i+1] = tbd;
    BTD[i+1] = btd;
  }
  for (int j = 0 ; j < m-2; j++){
    U = TBD[j]; /* top diagonal multiplication */
    for (int i = j + 2; i < m; i++){
      D = BTD[m-1-i]; /* bottom diagonal multiplication */
      double IM[i-j-1+3]; /* for computing the sucessive determinant of the submatrix*/
      /* [0, D, Db, P, P1, P2,...Pl] */
      IM[0] = 0;
      IM[1] = D;
      /* a = *(*(ptrR + i-1) + i - 1), b = *(*(ptrR + i) + i - 1) */
      /*alpha = *(*(ptrR + i-1) + i - 2), beta = *(*(ptrR + i) + i - 2) */
      IM[2] = D * *(*(ptrR + i) + i - 1); /* d*b */
      /* p = alpha b*d - beta * a * d */
      IM[3] = D * ( *(*(ptrR + i-1) + i - 2) * *(*(ptrR + i) + i - 1)  - *(*(ptrR + i) + i - 2) *  *(*(ptrR + i-1) + i - 1) );
      /*For any such off diagnoal other than imediate off diagonal these four values has to be used to compute the determinant of the resultant
       * matrix! */
      /* there will be more depends on the jth position as well!, so*/
      for (int k = 1; k < (i - j - 1); k++){
        /*complex loop has to be written */
        double gamma = 1.0;
        int l = 3 + k - 1;
        int z = 0, z1 = 1;
        IM[3 + k] = 0.0;
        for (; l > 0; l--){
          /* alpha' = *(*(ptrR + i-2 + z) + i - 3) */
          IM[3 + k] += IM[l] * gamma * *(*(ptrR + i-2 + z) + i - 3);
          z += 1;
          gamma = 0.0 - gamma * *(*(ptrR + i-2 + z) + i - 3 - z1);
          z1 += 1;
        }
      }
       *(*(ptrIR + i) + j) = (double) pow(-1, (i + j)) * (U * IM[i-j+1]) / det;
    }
  }

}/* need to review this code but this is hell out of complex matrix inversion for the given triangular matrix form to get faster and better numerical stable solutions 
    , not using the conventional matrix  inversion!*/



void NSing_UpTriag_Mat_Inv_Adj_Abdul_Mthd_Fl(float **ptrR, float **ptrIR, int m){
  float det = 1.0;
  for (int i = 0; i < m; i++)
    det *= *(*(ptrR + i) + i);
  for (int i = 0; i < m; i++) /* diagonal element was computed */
    *(*(ptrIR + i) + i) = *(*(ptrR + i) + i)  / det;
  for (int i = 0; i < m - 1; i++) /* for these indices just above the diagonal element  */
    *(*(ptrIR + i + 1) + i) = 0.0 - (   *(*(ptrR + i + 1) + i)   / (  *(*(ptrR + i) + i) * *(*(ptrR + i + 1) + i + 1) )  );
    /* rest of the element computed by newly devised algorithm, it has to be peer reviewed, i checked twice this algorithm */
  float tbd, btd;
  float TBD[m+1], BTD[m+1];/* TBD = [1,d1,d1d2,d1d2d3,...,d1d2d3...dm ]*/
  /*BTD = [1,dm,dmdm-1,dmdm-1dm-2,....,dmdm-1...d1]  */
  tbd = 1.0;
  btd = 1.0;
  TBD[0] = 1.0;
  BTD[0] = 1.0;
  float U,D;
  for (int i = 0; i < m; i++){
    tbd *= *(*(ptrR + i) + i);
    btd *= *(*(ptrR + m - 1 - i) + m - 1 - i)
    Tbd[i+1] = tbd;
    BTD[i+1] = btd;
  }
  for (int j = 0 ; j < m-2; j++){
    U = TBD[j]; /* top diagonal multiplication */
    for (int i = j + 2; i < m; i++){
      D = BTD[m-1-i]; /* bottom diagonal multiplication */
      float IM[i-j-1+3]; /* for computing the sucessive determinant of the submatrix*/
      /* [0, D, Db, P, P1, P2,...Pl] */
      IM[0] = 0;
      IM[1] = D;
      /* a = *(*(ptrR + i-1) + i - 1), b = *(*(ptrR + i) + i - 1) */
      /*alpha = *(*(ptrR + i-1) + i - 2), beta = *(*(ptrR + i) + i - 2) */
      IM[2] = D * *(*(ptrR + i) + i - 1); /* d*b */
      /* p = alpha b*d - beta * a * d */
      IM[3] = D * ( *(*(ptrR + i-1) + i - 2) * *(*(ptrR + i) + i - 1)  - *(*(ptrR + i) + i - 2) *  *(*(ptrR + i-1) + i - 1) );
      /*For any such off diagnoal other than imediate off diagonal these four values has to be used to compute the determinant of the resultant
       * matrix! */
      /* there will be more depends on the jth position as well!, so*/
      for (int k = 1; k < (i - j - 1); k++){
        float gamma = 1.0;
        int l = 3 + k - 1;
        int z = 0, z1 = 1;
        IM[3 + k] = 0.0;
        for (; l > 0; l--){
          /* alpha' = *(*(ptrR + i-2 + z) + i - 3) */
          IM[3 + k] += IM[l] * gamma * *(*(ptrR + i-2 + z) + i - 3);
          z += 1;
          gamma = 0.0 - gamma * *(*(ptrR + i-2 + z) + i - 3 - z1);
          z1 += 1;
        }
      }
       *(*(ptrIR + i) + j) = (float) pow(-1, (i + j)) * (U * IM[i-j+1]) / det;
    }
  }

}/* need to review this code but this is hell out of complex matrix inversion not using the conventional matrix  inversion!*/


