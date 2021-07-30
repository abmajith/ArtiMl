#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* lye denotes the number of edge layers will be less then one of the actual layer*/

/* nd is the pointer to the nodes in each layer i.e 
 * nd points -> [n0,n1,n2,....nlye]  n0 is the number of input nodes, nlye is the number of outputs nodes in the last layer*/

/* for each edge layer ni is the number of input nodes and ni+1 is the number of output nodes */



/* **ptrMw is an array of pointer to pointer of weights, *(*(ptrMw[l] + i) + j) is the weight wji from the ith node in layer l to jth node in
 * layer lplus1 */

/* ptrMw[l]  represents all possible lth layer edges */
/* *(ptrMw[l] + i) represents all possible edges from the node i in the layer lth layer to the lplu1th layer  */
/* when I say lth layer edges is different from lth layer later represents the node layer  */
/* *(*(ptrMw[l] + i) + j) is the edge from ith node in  lth layer  to the jth node in lplus1th layer*/


/* *ptrMt[] is the model threshhold for each node and for input layer there is no threshold! */
/* to reduce the space we dont have to store the threshold limits for the input layer!, because there is no thereshold */
/* ptrMt[0] -> [t1,t2,...tk0] t1 is the threshold for the node 1 in the layer 1 and so on*/
/* ptrMt[l] -> [t1,t2,...tkl] ti the theshold for the node i in the layer l and so on*/


/* ptrIn is the array to the pointers */
/* ptrIn[0] points -> [x1,x2,..xn] inputs for the input layer*/
/* ptrIn[lye] points -> [o1,o2,..om] outputs for the current model parameter (weigths)*/
/* note, this m should be equal to nlye */

/* ptrIn[0] -> [x1, x2, x3,...x[n0]] */
/* ptrIn[1] -> [o1,o2, o3,....o[n1]] */
/* ptrIn[2] -> [o1,o2,....o[n2]]*/
/*................ */
/* ptrIn[lye] -> [o1,o2,...,o[nlye]]*/

/* t is the point to the acutal output for the given input, and its points to  -> [t1,t2,....,t[nlye]] */

/* *del[] will be same as ptrIn except there is not array for the input nodes, since we dont need to calcuate such values i.e */
/* del[0] = [del1,del2,.....del[n1]] */
/* del[1] = [del1,del2,.....del[n2]] */


int NNlayer_Logfn_input_Db(double **ptrMw[], double *ptrMt[], int lye, int *nd, double *ptrIn[], double *t, double *del[], double eta ){

  double net;
  double de;

  /* we will compute the inputs for various layers, and assume that for the input layer its already written in the ptrIn*/
  for (int l = 0; l < lye; l++)
    for (int j = 0; j < *(nd + l + 1); j++){
      net = 0.0;
      for (int i = 0; i < *(nd + l); i++){
        net += (double) *(*(ptrMw[l] + i) + j) *  (double) *(ptrIn[l] + i);
      }
      *(ptrIn[l + 1] + j) = (double) 1.0 / (1.0 + (double) exp( *(ptrMt[l] + j) - net ) ); /*logistic functions */
    }

  /* computing the delta for each layer*/
  for (int j = 0; j < *(nd + lye + 1); j++) /* for output the calucation is supervised so!*/
    *(del[lye - 1] + j) = ( *(t + j) - *(ptrIn[lye] + j) )  *  *(ptrIn[lye] + j)   *  ( 1.0 - *(ptrIn[lye] + j) );
  for (int l = lye - 1; l > 0; l++)
    for (int j = 0; j < *(nd + l); j++){
      de = 0.0;
      for (int k = 0; k < *(nd + l + 1); k++)
        de += *(*(ptrMw[l] + j) + k) * *(del[l] + k);
      *(del[l - 1] + j) = de *  *(ptrIn[l] + j)   * ( 1.0 - *(ptrIn[l] + j)  );
    }

  /*updating the weigths now using the delta, and calucated input for the each layers */
  for (int l = 0; l < lye; l++)
    for (int j = 0; j < *(nd + l + 1); j++)
      for (int i = 0; i < *(nd + l); i++)
        *(*(ptrMw[l] + i) + j)  += eta * *(del[l] + j) * *(ptrIn[l] + i); /* wjinew = wji0ld + eta deltaj input at the ith node. */
  return 0;
}


int NNlayer_Logfn_input_Fl(float **ptrMw[], float *ptrMt[], int lye, int *nd, float *ptrIn[], float *t, float *del[], float eta ){

  float net;
  float de;

  /* we will compute the inputs for various layers, and assume that for the input layer its already written in the ptrIn*/
  for (int l = 0; l < lye; l++)
    for (int j = 0; j < *(nd + l + 1); j++){
      net = 0.0;
      for (int i = 0; i < *(nd + l); i++){
        net += (float) *(*(ptrMw[l] + i) + j) *  (float) *(ptrIn[l] + i);
      }
      *(ptrIn[l + 1] + j) = (float) 1.0 / (1.0 + (float) exp( *(ptrMt[l] + j) - net ) ); /*logistic functions */
    }

  /* computing the delta for each layer*/
  for (int j = 0; j < *(nd + lye + 1); j++) /* for output the calucation is supervised so!*/
    *(del[lye - 1] + j) = ( *(t + j) - *(ptrIn[lye] + j) )  *  *(ptrIn[lye] + j)   *  ( 1.0 - *(ptrIn[lye] + j) );
  for (int l = lye - 1; l > 0; l++)
    for (int j = 0; j < *(nd + l); j++){
      de = 0.0;
      for (int k = 0; k < *(nd + l + 1); k++)
        de += *(*(ptrMw[l] + j) + k) * *(del[l] + k);
      *(del[l - 1] + j) = de *  *(ptrIn[l] + j)   * ( 1.0 - *(ptrIn[l] + j)  );
    }

  /*updating the weigths now using the delta, and calucated input for the each layers */
  for (int l = 0; l < lye; l++)
    for (int j = 0; j < *(nd + l + 1); j++)
      for (int i = 0; i < *(nd + l); i++)
        *(*(ptrMw[l] + i) + j)  += eta * *(del[l] + j) * *(ptrIn[l] + i); /* wjinew = wji0ld + eta deltaj input at the ith node. */
  return 0;
}



