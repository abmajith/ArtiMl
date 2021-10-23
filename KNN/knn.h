#ifndef knn
#define knn

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <limits.h>


typedef struct Node {
  double value;
  int label;
  struct Node *left;
  struct Node *right;
} Node;

typedef struct LLnode {
  int lb_count;
  int label;
  struct LLnode *next;
} LLnode;



double biggervalue(Node *Tree);

Node* deletebiggernode(Node *Tree);

Node* addnode(Node *Tree, Node *New);

LLnode* update_linkedarray(LLnode *Ary, int x);

LLnode* find_label_array(Node *Tree, LLnode *Ary);

int getargmaxfromlist(LLnode *Ary);

int Normalize_KNN_vec_data_Db_prec(int vec_dim, int No_data_pts, double *data[], double mean[], double std[]);

int querry_KNN_maj_voting_exhaustive(int vec_dim, int No_data_pts, double *data[], int data_label[], int querry_data[], int k, int *querry_result_label);


#endif
