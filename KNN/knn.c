#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include "knn.h"


// Note: if the number of datas are too much, consider to use the fopen functions and process the data by chunck by chunk
int Normalize_KNN_vec_data_Db_prec(int vec_dim, int No_data_pts, double *data[], double mean[], double std[]){
  if (No_data_pts <= 0){
    printf("There are no data points there to normalize, exiting the program\n");
    exit(0);
  }
  if (vec_dim <= 0){
    printf("Input datas vector dimenstions are null, exiting the program\n");
    exit(0);
  }
  double Sum[vec_dim];
  double Sqsum[vec_dim];
  for (int i = 0; i < vec_dim; i++){
    Sum[i] = 0.0;
    Sqsum[i] = 0.0;
  }
  for (int j = 0; j < No_data_pts; j++){
    for (int i = 0; i < vec_dim; i++){
      Sum[i] += *(data[j] + i);
      Sqsum[i] += ( *(data[j] + i)  *  *(data[j] + i));
    }
  }
  for (int i = 0; i < vec_dim; i++){
    mean[i] = Sum[i] / No_data_pts;
    Sqsum[i] /= No_data_pts;
    std[i] = (double) sqrt(Sqsum[i] - ( mean[i] * mean[i] ) );
  }
  return 0;
}


double biggervalue(Node *Tree){
  if (Tree->right == NULL)
    return Tree->value;
  else 
    return biggervalue(Tree->right);
}


Node* deletebiggernode(Node *Tree){
  if ( (Tree->right == NULL) && (Tree->left == NULL) ){
    free(Tree);
    return NULL;
  } 
  if (Tree->right == NULL){
    Node *temp = Tree->left;
    free(Tree);
    return temp;
  }
  if (Tree->right != NULL) {
    Tree->right = deletebiggernode(Tree->right);
    return Tree;
  }
  return Tree;
}


Node* addnode(Node *Tree, Node *New){
  if (Tree->value >= New->value ){
    if (Tree->left == NULL){
      Tree->left = New;
    } else {
      addnode(Tree->left, New);
    }
  } else {
    if (Tree->right == NULL){
      Tree->right = New;
    } else {
      addnode(Tree->right, New);
    }
  }
  return Tree;
}


LLnode* update_linkedarray(LLnode *Ary, int x){
  if (Ary != NULL){
    if (Ary->label == x){
      Ary->lb_count += 1;
      return Ary;
    } else{
      Ary->next = update_linkedarray(Ary->next, x);
      return Ary;
    }
  }else {
    LLnode *new = NULL;
    new = (LLnode*) malloc(sizeof(LLnode));
    new->lb_count = 1;
    new->label = x;
    new->next = NULL;
    return new;
  }
}

LLnode* find_label_array(Node *Tree, LLnode *Ary){
  if (Tree != NULL){
    Ary = update_linkedarray(Ary, Tree->label);
  }
  if (Tree->left != NULL){
    Ary = find_label_array(Tree->left, Ary);
  }
  if (Tree->right != NULL){
    Ary = find_label_array(Tree->right, Ary);
  }
  return Ary;
}

int getargmaxfromlist(LLnode *Ary){
  int argmax;
  int max = 0;
  while (Ary != NULL){
    if (max < Ary->lb_count){
      max = Ary->lb_count;
      argmax = Ary->label;
    }
    Ary = Ary->next;
  }
  return argmax;
}


int querry_KNN_maj_voting_exhaustive(int vec_dim, int No_data_pts, double *data[], int data_label[], int querry_data[], int k, int *querry_result_label){
  if ( k <= 0){
    printf("Impossible to search for the %d neighbours, exiting the program\n",k);
    exit(0);
  }
  if ( No_data_pts <= 0){
    printf("There are no data points to serach for the nearest %d neighbours, exiting the program\n",k);
    exit(0);
  }
  if (vec_dim <= 0 ){
    printf("Input datas vector are dimenstions are null, exiting the program\n");
    exit(0);
  }
  if (k > No_data_pts){
    printf("Number of data points less than %d exiting the program\n",k);
    exit(0);
  }
  Node *tree = NULL;
  tree = (Node*) malloc(sizeof(Node));
  tree->value = DBL_MAX;
  tree->label = data_label[0];
  tree->left = NULL;
  tree->right = NULL;
  double cur_dis = 0.0;
  double cur_min;
  double cur_maxk;
  for (int j = 0; j < vec_dim; j++){
    cur_dis += ( *(data[0] + j) - querry_data[j] ) * ( *(data[0] + j) - querry_data[j] );
    cur_min = cur_dis;
    cur_maxk = cur_dis;
  }
  if (k == 1){
    int argmin = data_label[0];
    for (int i = 1; i < No_data_pts; i++){
      cur_dis = 0.0;
      for (int j = 0; j < vec_dim; j++){
        cur_dis += (  *(data[i] + j)  - querry_data[j] ) * (  *(data[i] + j) - querry_data[j] );
      }
      if (cur_dis < cur_min){
        cur_min = cur_dis;
        argmin = data_label[i];
      }
    }
    *querry_result_label = argmin;
    return 0;
  }
  tree->value = cur_dis;
  for (int i = 1; i < k; i++){
    cur_dis = 0.0;
    for (int j = 0; j < vec_dim; j++){
      cur_dis += (  *(data[i] + j)  - querry_data[j] ) * (  *(data[i] + j) - querry_data[j] );
    }
    if (cur_dis < cur_min){
      cur_min = cur_dis;
    }
    Node *newnode = NULL;
    newnode = (Node*) malloc(sizeof(Node));
    newnode->value = cur_dis;
    newnode->label = data_label[i];
    newnode->left = NULL;
    newnode->right = NULL;
    addnode(tree,newnode);
    if (cur_maxk < cur_dis){
      cur_maxk = cur_dis;
    }
  }

  for (int i = k; i < No_data_pts; i++){
    cur_dis = 0.0;
    for (int j = 0; j < vec_dim; j++){
      cur_dis += (  *(data[i] + j)  - querry_data[j] ) * (  *(data[i] + j) - querry_data[j] );
    }
    if (cur_dis < cur_min){
      cur_min = cur_dis;
    }
    if (cur_maxk > cur_dis){
      tree = deletebiggernode(tree);
      Node *newnode = NULL;
      newnode = (Node*) malloc(sizeof(Node));
      newnode->value = cur_dis;
      newnode->label = data_label[i];
      newnode->left = NULL;
      newnode->right = NULL;
      addnode(tree, newnode);
      cur_maxk = biggervalue(tree);
    }
  }

  LLnode *leftarray = NULL;
  leftarray = find_label_array(tree, leftarray);
  *querry_result_label = getargmaxfromlist(leftarray);
  return 0;
}
