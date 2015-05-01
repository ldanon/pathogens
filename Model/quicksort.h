#ifndef _QUICKSORT_H
#define _QUICKSORT_H

#include "network_lib.h"

static int partition_link_list(struct node_lis *x[], int i, int j);
void quicksort_link_list(struct node_lis *x[], int l, int r);
struct node_lis *QuickSortLinkList(struct node_lis *x[],int n);

static int partition_fitness_list(struct node_gra *x[], int i, int j);
void quicksort_fitness_list(struct node_gra *x[], int l, int r);
struct node_gra *QuickSortFitness(struct node_gra *net, int nnodes);

static int partition_degree_list(struct node_gra *x[], int i, int j);
void quicksort_degree_list(struct node_gra *x[], int l, int r);
struct node_gra *QuickSortDegree(struct node_gra *net, int nnodes);

static int partition_DegreeInv_list(struct node_gra *x[], int i, int j);
void quicksort_DegreeInv_list(struct node_gra *x[], int l, int r);
struct node_gra *QuickSortDegreeInv(struct node_gra *net, int nnodes);

static int partition_num_list(struct node_gra *x[], int i, int j);
void quicksort_num_list(struct node_gra *x[], int l, int r);
struct node_gra *QuickSortNum(struct node_gra *net, int nnodes);

#endif