#ifndef _BFS_LIB
#define _BFS_LIB

#include <stdio.h>
#include <stdlib.h>
#include "network_lib.h"
//#include "globals.h"


/////////ARBRE BFS
struct node_bfs{
    int d;
    struct node_gra *ref;    // node it represents 
    struct node_bfs *next;   // Next element in the list
    struct node_bfs *prev;   // Previous element in the list
    struct node_bfs *last;   // Last element in the list
    struct pred *pred;       // Predecessor for bfs algorithm
};

struct pred{
    struct node_bfs *ref;
    struct pred *next;
};



//BFS TREE
struct node_bfs *CreateHeaderList();
struct node_bfs *GetBFS(struct node_gra *node,struct node_bfs *list);

struct node_bfs *RenewQueue(struct node_bfs *list,struct node_bfs *lp,int *size,int d);
struct node_gra *Dequeue(struct node_bfs *list,int *size);
struct node_gra *DequeueOne(struct node_bfs *list,struct node_bfs *one,int *size);
void Enqueue(struct node_gra *node,struct node_bfs *predecessor,struct node_bfs *header,int *size,int dist);

void AddPredecessor(struct node_bfs *node,struct node_bfs *pred);
void ClearPredecessors(struct pred *p);
int CountPredecessors(struct node_bfs *node);

void ClearList(struct node_bfs *list,int *size);

void FPrintCurrentList(FILE *file,struct node_bfs *p); 
void PrintCurrentList(struct node_bfs *p); 

int IsGraphConnected(struct node_gra *p,int N);

#endif