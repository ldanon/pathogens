#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#include "bfs_lib.h"
#include <math.h>
#include "globals.h"
//FUNCIONS BFS
////////////////////////////////////////////////////////////////////////
//creates a node_bfs for use as header
struct node_bfs *CreateHeaderList() 
{
    struct node_bfs *temp;
    
    temp=(struct node_bfs *)malloc(sizeof(struct node_bfs));
    temp->d=-1;
    temp->ref=NULL;
    temp->next=NULL;
    temp->prev=NULL;
    temp->last=NULL;
    temp->pred=NULL;
    
    return temp;
}


////////////////////////////////////////////////////////////////////////
// find bfs node "*list" with reference to "*node"
struct node_bfs *GetBFS(struct node_gra *node,struct node_bfs *list)
{
    while(list!=NULL){
        if(list->ref==node) return list;
        list=list->next;
    }
    return list;
}



////////////////////////////////////////////////////////////////////////
void Enqueue(struct node_gra *node,struct node_bfs *predecessor,
             struct node_bfs *header,int *size,int dist)
{
    // *node= node you want to put in bfs list
    // *predecessor=
    // *header of predecessor list
    // *size = size of list
    // dist = distance from node.
    
    struct node_bfs *temp;
    
    if(node->state==0){    //if node hasn't been "queued" create *temp
        temp=(struct node_bfs *)malloc(sizeof(struct node_bfs));
        temp->d=dist;
        temp->next=NULL;
        temp->ref=node;
        temp->last=NULL;
        temp->pred=(struct pred *)malloc(sizeof(struct pred));
        (temp->pred)->ref=predecessor;
        (temp->pred)->next=NULL;
        *size+=1;
        node->state=1;
        
        // Actualitzem l'apuntador a l'ultim element de la llista
        // Refresh the pointer to the last element in the list
        //I enllacem el node al darrere de l'ultim
        
        if(header->last==NULL){
            header->next=header->last=temp;
            temp->prev=header;
        }else{
            temp->prev=header->last;
            (header->last)->next=temp;
            header->last=temp;
        }
    }
    //Si el node ja es dins la cua, mirem si l'actual es un predecessor seu
    //aixo passa si la distancia del cami actual es igual a la que tenim.
    else{
        temp=GetBFS(node,predecessor);
        if(temp!=NULL)
            if(temp->d==dist)
                AddPredecessor(temp,predecessor);
    }
}


////////////////////////////////////////////////////////////////////////
//Encuem tots els fills dels nodes q estan a distancia d
struct node_bfs *RenewQueue(struct node_bfs *list,struct node_bfs *lp,int *size,int d)
{
    struct node_lis *p;
    
    while((lp->next!=NULL)&&((lp->next)->d==d)){
        lp=lp->next;
        p=(lp->ref)->neig;		
        while(p->next!=NULL){
            p=p->next;
            Enqueue(p->ref,lp,list,size,d+1);
        }
    }
    return lp;//
}


////////////////////////////////////////////////////////////////////////
//add predecessor "*pred" to the end of the list of predecessors of "*node"
void AddPredecessor(struct node_bfs *node,struct node_bfs *pred)
{
    struct pred *p=node->pred;
    
    while(p->next!=NULL) p=p->next;
    
    p->next=(struct pred *)malloc(sizeof(struct pred));
    (p->next)->ref=pred;
    (p->next)->next=NULL;
}


////////////////////////////////////////////////////////////////////////
// clear all predecessors beginning at p and remove p from memory
void ClearPredecessors(struct pred *p)  
{                                  
    if(p->next!=NULL)
        ClearPredecessors(p->next);
    free(p);
}

////////////////////////////////////////////////////////////////////////
void ClearList(struct node_bfs *list,int *size)
{
    struct node_gra *temp=NULL;
    
    while(list->next!=NULL)
        temp=Dequeue(list,size);
}


////////////////////////////////////////////////////////////////////////
struct node_gra *Dequeue(struct node_bfs *list,int *size)
{
    struct node_bfs *temp; 
    struct node_gra *node;
    
    temp=list->next;            //save list->next in temp 
    list->next=temp->next;      //save temp->next in list->next
    if(list->next==NULL)        //if *list is last node in list
        list->last=NULL;          //remove *list and set last to NULL
    else                        
        (list->next)->prev=list;  //if not, set prev to list
    // that is point to previous in list
    
    *size-=1;
    
    node=temp->ref;             // return node_gra that original *list pointed to
    
    ClearPredecessors(temp->pred);
    free(temp);
    
    return node;
}


////////////////////////////////////////////////////////////////////////
struct node_gra *DequeueOne(struct node_bfs *list,struct node_bfs *one,int *size)
{
    struct node_bfs *temp;
    struct node_bfs *p;
    struct node_gra *node;
    
    temp=one->next;
    
    if(list->last==temp){
        p=list;
        while(p->next!=NULL){
            p->last=one;
            p=p->next;
        }
    }
    
    if(temp->next==NULL)  one->next=NULL;
    else  one->next=temp->next;
    
    *size-=1;
    
    node=temp->ref;
    
    ClearPredecessors(temp->pred);
    free(temp);
    
    return node;
}





////////////////////////////////////////////////////////////////////////
int CountPredecessors(struct node_bfs *node)
{
    struct pred *p=node->pred;
    int counter=0;
    
    if((p->ref)->ref!=NULL){  // If the predecessor is the header dont count
        while(p!=NULL){
            counter++;
            p=p->next;
        }
    }
    
    return counter;
}

////////////////////////////////////////////////////////////////////////
// print node_bfs list to file
void FPrintCurrentList(FILE *file,struct node_bfs *p) 
{
    while(p->next!=NULL){
        p=p->next;
        fprintf(file,"    %d (%d)\n",(p->ref)->num,p->d);
    }
}


void PrintCurrentList(struct node_bfs *p) 
{
    while(p->next!=NULL){
        p=p->next;
        printf("    %d (%d)\n",(p->ref)->num,p->d);
    }
}

int IsGraphConnected( node_gra *p,int N)
{
  //returns 1 if it is connected, and 0 if not.

  node_bfs *list,*lp;
  int *size,r1;
  int size_ant;
  int d;

  r1=0;
  size=&r1;

  ResetNodesState(p);
  p=p->next;
  list=CreateHeaderList();

  Enqueue(p,list,list,size,0);
  lp=list;
  d=0;
  do{
    size_ant=*size;
    lp=RenewQueue(list,lp,size,d);
    d++;
  }while(*size!=size_ant);

  ClearList(list,size);

  if(size_ant==N)
    return 1;
  else
    return 0;
}


