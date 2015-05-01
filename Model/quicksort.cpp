#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#include "quicksort.h"

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static int partition_link_list(struct node_lis *x[], int i, int j) {
    
    int k = j; 
    struct node_lis *v;
    struct node_lis *t;
    
    v= x[k];				//Seleccionem un element de la taula
    i--;
    while (i < j) {
        while (x[++i]->node < v->node);
        while (--j>i && x[j]->node > v->node);
        if (i < j) {
            t=x[i];	/* swap x[i] and x[j] */
            x[i] = x[j];
            x[j] = t;
            //printf("---%d\n",x[j]->label);
        }
    }
    x[k] = x[i];	/* swap x[i] and the pivot element */
    x[i] = v;
    return i;
}



//////////////////////////////////////////////////////////////////////
void quicksort_link_list(struct node_lis *x[], int l, int r) {
    //use : quicksort(list_array,1,n) where n is the number of elements
    
    if (r > l) {
        int k = partition_link_list(x, l, r); /*  */
        quicksort_link_list(x, l, k - 1);
        quicksort_link_list(x, k + 1, r);
    }
}



///////////////////////////////////////////////////////////////////////
struct node_lis *QuickSortLinkList(struct node_lis *x[],int n){
    int i;
    quicksort_link_list(x,1,n);
    for(i=1;i<n;i++){
        x[i]->next=x[i+1];
    }
    x[n]->next=NULL;
    return x[1];
}





















//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static int partition_fitness_list(struct node_gra *x[], int i, int j) {
    
    int k = j; 
    struct node_gra *v;
    struct node_gra *t;
    
    v=x[k];
    i--;
    while (i < j) {
        
        while (x[++i]->fitness < v->fitness);
        while (--j>i && x[j]->fitness > v->fitness);
        
        if (i < j) {
            t=x[i];	// swap x[i] and x[j] 
            x[i] = x[j];
            x[j] = t;
            //printf("---%d\n",x[j]->label);
        }
    }
    x[k] = x[i];	// swap x[i] and the pivot element 
    x[i] = v;
    return i;
}



//////////////////////////////////////////////////////////////////////
void quicksort_fitness_list(struct node_gra *x[], int l, int r) {
    //use : quicksort(list_array,1,n) where n is the number of elements
    
    if (r > l) {
        int k = partition_fitness_list(x, l, r);  
        quicksort_fitness_list(x, l, k - 1);
        quicksort_fitness_list(x, k + 1, r);
    }
}



///////////////////////////////////////////////////////////////////////
struct node_gra *QuickSortFitness(struct node_gra *net, int nnodes){
    
    int i=1;
    struct node_gra *p=net;
    struct node_gra *list_net[max_size];
    
    //list_net=malloc( sizeof(struct node_gra *) * nlinks);
    
    //Passem la xarxa sencera a un vector
    list_net[0]=p;
    while(p->next!=NULL){ //go through adja list
        p=p->next;
        list_net[i]=p;
        i++;
    }
    
    
    quicksort_fitness_list(list_net,1,i-1);
    
    for(i=0;i<nnodes;i++){
        list_net[i]->next=list_net[i+1];
    }
    list_net[nnodes]->next=NULL;
    return list_net[0];
}

























//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static int partition_degree_list(struct node_gra *x[], int i, int j) {
    
    int k = j; 
    struct node_gra *v;
    struct node_gra *t;
    
    v=x[k];
    i--;
    while (i < j) {
        
        while (x[++i]->nlinks > v->nlinks);
        while (--j>i && x[j]->nlinks < v->nlinks);
        
        if (i < j) {
            t=x[i];	// swap x[i] and x[j] 
            x[i] = x[j];
            x[j] = t;
            //printf("---%d\n",x[j]->label);
        }
    }
    x[k] = x[i];	// swap x[i] and the pivot element 
    x[i] = v;
    return i;
}



//////////////////////////////////////////////////////////////////////
void quicksort_degree_list(struct node_gra *x[], int l, int r) {
    //use : quicksort(list_array,1,n) where n is the number of elements
    
    if (r > l) {
        int k = partition_degree_list(x, l, r);  
        quicksort_degree_list(x, l, k - 1);
        quicksort_degree_list(x, k + 1, r);
    }
}



///////////////////////////////////////////////////////////////////////
struct node_gra *QuickSortDegree(struct node_gra *net, int nnodes){
    
    int i=1;
    struct node_gra *p=net;
    struct node_gra *list_net[max_size];
    
    //list_net=malloc( sizeof(struct node_gra *) * nlinks);
    
    //Passem la xarxa sencera a un vector
    list_net[0]=p;
    while(p->next!=NULL){ //go through adja list
        p=p->next;
        p->nlinks=CountNodeLinks(p);		//Guardar tambe la informació de tots els graus
        list_net[i]=p;
        i++;
    }
    
    
    quicksort_degree_list(list_net,1,i-1);
    
    for(i=0;i<nnodes;i++){
        list_net[i]->next=list_net[i+1];
    }
    list_net[nnodes]->next=NULL;
    return list_net[0];
}

































//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static int partition_DegreeInv_list(struct node_gra *x[], int i, int j) {
    
    int k = j; 
    struct node_gra *v;
    struct node_gra *t;
    
    v=x[k];
    i--;
    while (i < j) {
        
        while (x[++i]->nlinks < v->nlinks);
        while (--j>i && x[j]->nlinks > v->nlinks);
        
        if (i < j) {
            t=x[i];	// swap x[i] and x[j] 
            x[i] = x[j];
            x[j] = t;
            //printf("---%d\n",x[j]->label);
        }
    }
    x[k] = x[i];	// swap x[i] and the pivot element 
    x[i] = v;
    return i;
}



//////////////////////////////////////////////////////////////////////
void quicksort_DegreeInv_list(struct node_gra *x[], int l, int r) {
    //use : quicksort(list_array,1,n) where n is the number of elements
    
    if (r > l) {
        int k = partition_DegreeInv_list(x, l, r);  
        quicksort_DegreeInv_list(x, l, k - 1);
        quicksort_DegreeInv_list(x, k + 1, r);
    }
}



///////////////////////////////////////////////////////////////////////
struct node_gra *QuickSortDegreeInv(struct node_gra *net, int nnodes){
    
    int i=1;
    struct node_gra *p=net;
    struct node_gra *list_net[max_size];
    
    //list_net=malloc( sizeof(struct node_gra *) * nlinks);
    
    //Passem la xarxa sencera a un vector
    list_net[0]=p;
    while(p->next!=NULL){ //go through adja list
        p=p->next;
        p->nlinks=CountNodeLinks(p);		//Guardar tambe la informació de tots els graus
        list_net[i]=p;
        i++;
    }
    
    
    quicksort_DegreeInv_list(list_net,1,i-1);
    
    for(i=0;i<nnodes;i++){
        list_net[i]->next=list_net[i+1];
    }
    list_net[nnodes]->next=NULL;
    return list_net[0];
}











//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static int partition_num_list(struct node_gra *x[], int i, int j) {
    
    int k = j; 
    struct node_gra *v;
    struct node_gra *t;
    
    v=x[k];
    i--;
    while (i < j) {
        
        while (x[++i]->num < v->num);
        while (--j>i && x[j]->num > v->num);
        
        if (i < j) {
            t=x[i];	// swap x[i] and x[j] 
            x[i] = x[j];
            x[j] = t;
            //printf("---%d\n",x[j]->label);
        }
    }
    x[k] = x[i];	// swap x[i] and the pivot element 
    x[i] = v;
    return i;
}



//////////////////////////////////////////////////////////////////////
void quicksort_num_list(struct node_gra *x[], int l, int r) {
    //use : quicksort(list_array,1,n) where n is the number of elements
    
    if (r > l) {
        int k = partition_num_list(x, l, r);  
        quicksort_num_list(x, l, k - 1);
        quicksort_num_list(x, k + 1, r);
    }
}



///////////////////////////////////////////////////////////////////////
struct node_gra *QuickSortNum(struct node_gra *net, int nnodes){
    
    int i=1;
    struct node_gra *p=net;
    struct node_gra *list_net[max_size];
    
    //list_net=malloc( sizeof(struct node_gra *) * nlinks);
    
    //Passem la xarxa sencera a un vector
    list_net[0]=p;
    while(p->next!=NULL){ //go through adja list
        p=p->next;
        list_net[i]=p;
        i++;
    }
    
    
    quicksort_num_list(list_net,1,i-1);
    
    for(i=0;i<nnodes;i++){
        list_net[i]->next=list_net[i+1];
    }
    list_net[nnodes]->next=NULL;
    return list_net[0];
}


