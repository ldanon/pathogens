#ifndef _NETWORK_LIB
#define _NETWORK_LIB
#include "globals.h"
#include <stdio.h>
#include <stdlib.h>


///////////////// GRAF ///////////////////////////////////////
struct node_gra{
    struct node_lis *neig;  // header of the adjacency list
    
    int num;                // label of the node
    char nom[20];           // name of the node
    double x,y,z;			// Posicio del node en el graf
    
    int state;              // to control wether visited or not
    int state_in;           // to control wether visited or not
    int state_out;          // to control wether visited or not
    
    double fitness;
   
    int nlinks;				// Degree
    int nlinks_ini;         // Degree of the initial graph
    int nlinks_weight;      // Weighted degree
	
	struct strain *infected, *immune, *immune2, *queue;		

    struct node_gra *next;  // next node in the graph
};


struct node_lis{
    int node;               // label of the node the link points to
    int status;             // for weighted networks  
    struct node_gra *ref;   // pointer to the node the link points to
    double btw;             // link betweenness (rush) variable.
    struct node_lis *link;  // other side of the link

    //Random walks			

    struct node_lis *next;  // next node in the adjacency list
};






//RUTINES DE LINKS
int AddAdjacency(struct node_gra *node,int label,int field);
int AddAdjacencyWeight(struct node_gra *node,int label);		//Si el link ja exiteix, inc+1 el pes dl link
void RewireAdjacency(struct node_gra *root);		//Crea tots els links de la llista d'adjacencies
void SortAndRewireAdjacency(struct node_gra *root);
void LinkSymmetricAdjacency(struct node_gra *net);
void PrintAdjacency(struct node_lis *p);
void RemoveAdjacencyList(struct node_lis *p);
void CopyAdjacencyList(struct node_gra *n1,struct node_gra *n2);
void ClearAdjacencies(struct node_gra *p,int net[max_size]);


int ExistLink(struct node_gra *a,struct node_gra *b);		//1 si existeix, 0 no
int CountNodeLinks(struct node_gra *node);		//# total links d'1 node
int CountNodeOutLinks(struct node_gra *node);	//# total outlinks d'1 node
int CountNodeInLinks(struct node_gra *node);	//# total inlinks d'1 node
int CountNodeLinksWeight(struct node_gra *node);		//# total links d'1 node weighted

int CountTotalLinks(struct node_gra *p);				//# total links de tota la xarxa
int CountTotalInLinks(struct node_gra *p);			//# total inlinks de tota la xarxa
int CountTotalOutLinks(struct node_gra *p);		//# total outlinks de tota la xarxa
int CountTotalLinksWeight(struct node_gra *p);				//# total links de tota la xarxa weighted

void RemoveLink(struct node_gra *n1,struct node_gra *n2);
void RemoveDirectedLink(struct node_gra *n1,struct node_gra *n2);
void CreateOneLink(struct node_gra *ori,struct node_gra *des);
struct node_gra *ChangeOneLink(struct node_gra *ori,struct node_gra *des,int N,struct node_gra *root);

void SimmetrizeWeights(struct node_gra *p1);

// RUTINES DE GRAFS
struct node_gra *CreateHeaderGraph();
struct node_gra *CreateNodeGraph(struct node_gra *network,int label,int x,int y);
struct node_gra *CopyNetwork(struct node_gra *p1);
void CopyNodeInfo(struct node_gra *dest, struct node_gra *orig);
struct node_gra *GetNode(int label,struct node_gra *p); // does what it says
void RemoveGraph(struct node_gra *p);

int Exist(int node,struct node_gra *p);		//comprova si el #node apareix dins el graf
int CountNodes(struct node_gra *p);			//Retorna num. total nodes del graf
void ResetNodesState(struct node_gra *p);   //Posa a state=0 x tots els nodes. 
void ResetLinksState(struct node_gra *p);   //Posa a state=0 x tots els nodes.  

void CreateFastAccess(struct node_gra *llista[], struct node_gra *p); // crea llista d'acces rapid a nodes segons ordre xarxa
struct node_gra *GetNodeFastAccess(int label,struct node_gra *llista[]); 
void CreateFastAccessNum(struct node_gra *llista[], struct node_gra *p); // crea llista d'acces rapid a nodes segons num. del node
struct node_gra *GetNodeFastAccessNum(int label,struct node_gra *llista[]); 
int CreateFastAccessWeight(struct node_gra *llista[], struct node_gra *p); // crea llista d'acces rapid a nodes segons el seu degree weig
struct node_gra *GetNodeFastAccessWeight(int label,struct node_gra *llista[]); 

//IMPRESSIO DE RESULTATS
void PrintGraph(struct node_gra *p);
void PrintPajekGraphTranslation(struct node_gra *root, char filename[15]);
void PrintPajekGraphTranslationNom(struct node_gra *root, char filename[15]);

void PrintPajekPartitionInfo(struct node_gra *root, char filename[15]);
void PrintPajekGraphTranslationWeighted(struct node_gra *root, char filename[15]);


//CREACIO DE DIFERENTS TIPUS DE GRAFS
struct node_gra *LoadPajekGraphTranslation(char filename[15]);
struct node_gra *LoadPajekGraphTranslationNom(char filename[15]);
struct node_gra *LoadPajekGraphTranslationXYZ(char filename[15]);

void LoadPartitionFile(struct node_gra *network, char filename[15]);
struct node_gra * LoadPajekLinkList(char filename[100],int mode);
// Creates a network from a file with format:          Node_1 Node_2 <----------------------
// Mode = 0 creates only the directed links in the list
// Mode = 1 creates the nodes in the list and symetrizes the adjacency matrix
// S = num_nodes of the network
struct node_gra * LoadPajekLinkListWithWeight(char filename[100],int mode);
// Creates a network from a file with format:          Node_1 Node_2 1 <----------------------

struct node_gra *BuildNetNoLinksMetapops(struct node_gra *nodes_array[]);

struct node_gra *BuildNetNoLinks(struct node_gra *nodes_array[]);

struct node_gra *CreateMStar(int S,int m);									// S: #nodes, m: #center nodes
struct node_gra *CreateCompleteGraph(int S);								// S: #nodes
struct node_gra *CreateSemiCompleteGraph(int S,int k);						// k:node degree - between s/2 and s-1
struct node_gra *CreateCommunityGraph(int N_NOD,int N_COMS, double Z_out, int AV_DEG);
struct node_gra *CreateCommunityGraph2(int N_NOD,int N_COMS, int COM_SIZE, double Z_out, double AV_DEG);
struct node_gra *CreateCommunityGraph3(int N_NOD,int N_COMS, int COM_SIZE, double Z_out, double AV_DEG);
struct node_gra *CreateCommunityGraphWeighted(int N_NOD,int N_COMS, int MAX_ZIN, int MAX);
struct node_gra *CreateCommunityGraphWeightedNewman(int N_NOD,int N_COMS, double Z_out, double AV_DEG, float W_I);
struct node_gra *CreateHierarchicalCommunityGraph(int N_NOD,int N_COMS1,int N_COMS2, int Z_in1, int Z_in2, int AV_DEG);
struct node_gra *CreateHexagonalLattice(int N_FILS,int N_COLS);
struct node_gra *CreateSquareLattice(int N_FILS,int N_COLS);

struct node_gra *CreateScaleFree(int S,int m);							// S: #nodes and m: minimum links
struct node_gra *CreateScaleFreeWithLocal(int S,int m,double p);
struct node_gra *CreateSimonScaleFree(int S,double a);	
struct node_gra *CreatePoissonRandomGraph(int S,int L);			// S: #nodes and L: #links
struct node_gra *CreatePoissonRandomGraphDirected(int S,int L);// S: #nodes and L: #links
struct node_gra *CreateHierarchicalNetwork(int m,int z);		// m: number of levels, z: branching factor
struct node_gra *Create1DSmallWorld(int n,int k,double p);  // n #nodes, k #links for each node, p: prob. link
struct node_gra *CreateInternetModel(int n,double r);

void PositionGraph(struct node_gra *net,int nnodes);	//Col.loca els nodes en un cercle de 0 a 1
void IniNodeNames(struct node_gra *net);	//Col.loca els nodes en un cercle de 0 a 1


#endif


