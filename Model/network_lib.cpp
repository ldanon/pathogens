#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>

#include <math.h>
#include <string.h>
#include <time.h>



#include "network_lib.h"
#include "quicksort.h"
#include "rand_gen.h"
#include "globals.h"
#include "pathogens.h"
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// BUILD THE ELEMENTS OF THE NETWORK




int AddAdjacency(struct node_gra *node,int label,int field)
{
	// Returns 0 if no link has been created (because it already existed) 
	// or one otherwise

	struct node_lis *adja;

	if(node->num==label) // dont make link with itself
		return 0;
	else{
		adja=node->neig;
		while(adja->next!=NULL){
			adja=adja->next;
			if(adja->node==label) // Dont make links which already exist
				return 0;
		}

		adja->next=(struct node_lis *)malloc(sizeof(struct node_lis));
		(adja->next)->node=label;
		(adja->next)->status=field;
		(adja->next)->next=NULL;
		(adja->next)->ref=NULL;

		return 1;
	}
}



///////////////////////////////////////////////////////////////////////
int AddAdjacencyWeight(struct node_gra *node,int label)
{
	struct node_lis *adja;

	if(node->num==label) // Dont make a link with itself NO LOOPS!
		return 0;
	else{
		adja=node->neig;         //look through link list till 
		while(adja->next!=NULL){ //you get to end
			adja=adja->next;       
			if(adja->node==label){ // Dont place links which already exist
				adja->status++;      // instead increment by one the weight of the link
				/*  	printf("cstatus=%d\n",adja->status); */
				return 0;
			}
		}

		adja->next=(struct node_lis *)
			malloc(sizeof(struct node_lis)); //allocate memory for _NEXT_ in list

		(adja->next)->node=label; // label of node it points TO
		(adja->next)->status=1;   
		(adja->next)->next=NULL;  //set next to NULL, to stop from seg_faulting
		(adja->next)->ref=NULL;

		return 1;
	}
}


///////////////////////////////////////////////////////////////////////
void RewireAdjacency(struct node_gra *root)
{
	// sets up the adjacency list of the whole network again
	int i;
	struct node_gra *p=root;
	struct node_lis *adja;
	struct node_gra *vector[max_size];
	int counter=0;

	for(i=0;i<max_size;i++){ // initialise vector
		vector[i]=NULL;
	}

	while(p->next!=NULL){   // copy network to vector
		p=p->next;
		vector[p->num]=p;
	}

	p=root;

	while(p->next!=NULL){  // go through network
		p=p->next;           
		p->nlinks=0;
		counter++;           // count nodes
		adja=p->neig;        // get adjacency list of each node
		while(adja->next!=NULL){ //go through adja list
			adja=adja->next;       
			adja->ref=vector[adja->node]; // make ref pointer point to node with label node
			p->nlinks++;
		}
	}
}



///////////////////////////////////////////////////////////////////////
void SortAndRewireAdjacency(struct node_gra *root)
{
	// sets up the adjacency list of the whole network again
	// and sort the link lists so that they are in order

	int i;
	struct node_gra *p=root;
	struct node_lis *adja, *list_ra[max_size];

	struct node_gra *vector[max_size];
	int counter=0;

	for(i=0;i<max_size;i++){ // initialise vector
		vector[i]=NULL;
	}

	while(p->next!=NULL){   // copy network to vector
		p=p->next;
		vector[p->num]=p;
	}

	p=root;

	while(p->next!=NULL){  // go through network
		p=p->next;           
		counter++;           // count nodes

		adja=p->neig;        // get adjacency list of each node
		i=1;
		while(adja->next!=NULL){ //go through adja list
			adja=adja->next;       
			list_ra[i]=adja;
			i++;
			adja->ref=vector[adja->node]; // make ref pointer point to node with label node
		}
		(p->neig)->next=QuickSortLinkList(list_ra,i-1);	//Sort the link list
	}
}




void LinkSymmetricAdjacency(struct node_gra *net)
{
	struct node_gra *p=net,*p2;
	struct node_lis *l,*l2;

	//Firstly, initiallize all links
	while(p->next!=NULL) {
		p=p->next;
		l=p->neig;
		while(l->next!=NULL)
		{
			l=l->next;
			l->link=NULL;
		}
	}

	p=net;
	while(p->next!=NULL) {
		p=p->next;
		l=p->neig;
		while(l->next!=NULL)
		{
			l=l->next;
			if(l->link==NULL)
			{
				//The we link the two sides
				p2=l->ref;
				l2=p2->neig;
				while(l2->next!=NULL)
				{
					l2=l2->next;
					if(l2->ref->num == p->num) {
						l->link=l2;
						l2->link=l;
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////
void PrintAdjacency(struct node_lis *p)
{
	while(p->next!=NULL){
		printf(" %d #",(p->next)->node);
		p=p->next;
	}
}


///////////////////////////////////////////////////////////////////////
void RemoveAdjacencyList(struct node_lis *p)
{
	if(p==NULL) return;

	if(p->next!=NULL)     //see above
		RemoveAdjacencyList(p->next);

	free(p);
	p=NULL;
}


//////////////////////////////////////////////////////////////////////
void CopyAdjacencyList(struct node_gra *n1,struct node_gra *n2)
{
	struct node_lis *p1;

	p1=n1->neig;

	while(p1->next!=NULL){
		p1=p1->next;
		AddAdjacency(n2,p1->node,p1->status);
	}
}



////////////////////////////////////////////////////////////////////////
void ClearAdjacencies(struct node_gra *p,int net[max_size])
{
	//Elimina els links que van dirigits cap als nodes que tenen un 0
	//en la llista net[numnode]
	struct node_lis *nei,*temp;

	while(p->next!=NULL){
		p=p->next;
		nei=p->neig;
		while(nei->next!=NULL){
			if(net[(nei->next)->node]==0){
				temp=nei->next;
				if(temp->next==NULL)  nei->next=NULL;
				else  nei->next=temp->next;
				free(temp);
			}else{
				nei=nei->next;
			}
		}
	}

	return;
}


//////////////////////////////////////////////////////////////////////
// does a link exist between node a and b, return 1 if yes 0 if no
int ExistLink(struct node_gra *a,struct node_gra *b)
{
	struct node_lis *an=a->neig;

	while(an->next!=NULL){
		an=an->next;
		if(an->node==b->num){
			return 1;
		}
	}

	return 0;
}



//////////////////////////////////////////////////////////////////////
int CountNodeLinks(struct node_gra *node)
{
	struct node_lis *p=node->neig;
	int count=0;

	while(p->next!=NULL){
		p=p->next;
		count++;
	}
	node->nlinks=count;

	return count;
}


//////////////////////////////////////////////////////////////////////
int CountNodeOutLinks(struct node_gra *node)   //OUT=Status a 1
{
	struct node_lis *p=node->neig;
	int counter=0;

	while(p->next!=NULL){
		p=p->next;
		if(p->status==1)
			counter++;
	}

	return counter;
}


//////////////////////////////////////////////////////////////////////
int CountNodeInLinks(struct node_gra *node)   //IN=Status a 0
{
	struct node_lis *p=node->neig;
	int counter=0;

	while(p->next!=NULL){
		p=p->next;
		if(p->status==0)
			counter++;
	}

	return counter;
}


//////////////////////////////////////////////////////////////////////
int CountNodeLinksWeight(struct node_gra *node)
{
	struct node_lis *p=node->neig;
	int count=0;

	while(p->next!=NULL){
		p=p->next;
		count+=p->status;
	}
	node->nlinks_weight=count;

	return count;
}



///////////////////////////////////////////////////////////////////////
int CountTotalLinks(struct node_gra *p)
{
	int total=0;

	while(p->next!=NULL){
		p=p->next;
		p->nlinks=CountNodeLinks(p);
		total+=p->nlinks;

	}

	return total/2;
}


///////////////////////////////////////////////////////////////////////
int CountTotalInLinks(struct node_gra *p)
{
	int total=0;

	while(p->next!=NULL){
		p=p->next;
		total+=CountNodeInLinks(p);
	}

	return total/2;
}


///////////////////////////////////////////////////////////////////////
int CountTotalOutLinks(struct node_gra *p)
{
	int total=0;

	while(p->next!=NULL){
		p=p->next;
		total+=CountNodeOutLinks(p);
	}

	return total/2;
}


///////////////////////////////////////////////////////////////////////
int CountTotalLinksWeight(struct node_gra *p)
{
	int total=0;

	while(p->next!=NULL){
		p=p->next;
		total+=CountNodeLinksWeight(p);
	}

	return total/2;
}







///////////////////////////////////////////////////////////////////////
void RemoveLink(struct node_gra *n1,struct node_gra *n2)
{
	//Borrem el link en els 2 nodes (Unidrected!!)
	struct node_lis *nn1, *nn2;
	struct node_lis *temp1, *temp2;

	nn1=n1->neig;
	nn2=n2->neig;

	while((nn1->next)->ref!=n2){
		nn1=nn1->next;
	}
	temp1=nn1->next;
	nn1->next=temp1->next;
	free(temp1);

	while((nn2->next)->ref!=n1){
		nn2=nn2->next;
	}
	temp2=nn2->next;
	nn2->next=temp2->next;
	free(temp2);
}


///////////////////////////////////////////////////////////////////////
void RemoveDirectedLink(struct node_gra *n1,struct node_gra *n2)
{
	struct node_lis *nn1;
	struct node_lis *temp1;

	nn1=n1->neig;

	while((nn1->next)->ref!=n2){
		nn1=nn1->next;
	}

	temp1=nn1->next;
	nn1->next=temp1->next;
	free(temp1);
}


///////////////////////////////////////////////////////////////////////
void CreateOneLink(struct node_gra *ori,struct node_gra *des)
{
	struct node_lis *ne;
	int exist;

	exist=AddAdjacency(ori,des->num,0);
	exist=AddAdjacency(des,ori->num,0);

	// Create the physical link
	ne=ori->neig;
	do{
		ne=ne->next;
	}while(ne->node!=des->num);
	ne->ref=des;

	ne=des->neig;
	do{
		ne=ne->next;
	}while(ne->node!=ori->num);
	ne->ref=ori;
}


///////////////////////////////////////////////////////////////////////
struct node_gra *ChangeOneLink(struct node_gra *ori,struct node_gra *des,
	int N,struct node_gra *root)
{
	//N: nnodes total de la xarxa
	//Estreu el link i el col.loca aleatoriament entre el node origen i un nou node desti 
	//seleccionat random
	int targ;
	struct node_lis *ne;
	struct node_gra *des2;
	int exist;
	long seed;

	seed=-time(NULL);

	RemoveLink(ori,des);

	do{
		do{
			targ=(int) floor(ran2(&seed)*(double)N);
		}while((targ==ori->num)||(targ==des->num));
		des2=GetNode(targ,root);
		exist=AddAdjacency(ori,des2->num,0);
		exist=AddAdjacency(des2,ori->num,0);
	}while(exist==0);

	// Create the physical link
	ne=ori->neig;
	do{
		ne=ne->next;
	}while(ne->node!=des2->num);
	ne->ref=des2;
	ne=des2->neig;
	do{
		ne=ne->next;
	}while(ne->node!=ori->num);
	ne->ref=ori;

	return des2;
}


void SimmetrizeWeights(struct node_gra *p1)
{
	struct node_lis *l1,*l2;
	struct node_gra *p2;

	while(p1->next!=NULL){
		p1=p1->next;
		l1=p1->neig;
		while(l1->next!=NULL){
			l1=l1->next;
			p2=l1->ref;
			l2=(p2->neig)->next;
			while(l2->ref!=p1){
				l2=l2->next;
			}
			if(l1->status<l2->status)
				l2->status=l1->status;
		}
	}
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// BUILD THE ELEMENTS OF THE NETWORK


struct node_gra *CreateHeaderGraph()
{
	struct node_gra *temp;

	temp=(struct node_gra *)malloc(sizeof(struct node_gra));
	temp->num=-1;
	temp->neig=NULL;
	temp->next=NULL;
	temp->state_in=0;
	temp->state_out=0;

	return temp;
}



///////////////////////////////////////////////////////////////////////
struct node_gra *CreateNodeGraph(struct node_gra *network, int label,int x,int y)	//x i y no serveixen
{
	x=y;
	while(network->next!=NULL)    // find last node in network
		network=network->next;      // THEN add another node (network->next)

	network->next=(struct node_gra *)malloc(sizeof(struct node_gra));
	(network->next)->num=label;
	(network->next)->state=0;
	(network->next)->state_in=0;
	(network->next)->state_out=0;
	(network->next)->next=NULL;
	(network->next)->neig=(struct node_lis *)malloc(sizeof(struct node_lis));
	((network->next)->neig)->node=-1;
	((network->next)->neig)->next=NULL;

	return network->next;         
}


///////////////////////////////////////////////////////////////////////
void RemoveGraph(struct node_gra *p)
{
	struct node_gra *p2;

	if(p==NULL) return;

	while(p->next!=NULL)
	{
		p2=p;
		p=p->next;
		RemoveAdjacencyList(p2->neig);
		free(p2);
	}

	RemoveAdjacencyList(p->neig);
	free(p);
	p=NULL;
}



///////////////////////////////////////////////////////////////////////
struct node_gra *CopyNetwork(struct node_gra *p1)
{
	struct node_gra *root2,*p2;
	struct node_gra *last;

	root2=CreateHeaderGraph();
	last=root2;

	while(p1->next!=NULL){
		p1=p1->next;
		p2=CreateNodeGraph(last,p1->num,0,0);
		p2->fitness=p1->fitness;
		p2->state=p1->state;
		p2->state_in=p1->state_in;
		p2->state_out=p1->state_out;
		p2->x=p1->x;p2->y=p1->y;p2->z=p1->z;
		p2->nlinks_ini=p1->nlinks_ini;
		p2->nlinks=p1->nlinks;
		strcpy(p2->nom,p1->nom);
		last=p2;
		CopyAdjacencyList(p1,p2);
	}

	RewireAdjacency(root2);

	return root2;
}



///////////////////////////////////////////////////////////////////////
void CopyNodeInfo(struct node_gra *dest, struct node_gra *orig)
{
	dest->fitness=orig->fitness;
	dest->state=orig->state;
	dest->state_in=orig->state_in;
	dest->state_out=orig->state_out;
	dest->x=orig->x;dest->y=orig->y;dest->z=orig->z;
	dest->nlinks_ini=orig->nlinks_ini;
	dest->nlinks=orig->nlinks;
	strcpy(dest->nom,orig->nom);
}

///////////////////////////////////////////////////////////////////////
struct node_gra *GetNode(int label,struct node_gra *net) // does what it says
{              
	struct node_gra *p=net;

	while((p->next)->num!=label)
		p=p->next;

	return p->next;
}







///////////////////////////////////////////////////////////////////////
int Exist(int node,struct node_gra *p)
{

	while(p->next!=NULL){
		p=p->next;
		if(p->num==node)
			return 1;
	}

	return 0;
}


///////////////////////////////////////////////////////////////////////
int CountNodes(struct node_gra *p)
{
	int nodes=0;

	while(p->next!=NULL){
		nodes++;
		p=p->next;
	};

	return nodes;
}


///////////////////////////////////////////////////////////////////////
void ResetNodesState(struct node_gra *p) // reset node_gra state
{
	while(p->next!=NULL){
		p=p->next;
		p->state=0;
		p->state_in=0;
		p->state_out=0;
		p->fitness=0;
		//        p->nlinks=CountNodeLinks(p);		//Guardar tambe la informaci� de tots els graus
	}
}












///////////////////////////////////////////////////////////////////////
void CreateFastAccess(struct node_gra *llista[], struct node_gra *p)// crea llista d'acces rapid a nodes
{
	int i;

	//Funcio que proporciona acces rapid als nodes de la llista,
	//Es util sols en el cas de que tinguem els nodes de 0 al final ordenats
	i=0;
	while(p->next!=NULL){
		llista[i]=p;
		p=p->next;
		i++;
	}
}

///////////////////////////////////////////////////////////////////////
struct node_gra *GetNodeFastAccess(int label,struct node_gra *llista[])
{
	return llista[label];
}



///////////////////////////////////////////////////////////////////////
void CreateFastAccessNum(struct node_gra *llista[], struct node_gra *p)// crea llista d'acces rapid a nodes
{
	int i,nnod;

	nnod=CountNodes(p);

	for(i=0;i<nnod;i++) llista[i]=NULL;

	//Funcio que proporciona acces rapid als nodes de la llista,
	while(p->next!=NULL){
		p=p->next;
		llista[p->num]=p;
	}
}

///////////////////////////////////////////////////////////////////////
struct node_gra *GetNodeFastAccessNum(int label,struct node_gra *llista[])
{
	return llista[label];
}


///////////////////////////////////////////////////////////////////////
int CreateFastAccessWeight(struct node_gra *llista[], struct node_gra *p)// crea llista d'acces rapid a nodes
{
	int i=0,j;

	//Funcio que proporciona acces rapid als nodes de la llista en una proporcio depenent al seu degree weight, ENTRE 0 I NUM_LINKS_W-1
	while(p->next!=NULL){
		p=p->next;
		for(j=0;j<p->nlinks_weight;j++){
			llista[i]=p;
			i++;
		}
	}
	return i;//Num posicions creades dins la taula
}

///////////////////////////////////////////////////////////////////////
struct node_gra *GetNodeFastAccessWeight(int label,struct node_gra *llista[])
{
	return llista[label];
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void PrintGraph(struct node_gra *p)
{
	while(p->next!=NULL){
		printf(" %d --- ",(p->next)->num);
		PrintAdjacency((p->next)->neig);
		printf("\n");
		p=p->next;
	}
}



///////////////////////////////////////////////////////////////////////
void PrintPajekGraphTranslation(struct node_gra *root, char filename[15])
{
	struct node_gra *p=root;
	struct node_lis *li;
	FILE *fit;
	int trans[max_size];
	int i,co;

	//TANIS char command[30]="unix2dos ";
	//TANIS strcat(command,filename);

	for(i=0;i<max_size;i++)
		trans[i]=0;

	fit=fopen(filename,"w");

	fprintf(fit,"*Vertices %d\n",CountNodes(root));

	co=1;
	p=root;
	while(p->next!=NULL){
		p=p->next;
		if(trans[p->num]==0){
			trans[p->num]=co;
			co++;
		}
		fprintf(fit,"%d *%d*\n",trans[p->num],p->num+1);
	}

	fprintf(fit,"*Arcs\n");
	fprintf(fit,"*Edges\n");
	p=root;
	while(p->next!=NULL){
		p=p->next;
		li=p->neig;
		while(li->next!=NULL){
			li=li->next;
			/*        if(p->num<li->node) */
			fprintf(fit,"%d   %d   1\n",trans[p->num],trans[li->node]);
		}
	}

	fclose(fit);
}






///////////////////////////////////////////////////////////////////////
void PrintPajekGraphTranslationNom(struct node_gra *root, char filename[15])
{
	struct node_gra *p=root;
	struct node_lis *li;
	FILE *fit;
	int trans[max_size];
	int i,co;

	//TANIS char command[30]="unix2dos ";
	//TANIS strcat(command,filename);

	for(i=0;i<max_size;i++)
		trans[i]=0;

	fit=fopen(filename,"w");

	fprintf(fit,"*Vertices %d\n",CountNodes(root));

	co=1;
	p=root;
	while(p->next!=NULL){
		p=p->next;
		if(trans[p->num]==0){
			trans[p->num]=co;
			co++;
		}
		fprintf(fit,"%d %s   %.6f  %.6f  %.6f\n",trans[p->num],p->nom,p->x,p->y,p->z);
	}

	fprintf(fit,"*Arcs\n");
	fprintf(fit,"*Edges\n");
	p=root;
	while(p->next!=NULL){
		p=p->next;
		li=p->neig;
		while(li->next!=NULL){
			li=li->next;
			/*        if(p->num<li->node) */
			fprintf(fit,"%d   %d   1\n",trans[p->num],trans[li->node]);
		}
	}

	fclose(fit);
}







///////////////////////////////////////////////////////////////////////
void PrintPajekLinkListTranslation(struct node_gra *root, char filename[15])
{
	struct node_gra *p=root;
	struct node_lis *li;
	FILE *fit;
	int trans[max_size];
	int i,co;

	//TANIS char command[30]="unix2dos ";
	//TANIS strcat(command,filename);

	for(i=0;i<max_size;i++)
		trans[i]=0;

	fit=fopen(filename,"w");

	co=1;
	p=root;
	while(p->next!=NULL){
		p=p->next;
		if(trans[p->num]==0){
			trans[p->num]=co;
			co++;
		}
	}

	//    fprintf(fit,"*Arcs\n");
	//    fprintf(fit,"*Edges\n");
	p=root;
	while(p->next!=NULL){
		p=p->next;
		li=p->neig;
		while(li->next!=NULL){
			li=li->next;
			/*        if(p->num<li->node) */
			fprintf(fit,"%d   %d\n",trans[p->num],trans[li->node]);
		}
	}

	fclose(fit);
}



///////////////////////////////////////////////////////////////////////
void PrintPajekPartitionInfo(struct node_gra *root, char filename[15])
{
	struct node_gra *p=root;
	FILE *fit;

	fit=fopen(filename,"w");
	fprintf(fit,"*Vertices %d\n",CountNodes(root));

	QuickSortNum(p, CountNodes(root));
	while(p->next!=NULL){
		p=p->next;
		fprintf(fit,"%d\n",p->state+1);
	}

	fclose(fit);
}



///////////////////////////////////////////////////////////////////////
void PrintPajekGraphTranslationWeighted(struct node_gra *root, char filename[15])
{
	struct node_gra *p=root;
	struct node_lis *li;
	FILE *fit;
	int trans[max_size];
	int i,co;

	//TANIS char command[30]="unix2dos ";
	//TANIS strcat(command,filename);

	for(i=0;i<max_size;i++)
		trans[i]=0;

	fit=fopen(filename,"w");

	fprintf(fit,"*Vertices %d\n",CountNodes(root));

	co=1;
	p=root;
	while(p->next!=NULL){
		p=p->next;
		if(trans[p->num]==0){
			trans[p->num]=co;
			co++;
		}
		fprintf(fit,"%d %s   %.6f  %.6f  %.6f\n",trans[p->num],p->nom,p->x,p->y,p->z);
	}

	fprintf(fit,"*Arcs\n");
	fprintf(fit,"*Edges\n");
	p=root;
	while(p->next!=NULL){
		p=p->next;
		li=p->neig;
		while(li->next!=NULL){
			li=li->next;
			/*        if(p->num<li->node) */
			fprintf(fit,"%d   %d   %d\n",trans[p->num],trans[li->node],li->status);
		}
	}

	fclose(fit);
}






///////////////////////////////////////////////////////////////////////
struct node_gra *LoadPajekGraphTranslation(char filename[15]){
	FILE *inFile;
	struct node_gra *root=NULL;
	struct node_gra *list[max_size];
	int i,t1,t2,t3;
	int nnod;
	char stmp[10],line[255];
	struct node_gra *last_add=NULL;

	inFile=fopen(filename,"r");

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	last_add=root=CreateHeaderGraph();

	//Llegim el num. de nodes
	fgets(line,255,inFile);
	sscanf(line,"%s %d\n",stmp,&nnod);


	//Llegim els nodes de la llista
	for(i=0;i<nnod;i++){
		fgets(line,255,inFile);
	}


	//Llegim el num. de vertexs
	fgets(line,255,inFile);
	//Llegim el num. de arcs
	fgets(line,255,inFile);


	// Read the data and count nodes
	nnod=0;
	while(!feof(inFile)){
		fgets(line,255,inFile);
		sscanf(line,"%d %d %d\n",&t1,&t2,&t3);
		//printf("%d %d\n",t1,t2);

		if(list[t1]==NULL){
			list[t1]=CreateNodeGraph(last_add,t1,t1,t1);
			last_add=list[t1];
			nnod++;
		}
		if(list[t2]==NULL){
			list[t2]=CreateNodeGraph(last_add,t2,t2,t2);
			last_add=list[t2];
			nnod++;
		}
		AddAdjacency(list[t2],t1,t3);
		AddAdjacency(list[t1],t2,t3);
	}

	//printf("           Number of nodes: %d\n",nnod);

	SortAndRewireAdjacency(root);

	fclose(inFile);

	return root;
}








struct node_gra *BuildNetNoLinks(struct node_gra *nodes_array[]){

	int i;
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<N_NODE;i++) nodes_array[i]=NULL;

	nodes_array[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<N_NODE;i++){
		nodes_array[i]=CreateNodeGraph(root,i,i,i);
	}

	return root;
}
	
struct node_gra *BuildNetNoLinksMetapops(struct node_gra *nodes_array[]){

	int i;
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<N_NODE;i++) nodes_array[i]=NULL;

	nodes_array[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<N_NODE;i++){
		nodes_array[i]=CreateNodeGraph(nodes_array[i-1],i,0,0);
		nodes_array[i]->state=i%NUM_COM;
	}

	return root;
}




struct node_gra *LoadPajekGraphTranslationXYZ(char filename[15]){
	FILE *inFile;
	struct node_gra *root=NULL,*p;
	struct node_gra *list[max_size];
	float x,y,z;
	//    char noms[1000][20];
	int i,t1,t2,t3;
	int nnod;
	char stmp[30],line[255];
	struct node_gra *last_add=NULL;

	inFile=fopen(filename,"r");

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	last_add=root=CreateHeaderGraph();

	//Llegim el num. de nodes
	fgets(line,255,inFile);
	sscanf(line,"%s %d\n",stmp,&nnod);


	//Llegim els nodes de la llista
	for(i=0;i<nnod;i++){
		fgets(line,255,inFile);
		sscanf(line,"%d %s %f %f %f\n",&t1,stmp,&x,&y,&z);
	}


	//Llegim el num. de vertexs
	fgets(line,255,inFile);
	//Llegim el num. de arcs
	fgets(line,255,inFile);


	// Read the data and count nodes
	nnod=0;
	while(!feof(inFile)){
		fgets(line,255,inFile);
		sscanf(line,"%d %d %d\n",&t1,&t2,&t3);
		//printf("%d %d\n",t1,t2);

		if(list[t1]==NULL){
			list[t1]=CreateNodeGraph(last_add,t1,t1,t1);
			last_add=list[t1];
			nnod++;
		}
		if(list[t2]==NULL){
			list[t2]=CreateNodeGraph(last_add,t2,t2,t2);
			last_add=list[t2];
			nnod++;
		}
		AddAdjacency(list[t2],t1,t3);
		AddAdjacency(list[t1],t2,t3);
	}
	fclose(inFile);   

	SortAndRewireAdjacency(root);

	QuickSortNum(root,CountNodes(root));

	//Llegim dades addicionals del fitxer .net
	inFile=fopen(filename,"r");
	p=root;
	fgets(line,255,inFile);
	for(i=0;i<nnod;i++){
		fgets(line,255,inFile);
		sscanf(line,"%d %s %f %f %f\n",&t1,stmp,&x,&y,&z);
		p=p->next;
		p->x=x;
		p->y=y;
		p->z=z;
		strcpy(p->nom,stmp);
	}
	fclose(inFile);   

	return root;
}



struct node_gra * LoadPajekLinkList(char filename[100],int mode)
{
	FILE *inFile;
	struct node_gra *root=NULL;
	struct node_gra *list[max_size];
	struct node_gra *last_add=NULL;

	int i,t1,t2;
	int nnod;

	inFile=fopen(filename,"r");

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	last_add=root=CreateHeaderGraph();

	// Read the data and count nodes
	nnod=0;
	while(!feof(inFile)){
		fscanf(inFile,"%d %d\n",&t1,&t2);
		if(t1==0||t2==0){
			printf("ERROR: Zero in link list\n node1:%d node2:%d\n",t1,t2);
			printf("Use BuildNetworkFormatArxiv instead!!!\n");
			return NULL;
		}
		if(list[t1]==NULL){
			list[t1]=CreateNodeGraph(last_add,t1,t1,t1);
			last_add=list[t1];
			nnod++;
		}
		if(list[t2]==NULL){
			list[t2]=CreateNodeGraph(last_add,t2,t2,t2);
			last_add=list[t2];
			nnod++;
		}
		AddAdjacency(list[t1],t2,1);
		if(mode==1){
			AddAdjacency(list[t2],t1,1);
		}
		else if(mode!=0){
			printf("Error in BuildNetworkFormatTree mode: must be 0 or 1!!\n");
			return NULL;
		}
	}

	//    printf("           Number of nodes: %d\n",nnod);

	fclose(inFile);

	SortAndRewireAdjacency(root);

	return root;
}





struct node_gra * LoadPajekLinkListWithWeight(char filename[100],int mode)
{
	FILE *inFile;
	struct node_gra *root=NULL;
	struct node_gra *list[max_size];
	struct node_gra *last_add=NULL;
	int i,t1,t2;
	int nnod;

	inFile=fopen(filename,"r");

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	last_add=root=CreateHeaderGraph();

	// Read the data and count nodes
	nnod=0;
	while(!feof(inFile)){
		fscanf(inFile,"%d %d 1\n",&t1,&t2);
		if(t1==0||t2==0){
			printf("ERROR: Zero in link list\n node1:%d node2:%d\n",t1,t2);
			printf("Use BuildNetworkFormatArxiv instead!!!\n");
			return NULL;
		}
		if(list[t1]==NULL){
			list[t1]=CreateNodeGraph(last_add,t1,t1,t1);
			last_add=list[t1];
			nnod++;
		}
		if(list[t2]==NULL){
			list[t2]=CreateNodeGraph(last_add,t2,t2,t2);
			last_add=list[t2];
			nnod++;
		}
		AddAdjacency(list[t1],t2,1);
		if(mode==1){
			AddAdjacency(list[t2],t1,1);
		}
		else if(mode!=0){
			printf("Error in BuildNetworkFormatTree mode: must be 0 or 1!!\n");
			return NULL;
		}
	}

	//    printf("           Number of nodes: %d\n",nnod);

	fclose(inFile);

	SortAndRewireAdjacency(root);

	return root;
}



//Carrega un fitxer.clu i assigna als nodes la particio on 
//corresponen en la variable state
void LoadPartitionFile(struct node_gra *network, char filename[15]){
	FILE *inFile;
	struct node_gra *p;
	int t1,t3,i;
	char ts1[20];

	inFile=fopen(filename,"r");

	network=QuickSortNum(network,CountNodes(network));
	p=network;

	fscanf(inFile,"%s %i\n",ts1,&t3);
	for(i=0;i<t3;i++)
	{
		p=p->next;
		fscanf(inFile," %i",&t1);						//Agafem el valor de la particio on correspon
		p->state=t1;										//I li assignem el valor llegit	
	}

	fclose(inFile);
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateMStar(int S,int m)
{
	int i,j;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	list[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<S;i++){
		list[i]=CreateNodeGraph(list[i-1],i,i,0);
	}

	for(i=m;i<S;i++){
		for(j=0;j<m;j++){
			AddAdjacency(list[i],j,0);
			AddAdjacency(list[j],i,0);
		}
	}

	RewireAdjacency(root);

	return root;
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateCompleteGraph(int S)
{
	int i,j;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	list[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<S;i++){
		list[i]=CreateNodeGraph(list[i-1],i,i,i);
	}

	for(i=0;i<S;i++){
		for(j=i+1;j<S;j++){
			AddAdjacency(list[i],j,0);
			AddAdjacency(list[j],i,0);
		}
	}

	RewireAdjacency(root);

	return root;
}








///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateSemiCompleteGraph(int S,int k)   //k:node degree - between s/2 and s-1
{
	int i,j,pos;
	int degree_out,degree=k;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	if( degree < (S/2) ) degree = (S/2);
	if( degree > (S-1) ) degree = (S-1);
	degree_out=degree-(S/2)+1;

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	list[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<S;i++){
		list[i]=CreateNodeGraph(list[i-1],i,i,i);
	}


	for(i=0;i<S/2;i++){
		for(j=i+1;j<S/2;j++){
			AddAdjacency(list[i],j,0);
			AddAdjacency(list[j],i,0);

			AddAdjacency(list[i+S/2],j+S/2,0);
			AddAdjacency(list[j+S/2],i+S/2,0);
		}
	}

	for(i=0;i<S/2;i++){
		for(j=0;j<degree_out;j++){
			pos=(i+j)%(S/2) + S/2;
			AddAdjacency(list[i],pos,0);
			AddAdjacency(list[pos],i,0);
		}
	}

	RewireAdjacency(root);

	return root;
}






///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//int N_NOD,int N_COMS1,int N_COMS2, double Z_in1, double Z_in2, int AV_DEG
//z_in1=Num. de links cap a nodes de la comunitat superior
//z_in2=Num. de links cap a nodes de la comunitat meva
struct node_gra *CreateCommunityGraph(int N_NOD,int N_COMS, double Z_out, int AV_DEG)
{
	int COMS[max_size];
	int i,j;
	double Z_in=AV_DEG-Z_out;
	int COM_SIZE = N_NOD / N_COMS;
	double P_in,P_out;

	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++) list[i]=NULL;

	list[0]=CreateNodeGraph(root,1,1,1);
	for(i=1;i<N_NOD;i++){
		list[i]=CreateNodeGraph(list[i-1],i+1,i+1,i+1);
	}


	P_in=Z_in/(double)(COM_SIZE-1);
	P_out=Z_out/(double)(N_NOD-COM_SIZE);


	for(i=0;i<N_NOD;i++) {
		list[i]->state=i%N_COMS;
		COMS[i]=i%N_COMS;
	}


	for(i=0;i<N_NOD;i++){
		for(j=i+1;j<N_NOD;j++){
			if(COMS[i]==COMS[j]){
				if(ran2(&seed)<P_in){
					AddAdjacency(list[i],j+1,0);
					AddAdjacency(list[j],i+1,0);
				}
			}else {
				if(ran2(&seed)<P_out){
					AddAdjacency(list[i],j+1,0);
					AddAdjacency(list[j],i+1,0);
				}
			}
		}
	}

	RewireAdjacency(root);

	return root;
}





struct node_gra *CreateCommunityGraph3(int N_NOD,int N_COMS, int COM_SIZE,
				      double Z_out, double AV_DEG)
{
  
				/* N_NOD = number of nodes in net */
				/* N_COMS = number of communities */
 				/* COM_SIZE = size of community */
  				/* Z_out = av number of links pointing out */
  				/* AV_DEG = average degree of all nodes */
				
  int *COMS;
  double Z_in;
  double P_in,P_out;
  int i,j;
 
  struct node_gra **list;
  struct node_gra *root;
 
  COMS=(int *)calloc(N_NOD+1,sizeof(int)); 

  list=(struct node_gra **)calloc(N_NOD+1,sizeof(struct node_gra *));
  
  root=CreateHeaderGraph();
 
  for(i=0;i<N_NOD+1;i++){
    list[i]=NULL;
  }
  
  list[0]=CreateNodeGraph(root,0,0,0);
  for(i=1;i<N_NOD;i++){
    list[i]=CreateNodeGraph(list[i-1],i,i,i);
  }
  
  
  Z_in=AV_DEG-Z_out;
  
  P_in=(double)Z_in/(COM_SIZE-1);
  P_out=(double)Z_out/(N_NOD-COM_SIZE);
  
  
  for(i=0;i<N_NOD;i++){
    COMS[i]=i%N_COMS;
  }
  
  for(i=0;i<N_NOD;i++){
    for(j=i+1;j<N_NOD;j++){
      if(COMS[i]==COMS[j]){
	if(ran2(&seed)<P_in){
	  
	  AddAdjacency(list[i],j,0);
	  
	  AddAdjacency(list[j],i,0);
	}
      }
      else {
	if(ran2(&seed)<P_out){
	  
	  AddAdjacency(list[i],j,0);
	  
	  AddAdjacency(list[j],i,0);
	}
      }
    }
  }
  
  
  RewireAdjacency(root);
  
  return root;
}

struct node_gra *CreateCommunityGraph2(int N_NOD,int N_COMS, int COM_SIZE,
				      double Z_out, double AV_DEG)
{
  
				/* N_NOD = number of nodes in net */
				/* N_COMS = number of communities */
 				/* COM_SIZE = size of community */
  				/* Z_out = av number of links pointing out */
  				/* AV_DEG = average degree of all nodes */
				
  int *COMS;
  double Z_in;
  double P_in,P_out;
  int i,j;
 
  struct node_gra **list;
  struct node_gra *root;
  long seed;

  list=(struct node_gra **)calloc(N_NOD+1,sizeof(struct node_gra *));
  COMS=(int *)calloc(N_NOD+1,sizeof(int));
  root=CreateHeaderGraph();
 
  for(i=0;i<N_NOD+1;i++){
    list[i]=NULL;
  }
  
  list[0]=CreateNodeGraph(root,0,0,0);
  for(i=1;i<N_NOD;i++){
    list[i]=CreateNodeGraph(list[i-1],i,i,i);
  }
  
  
  Z_in=AV_DEG-Z_out;
  
  P_in=(double)Z_in/(COM_SIZE-1);
  P_out=(double)Z_out/(N_NOD-COM_SIZE);
  
  
  for(i=0;i<N_NOD;i++){
    COMS[i]=i%N_COMS;
    //    printf("%d %d\n",i,i%N_COMS);
  }
  
  for(i=0;i<N_NOD;i++){
    for(j=i+1;j<N_NOD;j++){
      if(COMS[i]==COMS[j]){
	if(ran2(&seed)<P_in){
	  
	  AddAdjacency(list[i],j,0);
	  
	  AddAdjacency(list[j],i,0);
	}
      }
      else {
	if(ran2(&seed)<P_out){
	  
	  AddAdjacency(list[i],j,0);
	  
	  AddAdjacency(list[j],i,0);
	}
      }
    }
  }
  
  
  RewireAdjacency(root);
  
  return root;
}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateHierarchicalCommunityGraph(int N_NOD,int N_COMS1,int N_COMS2, int Z_in1, int Z_in2, int AV_DEG)
{
	int COMS[max_size];
	int i,j;
	double Z_out=AV_DEG-Z_in1-Z_in2;

	int COM_SIZE1 = N_NOD / N_COMS1;
	int COM_SIZE2 = COM_SIZE1 / N_COMS2;

	double P_in1,P_in2,P_out;

	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++) list[i]=NULL;

	list[0]=CreateNodeGraph(root,1,1,1);
	for(i=1;i<N_NOD;i++){
		list[i]=CreateNodeGraph(list[i-1],i+1,i+1,i+1);
	}


	P_in2=Z_in2/(double)(COM_SIZE2-1);
	P_in1=Z_in1/(double)(COM_SIZE1-COM_SIZE2);
	P_out=Z_out/(double)(N_NOD-COM_SIZE1);


	for(i=0;i<N_NOD;i++) {
		list[i]->state=i/(COM_SIZE2);
		COMS[i]=i/(COM_SIZE2);
	}


	for(i=0;i<N_NOD;i++){
		for(j=i+1;j<N_NOD;j++){
			if(COMS[i]==COMS[j]){
				if(ran2(&seed)<P_in2){
					AddAdjacency(list[i],j+1,0);
					AddAdjacency(list[j],i+1,0);
				}
			}else if(COMS[i]/N_COMS1 == COMS[j]/N_COMS1){
				if(ran2(&seed)<P_in1){
					AddAdjacency(list[i],j+1,0);
					AddAdjacency(list[j],i+1,0);
				}
			}else
			{
				if(ran2(&seed)<P_out){
					AddAdjacency(list[i],j+1,0);
					AddAdjacency(list[j],i+1,0);
				}
			}
		}
	}

	RewireAdjacency(root);

	return root;
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateSquareLattice(int N_FILS,int N_COLS)
{
	int i;
	int N_NOD = N_FILS*N_COLS;

	int x,y,xa,xp,ya,yp,node;

	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++) list[i]=NULL;

	list[0]=CreateNodeGraph(root,1,1,1);
	for(i=1;i<N_NOD;i++){
		list[i]=CreateNodeGraph(list[i-1],i+1,i+1,i+1);
	}

	for(i=0;i<N_NOD;i++) {
		x=i%N_COLS;
		y=i/N_COLS;
		list[i]->y= (float) (0.5f + y) / (float) (N_FILS+1);
		list[i]->x= (float) (0.5f + x) / (float) (N_COLS+1);
		list[i]->z= 0.5;

		ya=y-1;if(ya==-1) ya=N_FILS-1;
		xa=x-1;if(xa==-1) xa=N_COLS-1;
		yp=y+1;if(yp==N_FILS) yp=0;
		xp=x+1;if(xp==N_COLS) xp=0;

		sprintf(list[i]->nom,"%i",i+1);

		node=y*N_COLS + xa;
		node=node+1;
		AddAdjacency(list[i],node,0);

		node=y*N_COLS + xp;
		node=node+1;
		AddAdjacency(list[i],node,0);

		node=ya*N_COLS + x;
		node=node+1;
		AddAdjacency(list[i],node,0);

		node=yp*N_COLS + x;
		node=node+1;
		AddAdjacency(list[i],node,0);
	}

	RewireAdjacency(root);

	return root;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateHexagonalLattice(int N_FILS,int N_COLS)
{
	int i;
	int N_NOD = N_FILS*N_COLS;

	int x,y,xa,xp,ya,yp,node;

	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++) list[i]=NULL;

	list[0]=CreateNodeGraph(root,1,1,1);
	for(i=1;i<N_NOD;i++){
		list[i]=CreateNodeGraph(list[i-1],i+1,i+1,i+1);
	}

	for(i=0;i<N_NOD;i++) {
		x=i%N_COLS;
		y=i/N_COLS;
		list[i]->y= (float) ( (float) y/(N_FILS+1)) * 0.866f;
		if(y%2==0) list[i]->x= (float) (0.5f + x) / (float) (N_COLS+1);
		else list[i]->x= (float) x / (float) (N_COLS+1);

		list[i]->z=0.5;

		ya=y-1;if(ya==-1) ya=N_FILS-1;
		xa=x-1;if(xa==-1) xa=N_COLS-1;
		yp=y+1;if(yp==N_FILS) yp=0;
		xp=x+1;if(xp==N_COLS) xp=0;

		sprintf(list[i]->nom,"%i",i+1);

		if(y%2==0) {
			node=ya*N_COLS + x;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=ya*N_COLS + xp;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=y*N_COLS + xa;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=y*N_COLS + xp;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=yp*N_COLS + x;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=yp*N_COLS + xp;
			node=node+1;
			AddAdjacency(list[i],node,0);
		}else{
			node=ya*N_COLS + x;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=ya*N_COLS + xa;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=y*N_COLS + xa;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=y*N_COLS + xp;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=yp*N_COLS + x;
			node=node+1;
			AddAdjacency(list[i],node,0);

			node=yp*N_COLS + xa;
			node=node+1;
			AddAdjacency(list[i],node,0);
		}
	}

	RewireAdjacency(root);

	return root;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateScaleFree(int S,int m)		
	//S: Num nodes, m:min.links node segons article 
{
	int i,j;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;
	int des;
	double dau,cum;
	int control;
	int norm;
	int conn[max_size];

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
		conn[i]=0;
		list[i]=NULL;
	}

	//Generem n nodes
	list[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<S;i++){
		list[i]=CreateNodeGraph(list[i-1],i,i,0);
	}

	//Creem els m nodes inicials que estaran tots linkats entre ells. 
	for(i=0;i<m;i++){
		AddAdjacency(list[i],m,0);
		AddAdjacency(list[m],i,0);
		conn[i]++;
		conn[m]++;
	}
	norm=2*m;		//norm=num. total de links que tenim fins ara creats

	for(i=m+1;i<S;i++){		//La resta de nodes els afegim 1 a 1
		for(j=0;j<m;j++){		//Per cada nou node hi afegim el nou link.
			do{
				dau=(double)norm*ran2(&seed);
				cum=0.0;
				des=-1;
				do{
					des++;
					cum+=conn[des];
				}while(cum<dau);
				control=AddAdjacency(list[i],des,0);		//Posem control pq no esrepeteixin links
				control=AddAdjacency(list[des],i,0);
			}while(control==0);
			conn[i]++;
			conn[des]++;
			norm+=2;
		}
	}

	RewireAdjacency(root);

	return root;
}












///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateScaleFreeWithLocal(int S,int m,double p)
{
	int i,j;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;
	int ori,des;
	double dau,cum;
	int control;
	int norm;
	int conn[max_size];
	double daus[max_size];
	int whois[max_size];
	int posic[max_size];
	int itemp;
	double dtemp;
	int loclin;
	int loccon;


	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
		conn[i]=0;
	}

	for(i=0;i<S;i++){
		list[i]=CreateNodeGraph(root,i,i,0);
		daus[i]=(double)ran2(&seed);
		posic[i]=i;
		whois[i]=i;
	}

	//order the nodes according to their random value
	for(i=0;i<S-1;i++){
		for(j=i+1;j<S;j++){
			if(daus[i]>daus[j]){
				dtemp=daus[i];
				daus[i]=daus[j];
				daus[j]=dtemp;
				itemp=whois[posic[i]];
				whois[posic[i]]=whois[posic[j]];
				whois[posic[j]]=itemp;
				itemp=posic[i];
				posic[i]=posic[j];
				posic[j]=itemp;
			}
		}
	}    

	// build the initial core of the network
	for(i=0;i<m;i++){
		AddAdjacency(list[i],m,0);
		AddAdjacency(list[m],i,0);
		conn[i]+=1;
		conn[m]+=1;
	}
	norm=2*m;

	// build the preferential structure
	for(i=m+1;i<S;i++){
		for(j=0;j<m;j++){
			dau=(double)ran2(&seed);
			if(dau>p){
				do{
					dau=(double)norm*ran2(&seed);
					cum=0.0;
					des=-1;
					do{
						des++;
						cum+=conn[des];
					}while(cum<dau);
					control=AddAdjacency(list[i],des,0);
					control=AddAdjacency(list[des],i,0);
				}while(control==0);
				conn[i]+=1;
				conn[des]+=1;
				norm+=2;
			}
		}
	}

	// add the local structure
	loclin=(int) floor(p*(m*S));
	for(i=0;i<loclin;i++){
		ori=(int) floor(S*ran2(&seed));
		loccon=1;
		do{
			des=(ori+loccon)%S;
			control=AddAdjacency(list[whois[ori]],whois[des],1);
			control=AddAdjacency(list[whois[des]],whois[ori],1);
			loccon++;
		}while(control==0);
	}

	RewireAdjacency(root);

	//TANIS PrintPajekGraph1DList(root,posic);

	return root;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateSimonScaleFree(int S,double a)
{
	int i;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;
	double dau,cum;
	int control;
	int norm;
	int conn[max_size];
	int nnod;
	int nori,ndes;

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
		conn[i]=0;
	}

	// Create all the nodes
	list[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<S;i++){
		list[i]=CreateNodeGraph(list[i-1],i,i,0);
	}


	// Connect the 2 starting nodes
	AddAdjacency(list[0],1,1);
	AddAdjacency(list[1],0,0);
	conn[0]+=1;
	conn[1]+=1;
	norm=2;
	nnod=2;

	//create adjacencies
	while(nnod<S){
		dau=(double)ran2(&seed);
		if(dau<a){   // add new node & attach to randomly selected node
			ndes=(int) floor(nnod*ran2(&seed)); // select destination
			AddAdjacency(list[nnod],ndes,0);
			AddAdjacency(list[ndes],nnod,1);
			conn[nnod]+=1;
			conn[ndes]+=1;
			norm+=2;
			nnod++;
		}
		else{   // add a link preferentially
			nori=(int) floor(nnod*ran2(&seed)); //origin at random
			// destination preferentially
			do{
				dau=(double)ran2(&seed)*norm;
				cum=0.0;
				ndes=-1;
				do{
					ndes++;
					cum+=(double)conn[ndes];
				}while(cum<dau);
			}while(nori==ndes);
			control=AddAdjacency(list[nori],ndes,1);
			control=AddAdjacency(list[ndes],nori,0);

			conn[ndes]+=control;
			conn[nori]+=control;
			norm+=2*control;
		}
	}

	RewireAdjacency(root);

	return root;
}









///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreatePoissonRandomGraph(int S,int L)
{
	int i,des1,des2;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;
	int control;

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	list[0]=CreateNodeGraph(root,0,0,0);
	for(i=1;i<S;i++){
		list[i]=CreateNodeGraph(list[i-1],i,i,i);
	}

	printf("Nodes creats\n");

	for(i=0;i<L;i++){			//Situem els L links
		do{
			des1=(int) floor((S-1)*ran2(&seed));			//Tria dos nodes aleatoriament 1 i (S-1)
			des2=(int) floor((S-1)*ran2(&seed));
			control=AddAdjacency(list[des1],des2,0);		//I els enlla�a
			control=AddAdjacency(list[des2],des1,0);
		}while(control==0);
	}

	printf ("Links creats\n");
	RewireAdjacency(root);

	return root;
}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreatePoissonRandomGraphDirected(int S,int L)
{
	int i,des1,des2;
	struct node_gra *list[max_size];
	struct node_gra *root=NULL;
	int control;

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	for(i=0;i<S;i++){
		list[i]=CreateNodeGraph(root,i,i,0);
	}

	for(i=0;i<L;i++){
		do{
			des1=(int) floor((S-1)*ran2(&seed));			//Tria dos nodes aleatoriament
			des2=(int) floor((S-1)*ran2(&seed));
			control=AddAdjacency(list[des1],des2,0);
		}while(control==0);
	}

	RewireAdjacency(root);

	return root;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateHierarchicalNetwork(int m,int z)	
	// m: number of levels, z: branching factor
{
	int i,j,temp;
	struct node_gra *root;
	struct node_gra *vector[max_size];
	int kk;
	int nodes,start,end;

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		vector[i]=NULL;
	}

	vector[0]=CreateNodeGraph(root,0,0,0);

	nodes=1;
	end=0;
	for(i=1;i<m;i++){
		start=end+1;
		nodes=nodes*z;
		end=start+nodes-1;
		for (j=start;j<=end;j++){
			vector[j]=CreateNodeGraph(root,j,j,0);
			temp=(j-1)/z;
			kk=AddAdjacency(vector[j],temp,0);
			kk=AddAdjacency(vector[temp],j,0);
		}
	}

	RewireAdjacency(root);

	return root;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *Create1DSmallWorld(int n,int k,double p)
	//n: num. de la xarxa
	//k: num .de veins posteriors als que m'enganxare
	//p: probabilitat que enlloc del vei m'enllasi amb un qualsevol
{
	int i,j;
	int dest;
	struct node_gra *root=NULL;
	struct node_gra *list[max_size];


	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	for(i=0;i<n;i++){
		list[i]=CreateNodeGraph(root,i,i,0);
	}

	for(i=0;i<n;i++){    
		for(j=0;j<k;j++){
			if(ran2(&seed)>p){		
				AddAdjacency(list[i],(i+j+1)%n,0);
				AddAdjacency(list[(i+j+1)%n],i,0);
			}
			else{
				dest=(int) floor(n*ran2(&seed));
				AddAdjacency(list[i],dest,0);
				AddAdjacency(list[dest],i,0);	
			}
		}
	}

	RewireAdjacency(root);

	return root;
}




struct node_gra *CreateCommunityGraphWeighted(int N_NOD,int N_COMS, int MAX_ZIN, int MAX)
{
	int COMS[max_size];
	int i,j,temp;
	int MAX_ZOUT=MAX-MAX_ZIN;
	int COM_SIZE = N_NOD / N_COMS;
	int MAX_W_IN,MAX_W_OUT;


	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++) list[i]=NULL;

	list[0]=CreateNodeGraph(root,1,1,1);
	for(i=1;i<N_NOD;i++){
		list[i]=CreateNodeGraph(list[i-1],i+1,i+1,i+1);
	}


	MAX_W_IN=(int) (MAX_ZIN/(double)(COM_SIZE-1));
	MAX_W_OUT=(int) (MAX_ZOUT/(double)(N_NOD-COM_SIZE));


	for(i=0;i<N_NOD;i++) {
		list[i]->state=i%N_COMS;
		COMS[i]=i%N_COMS;
	}


	for(i=0;i<N_NOD;i++){
		for(j=i+1;j<N_NOD;j++){
			if(COMS[i]==COMS[j]){
				temp=(int) (ran2(&seed)*MAX_W_IN);
				if(temp>0){
					AddAdjacency(list[i],j+1,temp);
					AddAdjacency(list[j],i+1,temp);
				}
			}else {
				temp=(int) (ran2(&seed)*MAX_W_OUT);
				if(temp>0){
					AddAdjacency(list[i],j+1,temp);
					AddAdjacency(list[j],i+1,temp);
				}
			}
		}
	}

	RewireAdjacency(root);

	return root;
}





struct node_gra *CreateCommunityGraphWeightedNewman(int N_NOD,int N_COMS, double Z_out, double AV_DEG, float W_IN)
{
	int COMS[max_size];
	int i,j;
	double Z_in=AV_DEG-Z_out;
	int COM_SIZE = N_NOD / N_COMS;
	double P_in,P_out;

	struct node_gra *list[max_size];
	struct node_gra *root=NULL;

	seed+=-1*time(NULL);

	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++) list[i]=NULL;

	list[0]=CreateNodeGraph(root,1,1,1);
	for(i=1;i<N_NOD;i++){
		list[i]=CreateNodeGraph(list[i-1],i+1,i+1,i+1);
	}


	P_in=Z_in/(double)(COM_SIZE-1);
	P_out=Z_out/(double)(N_NOD-COM_SIZE);


	for(i=0;i<N_NOD;i++) {
		list[i]->state=i%N_COMS;
		COMS[i]=i%N_COMS;
	}


	for(i=0;i<N_NOD;i++){
		for(j=i+1;j<N_NOD;j++){
			if(COMS[i]==COMS[j]){
				if(ran2(&seed)<P_in){
					AddAdjacency(list[i],j+1,(int)(W_IN * 10));		//Creem el graf igual pero ara amb pesos
					AddAdjacency(list[j],i+1,(int)(W_IN * 10));
				}
			}else {
				if(ran2(&seed)<P_out){
					AddAdjacency(list[i],j+1,10);
					AddAdjacency(list[j],i+1,10);
				}
			}
		}
	}

	RewireAdjacency(root);

	return root;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
struct node_gra *CreateInternetModel(int n,double r)
	// Atencio... parametres del model en:
	//  p->state ...... # usuaris del node (tamany node)
	//  p->state_out .. # links del nodes
	//  p->state_in ... # links q falten per enlla�ar
{
	int i,j,k,l;
	int dest,orig,ant;
	struct node_gra *root=NULL;
	struct node_gra *list[max_size];

	double a_b= 1 / 0.8;
	double d_b= 1 / 0.8; 
	/*    double alfa=0.035;
	double beta=0.03;
	double delta=0.04;
	double a_b= alfa / beta;
	double d_b= delta / beta;*/
	double a_t;
	double num_users,nous_users,min_users=5000.0f;
	double num_links,num_links_esperats,num_links_pendents;
	double num_nodes,num_nodes_links_pendents;
	double angle,dist,distx,disty,k_dist;
	double maxx,maxy,minx,miny;
	double prob,prob_rand;
	double max_dist=1;
	int final,control;

	FILE *fppf;

	seed+=-1*time(NULL);

	fppf=fopen("info_inet.txt","w");

	k_dist= 1 - pow((double) max_dist,(double) -0.5f);

	//Creem la cap�alera de la nostra xarxa
	root=CreateHeaderGraph();

	for(i=0;i<max_size;i++){
		list[i]=NULL;
	}

	//Creem els 2 nodes inicials amb 5000 users i els linkem
	list[0]=CreateNodeGraph(root,0,0,0);        //Node 0
	list[0]->state= (int) min_users;
	list[0]->state_out=1;
	list[0]->x=max_dist; maxx=max_dist; minx=max_dist;
	list[0]->y=max_dist; maxy=max_dist; miny=max_dist;

	list[1]=CreateNodeGraph(root,1,1,0);        //Node 1
	list[1]->state= (int) min_users;
	list[1]->state_out=1;
	angle=ran2(&seed)*6.2831853;
	dist=pow( 1 - ran2(&seed)*k_dist, -2.0); 
	list[1]->x = list[0]->x + dist * sin(angle);
	list[1]->y = list[0]->y + dist * cos(angle);
	ant=1;
	if(list[1]->x<minx) minx=list[1]->x;
	if(list[1]->x>maxx) maxx=list[1]->x;
	if(list[1]->y<miny) miny=list[1]->y;
	if(list[1]->y>maxy) maxy=list[1]->y;

	AddAdjacencyWeight(list[0],1);                  //Link
	AddAdjacencyWeight(list[1],0);

	num_users=2*min_users;
	num_nodes=2;
	num_links=1;

	fprintf(fppf,"NODE\t#users\tNous_us\t#Links\tNous_li\n");

	//Creem la resta de nodes de la xarxa seguint el model del marian
	for(i=2;i<n;i++){
		if(i%1000==999) printf("*"); else printf(".");
		fprintf(fppf,"%i\t",i);
		list[i]=CreateNodeGraph(root,i,i,0);
		list[i]->state_out=0;
		list[i]->state=(int) min_users;

		//Creem els nous usuaris i els assignem als nodes JA EXISTENTS
		nous_users= min_users * pow ((double) (i+1) / 2.0f, a_b) - num_users;
		if( nous_users <0 ) nous_users=0;

		fprintf(fppf,"%i\t%i\t",(int) num_users,(int) nous_users);
		for(j=0;j<nous_users;j++)
		{
			//Triem a quin node s'enganxara el nou usuari;
			dest=(int) (ran2(&seed)*num_users); l=0; k=0;
			while(l<dest)
			{
				l += list[k]->state;
				if(l<dest) k++;
			}

			list[k]->state++;
			num_users++;
		}


		//Robem min_users usuaris als nodes existents per a posar-los en el nou usuari
		for(j=0;j<min_users;j++)
		{
			dest=(int) (ran2(&seed)*i);
			if( list[dest]->state > 1) list[dest]->state--;
			else j--;
		}

		//Calculem la posici� x,y del nou node. Mirem si cal fer un salt o si continuem el mateix cami
		angle=ran2(&seed)*6.2831853;//2pi
		dist=pow( 1 - ran2(&seed)*k_dist, -2.0); 
		list[i]->x = list[ant]->x + dist * sin(angle);
		list[i]->y = list[ant]->y + dist * cos(angle);
		if(list[i]->x<minx) minx=list[i]->x;
		if(list[i]->x>maxx) maxx=list[i]->x;
		if(list[i]->y<miny) miny=list[i]->y;
		if(list[i]->y>maxy) maxy=list[i]->y;
		ant=i;
		if(ran2(&seed) > 0.99 ) ant=(int) (ran2(&seed)*i);


		num_links_esperats = pow ((double) (i+1) / 2.0f, d_b );

		//a_t = (2.0 * num_links_esperats) / (num_users);
		a_t = ((2.0 * num_links_esperats) - i) / ((num_users) - (min_users * i)); if(a_t<0) a_t=0;
		num_links_pendents=0;
		num_links=0;
		//Recalculem els valors dels bandwiths de tots els nodes
		for(j=0;j<=i;j++)
		{
			list[j]->state_in = (int) (1 + a_t * (list[j]->state - min_users)); //OJO AMB LA PART ENTERA I REAL!!
			list[j]->state_in -= list[j]->state_out; 
			if(list[j]->state_in<0) list[j]->state_in=0;
			num_links_pendents += list[j]->state_in;
			num_links += list[j]->state_out;
		}

		//Generem tots els nous links que facin falta. Aturarem quan no es pugui formar un link degut a:
		//O be que queden 0 o 1 links per enlla�ar, o be en queden mes pero tots pertanyen al mateix node
		final=0; control=0;
		fprintf(fppf,"%i\t%i\n",(int) num_links,(int) num_links_pendents);
		while(num_links_pendents>=2 && !final)
		{
			//Comprovem q hi hagi mes d'un node amb links per enlla�ar
			num_nodes_links_pendents=0; 
			for(j=0;j<=i;j++) { if(list[j]->state_in>0) num_nodes_links_pendents++; }
			if( num_nodes_links_pendents<2 ) final=1;
			control++; if( control>100 ) final=1;      //Aixo indica q durant 100 intents no hem pogut crear cap link


			if(!final) //Es pot fer l'enlla� pq hi ha 2 o mes nodes per enlla�ar.
			{
				//Seleccionem els dos nodes a canviar
				orig=(int) (ran2(&seed)*num_links_pendents)+1; l=0; k=0;
				while(l<orig)
				{
					l += list[k]->state_in;
					if(l<orig) k++;
				}
				dest=orig=k;

				while(dest==orig)
				{
					dest=(int) (ran2(&seed)*num_links_pendents)+1; l=0; k=0;
					while(l<dest)
					{
						l += list[k]->state_in;
						if(l<dest) k++;
					}
					dest=k;
				}

				//Calculem la probabilitat q es formi un link entre els dos nodes triats
				distx= (list[orig]->x - list[dest]->x) / (maxx-minx);
				disty= (list[orig]->y - list[dest]->y) / (maxy-miny);
				dist= sqrt( pow(distx,2) + pow(disty ,2));

				prob = (list[orig]->state * list[dest]->state) / ( 100.0f * num_users);  //OJO AMB EL VALOR DEL COST
				prob = dist / prob;
				prob = exp(-prob);

				prob_rand=ran2(&seed);
				if(prob_rand<prob)
				{
					control=0;
					//Fem el primer link
					AddAdjacencyWeight(list[orig],dest);                  //Link
					list[orig]->state_out++;
					list[orig]->state_in--;
					AddAdjacencyWeight(list[dest],orig);
					list[dest]->state_out++;
					list[dest]->state_in--;
					num_links_pendents-=2;

					//I mirem si es poden fer mes enlla�os amb probabilitat r
					prob_rand=ran2(&seed);
					while(list[orig]->state_in>0 && list[dest]->state_in>0 && prob_rand<r)
					{
						AddAdjacencyWeight(list[orig],dest);                  //Link
						list[orig]->state_out++;
						list[orig]->state_in--;
						AddAdjacencyWeight(list[dest],orig);
						list[dest]->state_out++;
						list[dest]->state_in--;
						prob_rand=ran2(&seed);
						num_links_pendents-=2;
					}
				}
			}
		}  
	}
	//Reescalem les distancies entre 0 i 1
	for(i=0;i<n;i++)
	{
		sprintf(list[i]->nom,"%i",list[i]->state);
		list[i]->x = (list[i]->x - minx) / (maxx-minx);
		list[i]->y = (list[i]->y - miny) / (maxy-miny);
	}

	RewireAdjacency(root);

	fclose(fppf);
	return root;
}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void PositionGraph(struct node_gra *net,int nnodes)
{
	//Serveix per a assignar valors de x i y als nodes del graf
	//Util per si es representa graficament
	struct node_gra *p=net;
	float angle;

	while(p->next!=NULL)
	{
		p=p->next;
		angle=(float)2.0*(float)3.14159*(float)p->num/(float)nnodes;
		p->x=(cos(angle)*0.49)+(0.5);
		p->y=(sin(angle)*0.49)+(0.5);
		p->z=0.5;
	}
}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void IniNodeNames(struct node_gra *net)
{
	//Serveix per a assignar valors de x i y als nodes del graf
	//Util per si es representa graficament
	struct node_gra *p=net;

	while(p->next!=NULL)
	{
		p=p->next;
		sprintf(p->nom,"*%i*",p->num);
	}
}


















