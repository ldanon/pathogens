#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "network_lib.h"
#include "bitops.cpp"
#include "rand_gen.h"
#include "globals.h"
#include "pathogens.h"




double CalculateAverage(char infname[50],int cutoff){
  
  int itmp;
  double dtmp,sum=0;
  FILE *inFile=fopen(infname,"r");
  
  while(!feof(inFile)){
    fscanf(inFile,"%d %lf\n",&itmp,&dtmp);
    if(itmp>cutoff)sum+=dtmp;
  }
  fclose(inFile);

  return (double)sum/(ITERATIONS-cutoff);

}

parameters *InitialiseParameters(double p_init,double p_infect, 
				 double p_recover,double p_mutate,
				 double p_loss_imm,double p_recom,
				 double gamma,
				 double k, double theta,double k2){
  
  parameters *par;
  par=(parameters *)malloc(sizeof(parameters));
  par->p_init=p_init;
  par->p_infect=p_infect;
  par->p_recover=p_recover;
  par->p_loss_imm=p_loss_imm;

  par->p_mutate=p_mutate;
  par->p_recom=p_recom;
  par->gamma=gamma;
  
  par->k=k;
  par->theta=theta;
  par->k2=k2;  	  
  
//  SetSusceptibilityTimesMatrix(par); //probably won't need.
  UpdateVulnerabilityVector(par);
  
  return par;
  
}

void SetSusceptibilityTimesMatrix(parameters *par){
	int i,j;
	int dist;
	int num_strains=NUM_STRAINS;
	for(i=0;i<num_strains;i++){
		for(j=0;j<num_strains;j++){
			dist=HammingDist(i,j);
			switch (dist){
			case 0:
				par->susc_times[i][j]=0;
				break;

			case 1:
				par->susc_times[i][j]=1;
				break;

			case 2:
				par->susc_times[i][j]=2;
				break;
	
			case 3:
				printf("Hamming Distance more than 2 problem!!\n");
				break;
			case 4:
				printf("Hamming Distance more than 3 problem!!\n");
				break;

			}
		}
	}

}


void UpdateVulnerabilityVector(parameters *par){
  int i;
  double gamma=par->gamma;
  double fraction;
  for(i=0;i<NUM_LOCI;i++){
    fraction=(double)i/NUM_LOCI;
    par->vulnerability[i]=1.0-pow(1-pow(fraction,(double)1/gamma),gamma);
  }
  par->vulnerability[NUM_LOCI]=1;
     
 
       for(i=0;i<=NUM_LOCI;i++){
       printf("gamma %lf vulnvector %d %lf\n",gamma,i,par->vulnerability[i]);
       } 

  
}

strain *CreateStrain(int num){

  strain *str;
  str=(strain *)malloc(sizeof(strain));
  str->timer=0;
  str->length=0;
  str->string=num;
  str->next=NULL;
  return str;
  
}

strain *CreateStrainTimestamp(int num, int T){
	strain *str;
	str=(strain *)malloc(sizeof(strain));
	str->timer=T;
	str->string=num;
	str->next=NULL;
	return str;
}

/*+++++++++++++++++ STRAINS ++++++++++++++++++++++++++ */

void StartInfectionStrains(struct node_gra *net,
			   double p_init,int start_strain){
  struct node_gra *n;
  
  n=net;
  while(n->next!=NULL){
    n=n->next;
    n->infected=CreateStrain(0);
    n->immune=CreateStrain(0);
    n->queue=CreateStrain(0);

    if(p_init>ran2(&seed)){
      n->queue->string=1;
      n->queue->next=CreateStrain(start_strain);
    }
  }
}

void StartInfection2Strains(struct node_gra *net,
			   double p_init,int start_strain1,int start_strain2){
  struct node_gra *n;
  n=net;
  while(n->next!=NULL){
    n=n->next;
    n->infected=CreateStrain(0);
    n->immune=CreateStrain(0);
    n->queue=CreateStrain(0);

    if(p_init>ran2(&seed)){
      n->queue->string=1;
      if(ran2(&seed)>0.5)
	 n->queue->next=CreateStrain(start_strain1);
      else n->queue->next=CreateStrain(start_strain2);
    }
  }
}

void StartInfection2COMS(struct node_gra *net,
			   double p_init){
  struct node_gra *n;
  
  //printf("starting infection in two coms\n");
  n=net;
  while(n->next!=NULL){
    n=n->next;
    n->infected=CreateStrain(0);
    n->immune=CreateStrain(0);
    n->queue=CreateStrain(0);
    
    if(n->num%NUM_COM==0){
      if(p_init>ran2(&seed)){
	n->queue->string=1;
	if(ran2(&seed)>0.5)
	  n->queue->next=CreateStrain(0);
	else n->queue->next=CreateStrain(3);
      }
    }
    else{
      if(p_init>ran2(&seed)){
		n->queue->string=1;
		if(ran2(&seed)>0.5)
	  		n->queue->next=CreateStrain(1);
	  else n->queue->next=CreateStrain(2);
      }

    }
  }
}


void StartInfectionDengue(struct node_gra *net,
							parameters *params,gsl_rng *rngen)
{
	struct node_gra *n;
	double rndm;
	//printf("starting infection in two coms\n");
	n=net;
	while(n->next!=NULL){
		n=n->next;
		n->infected=CreateStrain(0); // set up headers
		n->immune=CreateStrain(0);
		n->queue=CreateStrain(0);

//		if(n->num==1){
		if(params->p_init>ran2(&seed)){
			n->queue->string=1;
			rndm=ran2(&seed);
			
			if(rndm<0.25){
				n->queue->next=CreateStrain(0);
				n->queue->next->length=gsl_ran_gamma(rngen, params->k,params->theta);
			}
			else if(rndm<0.5){
				n->queue->next=CreateStrain(1);
				n->queue->next->length=gsl_ran_gamma(rngen, params->k,params->theta);
			}
			else if(rndm<0.75){
				n->queue->next=CreateStrain(2);
				n->queue->next->length=gsl_ran_gamma(rngen, params->k,params->theta);
			}
			else {
				n->queue->next=CreateStrain(3);
				n->queue->next->length=gsl_ran_gamma(rngen, params->k,params->theta);		}
		}
		
	}
}


void MutateStrain(strain *s,parameters *params){
  int i;
  for(i=0;i<NUM_LOCI;i++){
    if(params->p_mutate>ran2(&seed)){
      s->string=BitFlp(s->string,i);
    }
  }
}



int RecombineStrains(struct node_gra *n,parameters *params){
  
  int string1,string2;
  
  
  int i=0;
  
  string1=n->infected->next->string;
  string2=n->infected->next->next->string;
  
  //  printf("1 %d 2 %d\n",string1,string2);
  //  printf("before %d ",string1);
  
  for(i=0;i<NUM_LOCI;i++){
    if(BitTst(string1,i)!=BitTst(string2,i)){
      if(params->p_recom>ran2(&seed)){
	//printf("recombining! ");
	string1=BitFlp(string1,i);
      }
    }
  }
  return string1;
}
 

int InfectedWithStrain(struct node_gra *n,int num){
		/* returns 1 if agent is infected with strain  */
		/* and 0 if not NEW: also includes the queue*/
		/* if there is something in the queue, return 1 */

  strain *s;

  s=n->infected;
  while(s->next!=NULL){
    if(s->next->string==num)return 1;
    s=s->next;
  }
  return 0;
}

int InfectedWithAnything(struct node_gra *n){
		/* returns 1 if agent is infected with strain  */
		/* and 0 if not NEW: also includes the queue*/
		/* if there is something in the queue, return 1 */

	if(n->infected->next!=NULL) return 1;
	else  return 0;
}


int IsSomethingInQueue(struct node_gra *n){
  if(n->queue->next==NULL)return 0;
  else return 1;
}
 


int ImmuneToStrain(struct node_gra *n,int num){
		/* returns 1 if agent is immune to strain */
 				/* returns 0 if it is not */
  strain *s;
  
  s=n->immune;
  while(s->next!=NULL){
    if(s->string==num)return 1;
    s=s->next;
  }
  return 0;
}



void UpdateInfectionStrains(struct node_gra *net,parameters *params){
	struct node_gra *n;
	strain *s_q,*s_inf;

	n=net;
	while(n->next!=NULL){
		n=n->next;
		s_q=n->queue;
		if(s_q->next!=NULL){
			//    printf("n %d ",n->num);

			s_inf=n->infected;
			while(s_inf->next!=NULL){
				printf("--> %d ",s_inf->next->string);
				s_inf=s_inf->next;
			}
			//printf("\n");
			s_inf->next=s_q->next;
			s_inf->next->next=NULL;
		}
		n->queue->next=NULL;
	}
}


void UpdateInfectionStrainsDengue(struct node_gra *net,parameters *params){
	struct node_gra *n;
	strain *s_q,*s_inf;
	int i=0,j=0;

	n=net;
	while(n->next!=NULL){
		n=n->next;
		s_q=n->queue;
		i=0;

		while(s_q->next!=NULL){
			i++;		
			s_q=s_q->next;
		}
 		s_q=n->queue;
		
		
		
		if(i==1){
			s_inf=n->infected;
			while(s_inf->next!=NULL){
				s_inf=s_inf->next;
			}

			s_inf->next=s_q->next;
			s_inf->next->next=NULL;
			s_q->next=NULL;
		}
		
		if(i>1){


			j=(int)floor(ran2(&seed)*i);
			s_inf=n->infected;
			while(j>0){
				j--;
				s_q=s_q->next;		
			}
			s_inf->next=s_q->next;
			s_q->next=s_q->next->next;

			s_inf->next->next=NULL;
			
		}
		
		if(i>0)n->queue->next=ClearStrains(n->queue);
		
	}
}



void UpdateInfectionStrainsInstantImmunity(struct node_gra *net,
					   parameters *params){
  struct node_gra *n;
  strain *s_q,*s_inf;
  strain *s_imm;

  n=net;
  while(n->next!=NULL){
	  n=n->next;
	  s_q=n->queue;
	  s_imm=n->immune;
	  if(s_q->next!=NULL){
		  s_inf=n->infected;
		  while(s_inf->next!=NULL){
			  s_inf=s_inf->next;
		  }
		  s_inf->next=s_q->next;
		  s_inf->next->next=NULL;

		  while(s_imm->next!=NULL){	/* instant immunity */
			  s_imm=s_imm->next;
		  }
		  s_imm->next=CreateStrain(s_inf->next->string);

		  n->queue->next=NULL;
	  }
  }
}
 
 
void InfectNeighboursStrains(struct node_gra *n, parameters *params){
   /* picks a strain at random from infectious node, */
   /* checks whether neightbours are infected/immune */  
   /* if not, infect neighbours with probability  */
   
//  struct node_gra *n1;
  struct node_lis *l;
  strain *s,*s_ra[10],*s_ri[10];
  strain *s_imm,*s_q;
  int string;
  
  int i,j;
//  int found,found2;
  double vuln;
  
  l=n->neig;
  
  s=n->infected;
  i=0;
  while(s->next!=NULL){
    s=s->next;
    s_ra[i]=s;
    i++;
  }

  j=0;
  s_imm=n->immune;
  while(s_imm->next!=NULL){
    s_imm=s_imm->next;
    s_ri[j]=s_imm;
    j++;
  }
  
  
  if(i==0)printf("something is wrong %d\n",n->infected->string);

                   /* pick strain at random*/
  
  if(i==1)string=n->infected->next->string;
  

  if(i>1){
/*     printf(" %d more than one strain found, %d %d ", */
/* 	      i,       s_ra[0]->string,    s_ra[1]->string); */
    string=RecombineStrains(n,params);
/*     printf("recombinant %d\n",string); */
  }
  
/*   if(i>2){ */
/*     printf("marker i %d %d %d %d %d \n",i,s_ra[0]->string, */
/* 		s_ra[1]->string,s_ra[2]->string/\* ,s_ra[3]->string *\/); */
/*   } */


  while(l->next!=NULL){		/* go through neighbours */

    l=l->next;
    if(IsSomethingInQueue(l->ref)==0){
               /* if there is nothing in the queue */
      if(InfectedWithStrain(l->ref,string)==0){ 
	       /* if it's not already infected with strain*/
	vuln=CalculateVulnerabilityOfNode(l->ref,string,params); 
	                  //if strains share alele,0, else 1
	
	if((params->p_infect*vuln)>ran2(&seed)){
	
	  s_q=l->ref->queue;		     /* go to queue */
	  			
	  s_q->next=CreateStrain(string); /* copy strain to end of queue */
	  s_q->next->next=NULL;
	}
      
      }
    
    }
    
  }

}


void InfectStrainsMetapop(struct node_gra *n, parameters *params, 
struct node_gra *nodes_array[N_NODE]){
	/* picks a strain at random from infectious node, */
	/* checks whether random node is infected/immune */  
	/* if not, infect random node with probability  */

//	struct node_gra *n1;
	strain *s,*s_ra[10],*s_ri[10];
	strain *s_imm,*s_q;
	int string;

	int i,j;
//	int found2;
	double vuln;
	double prob_temp;

	s=n->infected;
	i=0;
	while(s->next!=NULL){
		s=s->next;
		s_ra[i]=s;
		i++;
	}

	j=0;
	s_imm=n->immune;
	while(s_imm->next!=NULL){
		s_imm=s_imm->next;
		s_ri[j]=s_imm;
		j++;
	}


	if(i==0)printf("something is wrong %d\n",n->infected->string);

	/* pick strain at random*/

	if(i==1)string=n->infected->next->string;


	if(i>1){
		/*     printf(" %d more than one strain found, %d %d ", */
		/* 	      i,       s_ra[0]->string,    s_ra[1]->string); */
		string=RecombineStrains(n,params);
		/*     printf("recombinant %d\n",string); */
	}

	/*   if(i>2){ */
	/*     printf("marker i %d %d %d %d %d \n",i,s_ra[0]->string, */
	/* 		s_ra[1]->string,s_ra[2]->string/\* ,s_ra[3]->string *\/); */
	/*   } */


	//while(l->next!=NULL){		/* go through neighbours */
	for(i=0;i<N_NODE;i++){	/* go through all nodes*/

		if((IsSomethingInQueue(nodes_array[i])==0)&&(i!=n->num)){
			/* if there is nothing in the queue */
			if(InfectedWithStrain(nodes_array[i],string)==0){ 
				/* if it's not already infected with strain*/
				vuln=CalculateVulnerabilityOfNode(nodes_array[i],string,params); 
				//if strains share alele,0, else 1
				if(n->state==nodes_array[i]->state) prob_temp=params->p_in;
				else prob_temp=params->p_out;
				if((prob_temp*params->p_infect*vuln)>ran2(&seed)){

					s_q=nodes_array[i]->queue;		     /* go to queue */

					s_q->next=CreateStrain(string); /* copy strain to end of queue */
					s_q->next->next=NULL;
				}

			}

		}

	}

}


void InfectStrainsNoStructureDengueStrainSpecificCrossImmunity(struct node_gra *n, parameters *params, 
									struct node_gra *nodes_array[N_NODE],int T,gsl_rng *rngen){
	/* picks a strain at random from infectious node, */
	/* checks whether random node is infected/immune */  
	/* if not, infect random node with probability  */
	/* specifically deals with four strains two of which have a longer */
	/* cross immune period than the other two */

	struct node_gra *n1;
	strain *s,*s_ra[10],*s_ri[10];
	strain *s_inf,*s_imm,*s_q;
	int string;

	int i,j;
	double vuln;
	double randnum;
	int ni=0;

	s=n->infected;
	i=0;
	while(s->next!=NULL){
		s=s->next;
		s_ra[i]=s;
		i++;
	}

	j=0;
	s_imm=n->immune;
	while(s_imm->next!=NULL){
		s_imm=s_imm->next;
		s_ri[j]=s_imm;
		j++;
	}


	if(i==0){
		printf("something is wrong %d\n",n->infected->string);
	}
	/* pick strain at random*/

	if(i==1)string=n->infected->next->string;


	if(i>1){
		printf("%d WRONG!. more than one strain found,T %d %d %d %d", n->num,T, 
				i, s_ra[0]->string,    s_ra[1]->string); 
//		string=RecombineStrains(n,params);
	}

	//if(i>2){ 
	//	printf("PROBLEM!!!!! marker i %d %d %d %d %d \n",i,s_ra[0]->string, 
	//		s_ra[1]->string,s_ra[2]->string ,s_ra[3]->string ); 
	//} 


	//while(l->next!=NULL){		/* go through neighbours */


	ni=gsl_ran_binomial(rngen,params->p_infect,N_NODE);
	while(ni>0){
		i=(int)floor(ran2(&seed)*N_NODE);

		if(i!=n->num){
			ni--;
			if(InfectedWithAnything(nodes_array[i])==0){ 
				/* if it's not already infected*/
								
				
				
//				if(n->state==nodes_array[i]->state) prob_temp=params->p_in; //THIS IS FOR INTRODUCING COMMUNITY STRUCTURE

				
				if(CalculateVulnerabilityOfNodeDengue(nodes_array[i],string,params,T)==1){
//					printf("string %d %d %lf %lf \n",i,string,vuln,randnum);
					s_q=nodes_array[i]->queue;		     /* go to queue */
					while(s_q->next!=NULL)s_q=s_q->next;	
					s_q->next=CreateStrain(string); /* copy strain to end of queue */
					s_q->next->length=gsl_ran_gamma(rngen,params->k,params->theta);
					s_q->next->timer=T;
					//s_q->next->timer=gsl_ran_gamma(rngen,params->k,params->theta);
					s_q->next->next=NULL;
//					printf("NODE %d \n created new strain %d at time T= %d, immunity wanes at T=%d\n",
//						          i,                    string,          T,          s_q->next->timer);
//					printf("blah\n");
				}

			}

		}

	}

}

void InfectStrainsNoStructureDengueStrainSpecificCrossImmunityAndTransmissionEnhancement(struct node_gra *n, parameters *params, 
									struct node_gra *nodes_array[N_NODE],int T,gsl_rng *rngen){
	/* picks a strain at random from infectious node, */
	/* checks whether random node is infected/immune */  
	/* if not, infect random node with probability  */
	/* specifically deals with four strains two of which have a longer */
	/* cross immune period than the other two */
	strain *s,*s_ra[10],*s_ri[10];
	strain *s_imm,*s_q;
	int string;

	int i,j;
	double vuln;

	double enhanced_p=params->p_infect;
	int ni=0;

	s=n->infected;
	i=0;
	while(s->next!=NULL){
		s=s->next;
		s_ra[i]=s;
		i++;
	}

	j=0;
	s_imm=n->immune;
	while(s_imm->next!=NULL){
		s_imm=s_imm->next;
		s_ri[j]=s_imm;
		j++;
	}


	if(i==0){
		printf("something is wrong %d\n",n->infected->string);
	}
	/* pick strain at random*/

	if(i==1)string=n->infected->next->string;


	if(i>1){
		printf("%d WRONG!. more than one strain found,T %d %d %d %d", n->num,T, 
				i, s_ra[0]->string,    s_ra[1]->string); 
	}

	//if(i>2){ 
	//	printf("PROBLEM!!!!! marker i %d %d %d %d %d \n",i,s_ra[0]->string, 
	//		s_ra[1]->string,s_ra[2]->string ,s_ra[3]->string ); 
	//} 

	if(n->immune->next!=NULL){
		if(string==1||string==0)
		enhanced_p=params->p_infect*params->enhancement;
	}
	

	if(params->p_mutate>ran2(&seed)){
		string=(int) floor(ran2(&seed)*NUM_STRAINS);
	}

		
	ni=gsl_ran_binomial(rngen,enhanced_p,N_NODE);
	
	while(ni>0){
		i=(int)floor(ran2(&seed)*N_NODE);

		if(i!=n->num){
			ni--;
			if(InfectedWithAnything(nodes_array[i])==0){ 
				/* if it's not already infected*/
				
				vuln=CalculateVulnerabilityOfNodeDengue(nodes_array[i],string,params,T);
				

//				randnum=ran2(&seed);
				if((vuln)>0){

//					printf("string %d %d %lf %lf \n",i,string,vuln,randnum);
					s_q=nodes_array[i]->queue;		     /* go to queue */
					while(s_q->next!=NULL)s_q=s_q->next;	
					s_q->next=CreateStrain(string); /* copy strain to end of queue */
					s_q->next->length=gsl_ran_gamma(rngen,params->k,params->theta);
					s_q->next->timer=T;

					//s_q->next->timer=gsl_ran_gamma(rngen,params->k,params->theta);
					s_q->next->next=NULL;

//					printf("NODE %d \n created new strain %d at time T= %d, immunity wanes at T=%d\n",
//						          i,                    string,          T,          s_q->next->timer);
//					printf("blah\n");
				}

			}

		}

	}

}

void InfectStrainsNoStructureDengueStrainSpecificCrossImmunityAndSusceptibilityEnhancement(struct node_gra *n, parameters *params, 
									struct node_gra *nodes_array[N_NODE],int T,gsl_rng *rngen){
	/* picks a strain at random from infectious node, */
	/* checks whether random node is infected/immune */  
	/* if not, infect random node with probability  */
	/* specifically deals with four strains two of which have a longer */
	/* cross immune period than the other two */


	strain *s,*s_ra[10],*s_ri[10];
	strain *s_imm,*s_q;
	int string;

	int i,j;
	double vuln;
	
	double enhanced_p=params->p_infect;
	int ni=0;

	s=n->infected;
	i=0;
	while(s->next!=NULL){
		s=s->next;
		s_ra[i]=s;
		i++;
	}

	j=0;
	s_imm=n->immune;
	while(s_imm->next!=NULL){
		s_imm=s_imm->next;
		s_ri[j]=s_imm;
		j++;
	}


	if(i==0){
		printf("something is wrong %d\n",n->infected->string);
	}
	/* pick strain at random*/

	if(i==1)string=n->infected->next->string;


	if(i>1){
		printf("%d WRONG!. more than one strain found,T %d %d %d %d", n->num,T, 
				i, s_ra[0]->string,    s_ra[1]->string); 
	}

	//if(i>2){ 
	//	printf("PROBLEM!!!!! marker i %d %d %d %d %d \n",i,s_ra[0]->string, 
	//		s_ra[1]->string,s_ra[2]->string ,s_ra[3]->string ); 
	//} 


	//while(l->next!=NULL){		/* go through neighbours */

	ni=gsl_ran_binomial(rngen,params->p_infect,N_NODE);

	while(ni>0){
		i=(int)floor(ran2(&seed)*N_NODE);

		if(i!=n->num){
			ni--;
			if(InfectedWithAnything(nodes_array[i])==0){ 
				/* if it's not already infected*/
				
				if(CalculateVulnerabilityOfNodeDengue(nodes_array[i],string,params,T)==1){

//					printf("string %d %d %lf %lf \n",i,string,vuln,randnum);
					s_q=nodes_array[i]->queue;		     /* go to queue */
					while(s_q->next!=NULL)s_q=s_q->next;	
					s_q->next=CreateStrain(string); /* copy strain to end of queue */
					s_q->next->length=gsl_ran_gamma(rngen,params->k,params->theta);
					s_q->next->timer=T;

					//s_q->next->timer=gsl_ran_gamma(rngen,params->k,params->theta);
					s_q->next->next=NULL;

//					printf("NODE %d \n created new strain %d at time T= %d, immunity wanes at T=%d\n",
//						          i,                    string,          T,          s_q->next->timer);
//					printf("blah\n");
				}

			}

		}

	}

}

void InfectStrainsNoStructureDengue(struct node_gra *n, parameters *params, 
									struct node_gra *nodes_array[N_NODE],int T,gsl_rng *rngen){
	/* picks a strain at random from infectious node, */
	/* checks whether random node is infected/immune */  
	/* if not, infect random node with probability  */
	/* specifically deals with four strains two of which have a longer */
	/* cross immune period than the other two */

	struct node_gra *n1;
	strain *s,*s_ra[10],*s_ri[10];
	strain *s_inf,*s_imm,*s_q;
	int string;

	int i,j;
	int found,found2;
	double vuln;
	double prob_temp;
	double randnum;


	s=n->infected;
	i=0;
	while(s->next!=NULL){
		s=s->next;
		s_ra[i]=s;
		i++;
	}

	j=0;
	s_imm=n->immune;
	while(s_imm->next!=NULL){
		s_imm=s_imm->next;
		s_ri[j]=s_imm;
		j++;
	}


	if(i==0){
		printf("something is wrong %d\n",n->infected->string);
	}
	/* pick strain at random*/

	if(i==1)string=n->infected->next->string;


	if(i>1){
		printf("%d WRONG!. more than one strain found,T %d %d %d %d", n->num,T, 
				i, s_ra[0]->string,    s_ra[1]->string); 
//		string=RecombineStrains(n,params);
	}

	//if(i>2){ 
	//	printf("PROBLEM!!!!! marker i %d %d %d %d %d \n",i,s_ra[0]->string, 
	//		s_ra[1]->string,s_ra[2]->string ,s_ra[3]->string ); 
	//} 


	//while(l->next!=NULL){		/* go through neighbours */


	for(i=0;i<N_NODE;i++){
		/* go through all nodes*/

// more efficient to pick a number from a binomial distribution with mean p_infect*number of nodes
// then pick that number of nodes randomly from nodes_array[i]. (int)(floor(ran2(&seed))*N_NODE)).

		if(i!=n->num){

			if(InfectedWithAnything(nodes_array[i])==0){ 
				/* if it's not already infected*/
				vuln=CalculateVulnerabilityOfNodeDengue(nodes_array[i],string,params,T); 

//				if(n->state==nodes_array[i]->state) prob_temp=params->p_in; //THIS IS FOR INTRODUCING COMMUNITY STRUCTURE

//				else prob_temp=params->p_out;
				//randnum=ran2(&seed);
				if((vuln)==1){
//					printf("string %d %d %lf %lf \n",i,string,vuln,randnum);
					s_q=nodes_array[i]->queue;		     /* go to queue */
					while(s_q->next!=NULL)s_q=s_q->next;	
					s_q->next=CreateStrain(string); /* copy strain to end of queue */
					if(string==0||string==2)s_q->next->timer=T+gsl_ran_gamma(rngen,params->k,params->theta);
					else s_q->next->timer=T+gsl_ran_gamma(rngen,params->k2,params->theta);
					s_q->next->next=NULL;
//					printf("NODE %d \n created new strain %d at time T= %d, immunity wanes at T=%d\n",
//						          i,                    string,          T,          s_q->next->timer);
//					printf("blah\n");
				}

			}

		}

	}

}

double AsymptoticHeaviside(int k, int x, int sharp){
	
	//asymptotic function approaching the Heaviside step function for large values of slope
	
	double ans;
	ans=1.0/double(1+exp((double)-sharp*(k-x)));
	return ans;
	
}

int InfectNeighboursStrainsCountCrossing(struct node_gra *n,
										 parameters *params){
											 /* picks a strain at random from infectious node, */
											 /* checks whether neightbours are infected/immune */  
											 /* if not, infect neighbours with probability  */

  struct node_gra *n1;
  struct node_lis *l;
  strain *s,*s_ra[10],*s_ri[10];
  strain *s_inf,*s_imm,*s_q;
  int string;
  
  int i,j;
  int found,found2;
  double vuln;
  int cross12,cross21;
  
  l=n->neig;
  
  s=n->infected;
  i=0;
  while(s->next!=NULL){
    s=s->next;
    s_ra[i]=s;
    i++;
  }

  j=0;
  s_imm=n->immune;
  while(s_imm->next!=NULL){
    s_imm=s_imm->next;
    s_ri[j]=s_imm;
    j++;
  }
  
  
  if(i==0)printf("something is wrong %d\n",n->infected->string);

  
  if(i==1)string=n->infected->next->string;
  

  if(i>1){
/*     printf(" %d more than one strain found, %d %d ", */
/* 	      i,       s_ra[0]->string,    s_ra[1]->string); */
    string=RecombineStrains(n,params);
/*     printf("recombinant %d\n",string); */
  }
  
/*   if(i>2){ */
/*     printf("marker i %d %d %d %d %d \n",i,s_ra[0]->string, */
/* 		s_ra[1]->string,s_ra[2]->string/\* ,s_ra[3]->string *\/); */
/*   } */


  while(l->next!=NULL){		/* go through neighbours */

    l=l->next;
    if(IsSomethingInQueue(l->ref)==0){
               /* if there is nothing in the queue */
      if(InfectedWithStrain(l->ref,string)==0){ 
	       /* if it's not already infected with strain*/
	vuln=CalculateVulnerabilityOfNode(l->ref,string,params); 
	                  //if strains share alele,0, else 1
	
	if((params->p_infect*vuln)>ran2(&seed)){
	
	  s_q=l->ref->queue;		     /* go to queue */
	  			
	  s_q->next=CreateStrain(string); /* copy strain to end of queue */
	  s_q->next->next=NULL;
	}
      
      }
    
    }
    
  }
return 0;
}




int CalculateVulnerabilityOfNodeDengue(struct node_gra *n, int s1,parameters *params, int T){

	strain *s_imm;
	double vuln;

	


	s_imm=n->immune;


	if(s_imm->next==NULL){
		//printf("fully susceptible, virgin immune system\n");
		return 1;
	}

	else if(s_imm->next->next==NULL){ // if there is only one in the immune pile

//		printf("T= %d     has been exposed to %d \n",T,s_imm->next->string);

		if(s1==s_imm->next->string){
			//printf("trying to infect with strain that the host is already fully immune to\n");
			return 0;
		}

		else {				
				vuln=1.0/params->vulnerability[HammingDist(s1,s_imm->next->string)];
				//			vu
				if((s_imm->next->timer + (int)vuln*(s_imm->next->length)) > T){
				//printf("trying to infect before short term cross immunity has waned\n");
					return 0;
				}
				else {
				//printf("cross immunity has waned, secondary infection possible\n");
				return 1;	
				}
			 }
			
		}
	
	
	else if(s_imm->next->next->next!=NULL){
		printf("has been exposed to 3 things and is VERY VERY WRONG\n");
		return -1;
	}
	else {
		//printf("has been exposed to %d and %d so is fully immune\n",s_imm->next->string,s_imm->next->next->string);
		return 0;
	}
}

double CalculateVulnerabilityOfNode(struct node_gra *n,
			            int s1,parameters *params){
			/* calculate vulnerability of node n to strain s */
  int i;
  int same;
  double fraction;
  double vuln,vuln_min;
  double gamma=params->gamma;
  strain *s_imm;
  int s2;
  
  s_imm=n->immune;
  
  vuln_min=1;
/*   printf("vuln %d %d ",n->num,s1); */
  
  while(s_imm->next!=NULL){
    s_imm=s_imm->next;
    same=0;
    s2=s_imm->string;
/*     printf("s2 %d ",s2); */
    for(i=0;i<NUM_LOCI;i++){
      if(BitTst(s1,i)==BitTst(s2,i)){
	same++;
/* 	printf("found same alele "); */
      }
    }
/*     if(same!=0){ */
/* /\*       printf("same %d ",same); *\/ */
/*       fraction=(double)same/NUM_LOCI; */
/*       vuln=pow(1-pow(fraction,(double)1/gamma),gamma); */
/* /\*       printf("vuln %lf, gamma %lf ",vuln,gamma); *\/ */
/*     } */
/*     else vuln=1; */

    vuln=params->vulnerability[same];
    
    if(vuln<vuln_min){
      vuln_min=vuln;
    }
  }
/*   printf("returning %lf\n",vuln_min); */
  return vuln_min;
}

double CalculateVulnerabilityOfNodeTest(struct node_gra *n,
			            int s1,parameters *params){
			/* calculate vulnerability of node n to strain s */
  int i;
  int same;
  double fraction;
  double vuln,vuln_min;
  double gamma=params->gamma;
  strain *s_imm;
  int s2;
  
  //printf("infecting %d ",s1);
  s_imm=n->immune;
  while(s_imm->next!=NULL){
    s_imm=s_imm->next;
    same=0;
    s2=s_imm->string;
    //printf("string in node %d ",s2);
    for(i=0;i<NUM_LOCI;i++){
      if(BitTst(s1,i)==BitTst(s2,i)){
	//printf("returning 0\n");
	return 0.0;
      }
    }
  }
  //  printf("returning 1\n");

  return 1.0;
  
}
void RunStrainsOnNetwork(struct node_gra *net,parameters *params,
			 int iterations){
  int i;
  int count;
  struct node_gra *n,*n1;
  double prob;
  
  strain *s;
  
  int com;
  int start_strain=2;
  int ctr;

  int inf_1,inf_2,inf_3,inf_4,not_inf;
  int imm_1,imm_2,imm_3,imm_4,not_imm;
 
  FILE *inf1,*inf2,*inf3,*inf4,*notinf;
  FILE *imm1,*imm2,*imm3,*imm4,*notimm;

  int check,check2;
  
  inf1=fopen("inf1.dat","w");
  inf2=fopen("inf2.dat","w");
  inf3=fopen("inf3.dat","w");
  inf4=fopen("inf4.dat","w");
  notinf=fopen("not_inf.dat","w");

  imm1=fopen("imm1.dat","w");
  imm2=fopen("imm2.dat","w");
  imm3=fopen("imm3.dat","w");
  imm4=fopen("imm4.dat","w");
  notimm=fopen("not_imm.dat","w");

  n=net;
  
  StartInfection2Strains(n,params->p_init,1,2); /* start with strain1*/
                                                  /* and strain2 */
  
  UpdateInfectionStrains(n,params);
/*   while(n->next!=NULL){ */
/*     n=n->next; */
/*     if(n->infected->next!=NULL)printf("initial infection there %d %d\n", */
/* 				      n->num,n->infected->next->string); */
/*   } */
  
  for(i=0;i<iterations;i++){
    
    
    //printf("timestep %d %d\n",i,i%1000);
    
    inf_1=inf_2=inf_3=inf_4=not_inf=0;
    imm_1=imm_2=imm_3=imm_4=not_imm=0;
    check=0;

    n=net;
    while(n->next!=NULL){	/* go throught network */
      
      n=n->next;
      s=n->immune;
      if(s->next==NULL)not_imm++;
      check2=0;
      while(s->next!=NULL){
	check++;
	check2++;
	switch(s->next->string){
	case 0:
	  imm_1++;
	  break;
	case 1:
	  imm_2++;
	  break;
	case 2:
	  imm_3++;
	  break;
	case 3:
	  imm_4++;
	  break;
	}
	s=s->next;
      }

      ctr=0;
      s=n->infected;
      if(s->next==NULL)not_inf++;
      
      while(s->next!=NULL){
	switch(s->next->string){
	  
	case 0:
	  inf_1++;
	  break;
	case 1:
	  inf_2++;
	  break;
	case 2:
	  inf_3++;
	  break;
	case 3:
	  inf_4++;
	  break;
	}

	MutateStrain(s->next,params);
	s=s->next;
	ctr++;
      }
      
      if(ctr>0){
	
	InfectNeighboursStrains(n,params);
      }
      //      if(ctr>2)printf("timestep %d check %d check2 %d\n",i,ctr,check2);
  
    
    }
    
    fprintf(inf1,"%d %d\n",i,inf_1);
    fprintf(inf2,"%d %d\n",i,inf_2);
    fprintf(inf3,"%d %d\n",i,inf_3);
    fprintf(inf4,"%d %d\n",i,inf_4);
    fprintf(notinf,"%d %d\n",i,not_inf);
    
    fprintf(imm1,"%d %d\n",i,imm_1);
    fprintf(imm2,"%d %d\n",i,imm_2);
    fprintf(imm3,"%d %d\n",i,imm_3);
    fprintf(imm4,"%d %d\n",i,imm_4);
    fprintf(notimm,"%d %d\n",i,not_imm);

    if(check==(imm_1+imm_2+imm_3+imm_4))
      printf("checkbig  %d %d\n",
	     check,imm_1+imm_2+imm_3+imm_4);

    RecoverInfectedLoseImmunityStrains(net,params);
    UpdateInfectionStrains(net,params); 

  }

  fclose(inf1);
  fclose(inf2);
  fclose(inf3);
  fclose(inf4);

  fclose(imm1);
  fclose(imm2);
  fclose(imm3);
  fclose(imm4);
  
}

void RunStrainsOnComsNetwork(struct node_gra *net,parameters *params,
			     int iterations){
  int i;
  int count;
  struct node_gra *n,*n1;
  double prob;
  
  strain *s;
  
  int com;
  int start_strain=2;
  int ctr;

  int inf1_1,inf1_2,inf1_3,inf1_4;
  int inf2_1,inf2_2,inf2_3,inf2_4;

  int inf1,inf2,inf3,inf4;

  int imm1_1,imm1_2,imm1_3,imm1_4;
  int imm2_1,imm2_2,imm2_3,imm2_4;

  int coinf0,coinf1,coinf2,coinf3;

  int not_inf,not_imm;

  double discordance1, diversity1, discordance2, diversity2;
  double discordance_tot, diversity_tot;
  double disc_sum,div_sum;
 
  FILE *div_totf=fopen("diversity_tot.dat","w");
  FILE *disc_totf=fopen("discordance_tot.dat","w");
  
  FILE *div1f=fopen("diversity1.dat","w");
  FILE *disc1f=fopen("discordance1.dat","w");
  
  FILE *div2f=fopen("diversity2.dat","w");
  FILE *disc2f=fopen("discordance2.dat","w");

  FILE *inf11=fopen("inf11.dat","w");
  FILE *inf12=fopen("inf12.dat","w");
  FILE *inf13=fopen("inf13.dat","w");
  FILE *inf14=fopen("inf14.dat","w");

  FILE *inf21=fopen("inf21.dat","w");
  FILE *inf22=fopen("inf22.dat","w");
  FILE *inf23=fopen("inf23.dat","w");
  FILE *inf24=fopen("inf24.dat","w");

  FILE *imm11=fopen("imm11.dat","w");
  FILE *imm12=fopen("imm12.dat","w");
  FILE *imm13=fopen("imm13.dat","w");
  FILE *imm14=fopen("imm14.dat","w");

  FILE *imm21=fopen("imm21.dat","w");
  FILE *imm22=fopen("imm22.dat","w");
  FILE *imm23=fopen("imm23.dat","w");
  FILE *imm24=fopen("imm24.dat","w");

  FILE *notimm=fopen("not_imm.dat","w");
  FILE *notinf=fopen("not_inf.dat","w");

  FILE *coinf0f=fopen("co_inf_0.dat","w");
  FILE *coinf1f=fopen("co_inf_1.dat","w");
  FILE *coinf2f=fopen("co_inf_2.dat","w");
  FILE *coinf3f=fopen("co_inf_3.dat","w");
 
  FILE *trackinf1f=fopen("node_infect1.dat","w");
  FILE *trackinf2f=fopen("node_infect2.dat","w");
  FILE *trackinf3f=fopen("node_infect3.dat","w");
  FILE *trackinf4f=fopen("node_infect4.dat","w");

  FILE *trackimm1f=fopen("node_immune1.dat","w");
  FILE *trackimm2f=fopen("node_immune2.dat","w");
  FILE *trackimm3f=fopen("node_immune3.dat","w");
  FILE *trackimm4f=fopen("node_immune4.dat","w");

  n=net;
  
/*   StartInfection2Strains(n,params->p_init,1,2); */ /* start with strain1*/
                                                  /* and strain2 */
  
   StartInfection2COMS(n,params->p_init); 
/*   UpdateInfectionStrainsInstantImmunity(n,params); */
    UpdateInfectionStrains(n,params);
/*   while(n->next!=NULL){ */
/*     n=n->next; */
/*     if(n->infected->next!=NULL)printf("initial infection there %d %d\n", */
/* 				      n->num,n->infected->next->string); */
/*   } */

  for(i=0;i<iterations;i++){
    
    //if(i%2==0)printf("timestep %i\n",i);

    inf1_1=inf1_2=inf1_3=inf1_4=inf2_1=inf2_2=inf2_3=inf2_4=0;

    imm1_1=imm1_2=imm1_3=imm1_4=imm2_1=imm2_2=imm2_3=imm2_4=0;
    
    coinf0=coinf1=coinf2=coinf3=0;
    
    not_inf=not_imm=0;

    n=net;
    //PrintStrainClusterFile(net);
    while(n->next!=NULL){	/* go throught network */
      n=n->next;
      ctr=0;
      if(n->num==TRACK)PrintNodeStatus2(trackinf1f,trackinf2f,
					trackinf3f,trackinf4f,
					trackimm1f,trackimm2f,
					trackimm3f,trackimm4f,
					i,n);


      if((n->num)%NUM_COM==0){ /* COMMUNITY 1 */
	s=n->immune;
	if(s->next==NULL)not_imm++;
	while(s->next!=NULL){
	  switch(s->next->string){
	  case 0:
	    imm1_1++;
	    break;
	  case 1:
	    imm1_2++;
	    break;
	  case 2:
	    imm1_3++;
	    break;
	  case 3:
	    imm1_4++;
	    break;
	  }
	  s=s->next;
	}
	
	s=n->infected;
	if(s->next==NULL)not_inf++;
	
	while(s->next!=NULL){	
	  switch(s->next->string){
	    
	  case 0:
	    inf1_1++;
	    break;
	  case 1:
	    inf1_2++;
	    break;
	  case 2:
	    inf1_3++;
	    break;
	  case 3:
	    inf1_4++;
	    break;
	  }
	  
	  MutateStrain(s->next,params);
	  s=s->next;
	  ctr++;
	}
      }

      else{                       /* COMMUNITY 2 */
	s=n->immune;
	if(s->next==NULL)not_imm++;
	
	while(s->next!=NULL){
	  switch(s->next->string){
	  case 0:
	    imm2_1++;
	    break;
	  case 1:
	    imm2_2++;
	    break;
	  case 2:
	    imm2_3++;
	    break;
	  case 3:
	    imm2_4++;
	    break;
	  }
	  s=s->next;
	}

	s=n->infected;
	if(s->next==NULL)not_inf++;
	
	while(s->next!=NULL){
	  switch(s->next->string){
	    
	  case 0:
	    inf2_1++;
	    break;
	  case 1:
	    inf2_2++;
	    break;
	  case 2:
	    inf2_3++;
	    break;
	  case 3:
	    inf2_4++;
	    break;
	  }
	  
	  MutateStrain(s->next,params);
	  s=s->next;
	  ctr++;
	}
      }
      if(ctr>0){
	InfectNeighboursStrains(n,params);
	coinf0++;
      }
      
      if(ctr>1)coinf1++;
      if(ctr>2)coinf2++;
      if(ctr>3)coinf3++;
    
    }
    fprintf(coinf0f,"%d %d\n",i,coinf0);
    fprintf(coinf1f,"%d %d\n",i,coinf1);
    fprintf(coinf2f,"%d %d\n",i,coinf2);
    fprintf(coinf3f,"%d %d\n",i,coinf3);
    
    fprintf(inf11,"%d %d\n",i,inf1_1);
    fprintf(inf12,"%d %d\n",i,inf1_2);
    fprintf(inf13,"%d %d\n",i,inf1_3);
    fprintf(inf14,"%d %d\n",i,inf1_4);

    fprintf(inf21,"%d %d\n",i,inf2_1);
    fprintf(inf22,"%d %d\n",i,inf2_2);
    fprintf(inf23,"%d %d\n",i,inf2_3);
    fprintf(inf24,"%d %d\n",i,inf2_4);

    
    fprintf(imm11,"%d %d\n",i,imm1_1);
    fprintf(imm12,"%d %d\n",i,imm1_2);
    fprintf(imm13,"%d %d\n",i,imm1_3);
    fprintf(imm14,"%d %d\n",i,imm1_4);

    fprintf(imm21,"%d %d\n",i,imm2_1);
    fprintf(imm22,"%d %d\n",i,imm2_2);
    fprintf(imm23,"%d %d\n",i,imm2_3);
    fprintf(imm24,"%d %d\n",i,imm2_4);

    fprintf(notimm,"%d %d\n",i,not_imm);
    fprintf(notinf,"%d %d\n",i,not_inf);
    
    diversity1=MeasureDiversity(inf1_1,inf1_2,inf1_3,inf1_4);
    discordance1=MeasureDiscordance(inf1_1,inf1_2,inf1_3,inf1_4);
    
    diversity2=MeasureDiversity(inf2_1,inf2_2,inf2_3,inf2_4);
    discordance2=MeasureDiscordance(inf2_1,inf2_2,inf2_3,inf2_4);
    
    fprintf(div1f,"%d %lf\n",i,diversity1);
    fprintf(disc1f,"%d %lf\n",i,discordance1);
    
    fprintf(div2f,"%d %lf\n",i,diversity2);
    fprintf(disc2f,"%d %lf\n",i,discordance2);


    inf1=inf1_1+inf2_1;
    inf2=inf1_2+inf2_2;
    inf3=inf1_3+inf2_3;
    inf4=inf1_4+inf2_4;
    
    if((inf1+inf2+inf3+inf4)==0){
      printf("infection died, returning\n");
      return ;
    }

    diversity_tot=MeasureDiversity(inf1,inf2,inf3,inf4);
    discordance_tot=MeasureDiscordance(inf1,inf2,inf3,inf4);
    
    fprintf(div_totf,"%d %lf\n",i,diversity_tot);
    fprintf(disc_totf,"%d %lf\n",i,discordance_tot);
    

    RecoverInfectedLoseImmunityStrains(net,params);
/*     UpdateInfectionStrainsInstantImmunity(net,params);  */
    UpdateInfectionStrains(net,params); 
#if defined( VTKNET )
	printf("REAL: %d\n",i);
	vtknet_draw_strains(net);
	getchar();
	//	Sleep(35);
#endif
  }

  fclose(inf11);
  fclose(inf12);
  fclose(inf13);
  fclose(inf14);

  fclose(inf21);
  fclose(inf22);
  fclose(inf23);
  fclose(inf24);

  fclose(imm11);
  fclose(imm12);
  fclose(imm13);
  fclose(imm14);
  
  fclose(imm21);
  fclose(imm22);
  fclose(imm23);
  fclose(imm24);

  fclose(notimm);
  fclose(notinf);
  
  fclose(div1f);
  fclose(disc1f);

  fclose(div2f);
  fclose(disc2f);

  fclose(div_totf);
  fclose(disc_totf);

  fclose(coinf1f);
  fclose(coinf2f);
  fclose(coinf3f);
  
  fclose(trackinf1f);
  fclose(trackinf2f);
  fclose(trackinf3f);
  fclose(trackinf4f);

  fclose(trackimm1f);
  fclose(trackimm2f);
  fclose(trackimm3f);
  fclose(trackimm4f);

}

void RunStrainsOnMetapop(struct node_gra *net,parameters *params,
						  int iterations, struct node_gra *nodes_array[N_NODE]){
  int i;
  int count;
  struct node_gra *n,*n1;
  double prob;
  
  strain *s;
  
  int com;
  int start_strain=2;
  int ctr;

  int inf1_1,inf1_2,inf1_3,inf1_4;
  int inf2_1,inf2_2,inf2_3,inf2_4;

  int inf1,inf2,inf3,inf4;

  int imm1_1,imm1_2,imm1_3,imm1_4;
  int imm2_1,imm2_2,imm2_3,imm2_4;

  int coinf0,coinf1,coinf2,coinf3;

  int not_inf,not_imm;

  double discordance1, diversity1, discordance2, diversity2;
  double discordance_tot, diversity_tot;
  double disc_sum,div_sum;
 
  FILE *div_totf=fopen("diversity_tot.dat","w");
  FILE *disc_totf=fopen("discordance_tot.dat","w");
  
  FILE *div1f=fopen("diversity1.dat","w");
  FILE *disc1f=fopen("discordance1.dat","w");
  
  FILE *div2f=fopen("diversity2.dat","w");
  FILE *disc2f=fopen("discordance2.dat","w");

  FILE *inf11=fopen("inf11.dat","w");
  FILE *inf12=fopen("inf12.dat","w");
  FILE *inf13=fopen("inf13.dat","w");
  FILE *inf14=fopen("inf14.dat","w");

  FILE *inf21=fopen("inf21.dat","w");
  FILE *inf22=fopen("inf22.dat","w");
  FILE *inf23=fopen("inf23.dat","w");
  FILE *inf24=fopen("inf24.dat","w");

  FILE *imm11=fopen("imm11.dat","w");
  FILE *imm12=fopen("imm12.dat","w");
  FILE *imm13=fopen("imm13.dat","w");
  FILE *imm14=fopen("imm14.dat","w");

  FILE *imm21=fopen("imm21.dat","w");
  FILE *imm22=fopen("imm22.dat","w");
  FILE *imm23=fopen("imm23.dat","w");
  FILE *imm24=fopen("imm24.dat","w");

  FILE *notimm=fopen("not_imm.dat","w");
  FILE *notinf=fopen("not_inf.dat","w");

  FILE *coinf0f=fopen("co_inf_0.dat","w");
  FILE *coinf1f=fopen("co_inf_1.dat","w");
  FILE *coinf2f=fopen("co_inf_2.dat","w");
  FILE *coinf3f=fopen("co_inf_3.dat","w");
 
  FILE *trackinf1f=fopen("node_infect1.dat","w");
  FILE *trackinf2f=fopen("node_infect2.dat","w");
  FILE *trackinf3f=fopen("node_infect3.dat","w");
  FILE *trackinf4f=fopen("node_infect4.dat","w");

  FILE *trackimm1f=fopen("node_immune1.dat","w");
  FILE *trackimm2f=fopen("node_immune2.dat","w");
  FILE *trackimm3f=fopen("node_immune3.dat","w");
  FILE *trackimm4f=fopen("node_immune4.dat","w");

  n=net;
  
/*   StartInfection2Strains(n,params->p_init,1,2); */ /* start with strain1*/
                                                  /* and strain2 */
  
   StartInfection2COMS(n,params->p_init); 
/*   UpdateInfectionStrainsInstantImmunity(n,params); */
   UpdateInfectionStrains(n,params);
/*   while(n->next!=NULL){ */
/*     n=n->next; */
/*     if(n->infected->next!=NULL)printf("initial infection there %d %d\n", */
/* 				      n->num,n->infected->next->string); */
/*   } */

  for(i=0;i<iterations;i++){
    
    if(i%200==0)printf("timestep %i\n",i);

    inf1_1=inf1_2=inf1_3=inf1_4=inf2_1=inf2_2=inf2_3=inf2_4=0;

    imm1_1=imm1_2=imm1_3=imm1_4=imm2_1=imm2_2=imm2_3=imm2_4=0;
    
    coinf0=coinf1=coinf2=coinf3=0;
    
    not_inf=not_imm=0;

    n=net;
    //PrintStrainClusterFile(net);
    while(n->next!=NULL){	/* go throught network */
      n=n->next;
      ctr=0;
      if(n->num==TRACK)PrintNodeStatus2(trackinf1f,trackinf2f,
					trackinf3f,trackinf4f,
					trackimm1f,trackimm2f,
					trackimm3f,trackimm4f,
					i,n);


      if((n->num)%NUM_COM==0){ /* COMMUNITY 1 */
	s=n->immune;
	if(s->next==NULL)not_imm++;
	while(s->next!=NULL){
	  switch(s->next->string){
	  case 0:
	    imm1_1++;
	    break;
	  case 1:
	    imm1_2++;
	    break;
	  case 2:
	    imm1_3++;
	    break;
	  case 3:
	    imm1_4++;
	    break;
	  }
	  s=s->next;
	}
	
	s=n->infected;
	if(s->next==NULL)not_inf++;
	
	while(s->next!=NULL){	
	  switch(s->next->string){
	    
	  case 0:
	    inf1_1++;
	    break;
	  case 1:
	    inf1_2++;
	    break;
	  case 2:
	    inf1_3++;
	    break;
	  case 3:
	    inf1_4++;
	    break;
	  }
	  
	  MutateStrain(s->next,params);
	  s=s->next;
	  ctr++;
	}
      }

      else{                       /* COMMUNITY 2 */
	s=n->immune;
	if(s->next==NULL)not_imm++;
	
	while(s->next!=NULL){
	  switch(s->next->string){
	  case 0:
	    imm2_1++;
	    break;
	  case 1:
	    imm2_2++;
	    break;
	  case 2:
	    imm2_3++;
	    break;
	  case 3:
	    imm2_4++;
	    break;
	  }
	  s=s->next;
	}

	s=n->infected;
	if(s->next==NULL)not_inf++;
	
	while(s->next!=NULL){
	  switch(s->next->string){
	    
	  case 0:
	    inf2_1++;
	    break;
	  case 1:
	    inf2_2++;
	    break;
	  case 2:
	    inf2_3++;
	    break;
	  case 3:
	    inf2_4++;
	    break;
	  }
	  
	  MutateStrain(s->next,params);
	  s=s->next;
	  ctr++;
	}
      }
      if(ctr>0){
		  InfectStrainsMetapop(n,params,nodes_array);
		  coinf0++;
      }
      
      if(ctr>1)coinf1++;
      if(ctr>2)coinf2++;
      if(ctr>3)coinf3++;
    
    }
    fprintf(coinf0f,"%d %d\n",i,coinf0);
    fprintf(coinf1f,"%d %d\n",i,coinf1);
    fprintf(coinf2f,"%d %d\n",i,coinf2);
    fprintf(coinf3f,"%d %d\n",i,coinf3);
    
    fprintf(inf11,"%d %d\n",i,inf1_1);
    fprintf(inf12,"%d %d\n",i,inf1_2);
    fprintf(inf13,"%d %d\n",i,inf1_3);
    fprintf(inf14,"%d %d\n",i,inf1_4);

    fprintf(inf21,"%d %d\n",i,inf2_1);
    fprintf(inf22,"%d %d\n",i,inf2_2);
    fprintf(inf23,"%d %d\n",i,inf2_3);
    fprintf(inf24,"%d %d\n",i,inf2_4);

    
    fprintf(imm11,"%d %d\n",i,imm1_1);
    fprintf(imm12,"%d %d\n",i,imm1_2);
    fprintf(imm13,"%d %d\n",i,imm1_3);
    fprintf(imm14,"%d %d\n",i,imm1_4);

    fprintf(imm21,"%d %d\n",i,imm2_1);
    fprintf(imm22,"%d %d\n",i,imm2_2);
    fprintf(imm23,"%d %d\n",i,imm2_3);
    fprintf(imm24,"%d %d\n",i,imm2_4);

    fprintf(notimm,"%d %d\n",i,not_imm);
    fprintf(notinf,"%d %d\n",i,not_inf);
    
    diversity1=MeasureDiversity(inf1_1,inf1_2,inf1_3,inf1_4);
    discordance1=MeasureDiscordance(inf1_1,inf1_2,inf1_3,inf1_4);
    
    diversity2=MeasureDiversity(inf2_1,inf2_2,inf2_3,inf2_4);
    discordance2=MeasureDiscordance(inf2_1,inf2_2,inf2_3,inf2_4);
    
    fprintf(div1f,"%d %lf\n",i,diversity1);
    fprintf(disc1f,"%d %lf\n",i,discordance1);
    
    fprintf(div2f,"%d %lf\n",i,diversity2);
    fprintf(disc2f,"%d %lf\n",i,discordance2);


    inf1=inf1_1+inf2_1;
    inf2=inf1_2+inf2_2;
    inf3=inf1_3+inf2_3;
    inf4=inf1_4+inf2_4;
    
    if((inf1+inf2+inf3+inf4)==0){
      printf("infection died, returning\n");
      return ;
    }

    diversity_tot=MeasureDiversity(inf1,inf2,inf3,inf4);
    discordance_tot=MeasureDiscordance(inf1,inf2,inf3,inf4);
    
    fprintf(div_totf,"%d %lf\n",i,diversity_tot);
    fprintf(disc_totf,"%d %lf\n",i,discordance_tot);
    

    RecoverInfectedLoseImmunityStrains(net,params);
/*     UpdateInfectionStrainsInstantImmunity(net,params);  */
    UpdateInfectionStrains(net,params); 
#if defined( VTKNET )
	printf("REAL: %d\n",i);
	vtknet_draw_strains(net);
	getchar();
	//	Sleep(35);
#endif
  }

  fclose(inf11);
  fclose(inf12);
  fclose(inf13);
  fclose(inf14);

  fclose(inf21);
  fclose(inf22);
  fclose(inf23);
  fclose(inf24);

  fclose(imm11);
  fclose(imm12);
  fclose(imm13);
  fclose(imm14);
  
  fclose(imm21);
  fclose(imm22);
  fclose(imm23);
  fclose(imm24);

  fclose(notimm);
  fclose(notinf);
  
  fclose(div1f);
  fclose(disc1f);

  fclose(div2f);
  fclose(disc2f);

  fclose(div_totf);
  fclose(disc_totf);

  fclose(coinf0f);
  fclose(coinf1f);
  fclose(coinf2f);
  fclose(coinf3f);
  
  fclose(trackinf1f);
  fclose(trackinf2f);
  fclose(trackinf3f);
  fclose(trackinf4f);

  fclose(trackimm1f);
  fclose(trackimm2f);
  fclose(trackimm3f);
  fclose(trackimm4f);

}


void RunStrainsOnMetapopDengue(struct node_gra *net,parameters *params,
						  int iterations, struct node_gra *nodes_array[N_NODE],
						  int inf_data[NUM_STRAINS+1][ITERATIONS],int imm_data[NUM_STRAINS+1][ITERATIONS]){
  int i;
  int count;
  struct node_gra *n,*n1;
  double prob;
  
  strain *s;

  strain *q; //for debugging
  
  int com;
  int start_strain=2;
  int ctr;

  int inf1_1,inf1_2,inf1_3,inf1_4;
  int inf2_1,inf2_2,inf2_3,inf2_4;

  int inf1,inf2,inf3,inf4;

  int imm1_1,imm1_2,imm1_3,imm1_4;
  int imm2_1,imm2_2,imm2_3,imm2_4;

  int coinf0,coinf1,coinf2,coinf3;

  int not_inf,not_imm;

  //double discordance1, diversity1, discordance2, diversity2;
  //double discordance_tot, diversity_tot;
  //double disc_sum,div_sum;
 
  //FILE *div_totf=fopen("diversity_tot.dat","w");
  //FILE *disc_totf=fopen("discordance_tot.dat","w");
  //
  //FILE *div1f=fopen("diversity1.dat","w");
  //FILE *disc1f=fopen("discordance1.dat","w");
  //
  //FILE *div2f=fopen("diversity2.dat","w");
  //FILE *disc2f=fopen("discordance2.dat","w");

  //FILE *inf11=fopen("inf11.dat","w");
  //FILE *inf12=fopen("inf12.dat","w");
  //FILE *inf13=fopen("inf13.dat","w");
  //FILE *inf14=fopen("inf14.dat","w");

  //FILE *inf21=fopen("inf21.dat","w");
  //FILE *inf22=fopen("inf22.dat","w");
  //FILE *inf23=fopen("inf23.dat","w");
  //FILE *inf24=fopen("inf24.dat","w");

  //FILE *imm11=fopen("imm11.dat","w");
  //FILE *imm12=fopen("imm12.dat","w");
  //FILE *imm13=fopen("imm13.dat","w");
  //FILE *imm14=fopen("imm14.dat","w");

  //FILE *imm21=fopen("imm21.dat","w");
  //FILE *imm22=fopen("imm22.dat","w");
  //FILE *imm23=fopen("imm23.dat","w");
  //FILE *imm24=fopen("imm24.dat","w");

	FILE *notimm=fopen("not_imm.dat","w");
	FILE *notinf=fopen("not_inf.dat","w");
 /*   
  FILE *coinf0f=fopen("co_inf_0.dat","w");
  FILE *coinf1f=fopen("co_inf_1.dat","w");
  FILE *coinf2f=fopen("co_inf_2.dat","w");
  FILE *coinf3f=fopen("co_inf_3.dat","w");
 
  FILE *trackinf1f=fopen("node_infect1.dat","w");
  FILE *trackinf2f=fopen("node_infect2.dat","w");
  FILE *trackinf3f=fopen("node_infect3.dat","w");
  FILE *trackinf4f=fopen("node_infect4.dat","w");

  FILE *trackimm1f=fopen("node_immune1.dat","w");
  FILE *trackimm2f=fopen("node_immune2.dat","w");
  FILE *trackimm3f=fopen("node_immune3.dat","w");
  FILE *trackimm4f=fopen("node_immune4.dat","w");*/

	gsl_rng * rngen = gsl_rng_alloc (gsl_rng_default);
	
	gsl_rng_set (rngen, -1*time(NULL)); // setting up the random number generator.


	n=net;
  
	StartInfectionDengue(net,params,rngen);
	UpdateInfectionStrainsDengue(net,params);

 
   for(i=0;i<iterations;i++){
		//if(i%1==0)printf("timestep %i\n",i);

	   inf1_1=inf1_2=inf1_3=inf1_4=inf2_1=inf2_2=inf2_3=inf2_4=0;

	   imm1_1=imm1_2=imm1_3=imm1_4=imm2_1=imm2_2=imm2_3=imm2_4=0;

	   coinf0=coinf1=coinf2=coinf3=0;

	   not_inf=not_imm=0;

	   n=net;
	   while(n->next!=NULL){	
		   n=n->next;
		   ctr=0;

		s=n->immune;
		if(s->next==NULL)not_imm++;
		while(s->next!=NULL){
			switch(s->next->string){
				case 0:
					imm1_1++;
					break;
				case 1:
					imm1_2++;
					break;
				case 2:
					imm1_3++;
					break;
				case 3:
					imm1_4++;
					break;
			}
			s=s->next;
		}

		s=n->infected;

		if(s->next==NULL)not_inf++;

		while(s->next!=NULL){
			switch(s->next->string){
				case 0:
					inf1_1++;
					break;
				case 1:
					inf1_2++;
					break;
				case 2:
					inf1_3++;
					break;
				case 3:
					inf1_4++;
					break;
			}
			s=s->next;
			ctr++;
		}
		if(ctr>0){
		   
		   //InfectStrainsNoStructureDengueStrainSpecificCrossImmunity(n,params,nodes_array,i,rngen);
		   InfectStrainsNoStructureDengueStrainSpecificCrossImmunityAndTransmissionEnhancement(n,params,nodes_array,i,rngen);
		   coinf0++;
		}

		if(ctr>1){
			coinf2++;
			printf("ERROR, more than two strains in one node\n");
		}
		if(ctr>2)
;
		if(ctr>3)coinf3++;

	   }
	//fprintf(coinf0f,"%d %d\n",i,coinf0);
	//fprintf(coinf1f,"%d %d\n",i,coinf1);
	//fprintf(coinf2f,"%d %d\n",i,coinf2);
	//fprintf(coinf3f,"%d %d\n",i,coinf3);

	   //fprintf(inf11,"%d %d\n",i,inf1_1);
	   //fprintf(inf12,"%d %d\n",i,inf1_2);
	   //fprintf(inf13,"%d %d\n",i,inf1_3);
	   //fprintf(inf14,"%d %d\n",i,inf1_4);

	//fprintf(inf21,"%d %d\n",i,inf2_1);
	//fprintf(inf22,"%d %d\n",i,inf2_2);
	//fprintf(inf23,"%d %d\n",i,inf2_3);
	//fprintf(inf24,"%d %d\n",i,inf2_4);

    //fprintf(imm11,"%d %d\n",i,imm1_1);
    //fprintf(imm12,"%d %d\n",i,imm1_2);
    //fprintf(imm13,"%d %d\n",i,imm1_3);
    //fprintf(imm14,"%d %d\n",i,imm1_4);

    //fprintf(imm21,"%d %d\n",i,imm2_1);
    //fprintf(imm22,"%d %d\n",i,imm2_2);
    //fprintf(imm23,"%d %d\n",i,imm2_3);
    //fprintf(imm24,"%d %d\n",i,imm2_4);

	   inf_data[1][i]=inf1_1;
	   inf_data[2][i]=inf1_2;
	   inf_data[3][i]=inf1_3;
	   inf_data[4][i]=inf1_4;

	   imm_data[1][i]=imm1_1;
	   imm_data[2][i]=imm1_2;
	   imm_data[3][i]=imm1_3;
	   imm_data[4][i]=imm1_4;
	   
	   fprintf(notimm,"%d %d\n",i,not_imm);
	fprintf(notinf,"%d %d\n",i,not_inf);

//	diversity1=MeasureDiversity(inf1_1,inf1_2,inf1_3,inf1_4);
//	discordance1=MeasureDiscordance(inf1_1,inf1_2,inf1_3,inf1_4);

//	diversity2=MeasureDiversity(inf2_1,inf2_2,inf2_3,inf2_4);
//	discordance2=MeasureDiscordance(inf2_1,inf2_2,inf2_3,inf2_4);

	//fprintf(div1f,"%d %lf\n",i,diversity1);
	//fprintf(disc1f,"%d %lf\n",i,discordance1);

	//fprintf(div2f,"%d %lf\n",i,diversity2);
	//fprintf(disc2f,"%d %lf\n",i,discordance2);


	inf1=inf1_1+inf2_1;
	inf2=inf1_2+inf2_2;
	inf3=inf1_3+inf2_3;
	inf4=inf1_4+inf2_4;

	if((inf1+inf2+inf3+inf4)==0){
		printf("infection died\n");

		fclose(notinf);
		fclose(notimm);
		ClearInfection(net);
		return ;
	}

	//diversity_tot=MeasureDiversity(inf1,inf2,inf3,inf4);
	//discordance_tot=MeasureDiscordance(inf1,inf2,inf3,inf4);


	
	RecoverInfectedDengue(net,params);
    //UpdateInfectionStrainsInstantImmunity(net,params);  

	UpdateInfectionStrainsDengue(net,params); 
	



  }
  ClearInfection(net);

}


void RecoverInfectedLoseImmunityStrains(struct node_gra *net,
					parameters *params){

  struct node_gra *n;

  strain *s_inf,*s_imm;
  strain *s_tmp;
  n=net;

  while(n->next!=NULL){  

    n=n->next;
       				/* lose immunity to a random strain */
    s_imm=n->immune;

	while(s_imm->next!=NULL){
		if(params->p_loss_imm>ran2(&seed)){
			s_tmp=s_imm->next;
			s_imm->next=s_imm->next->next;
			free(s_tmp);
		}
      
      if(s_imm->next!=NULL)s_imm=s_imm->next;
      
    }

    
    if(n->infected->next!=NULL){/* if node is infected, with prob.
		 p_recover, it recovers from and becomes immune to a strain*/
      
      s_inf=n->infected;
      
      while(s_inf->next!=NULL){

	if(params->p_recover>ran2(&seed)){
	  s_imm=n->immune;
	  while(s_imm->next!=NULL)s_imm=s_imm->next; /* go to end of list */

	  s_imm->next=s_inf->next;       /* put infected strain into immune */
	  s_inf->next=s_inf->next->next; /* remove from infected list */
	  s_imm->next->next=NULL;        /* end immune list */
	  n->immune->string=1;
	}
	
	if(s_inf->next!=NULL)s_inf=s_inf->next;

      }
      
      if(n->infected->next==NULL)n->infected->string=0;
    }
    
  }
}




void RecoverInfectedDengue(struct node_gra *net,
					parameters *params){

  struct node_gra *n;
  strain *s_inf,*s_imm;
  strain *s_tmp;
  n=net;

  while(n->next!=NULL){  

    n=n->next;
       				
	if(params->p_loss_imm>ran2(&seed)){
		if(n->immune->next!=NULL){
			n->immune->next=ClearStrains(n->immune);
		}
		if(n->infected->next!=NULL){
			n->infected->next=ClearStrains(n->infected);
		}
		if(n->queue->next!=NULL){
			n->queue->next=ClearStrains(n->queue);
		}
	}
	else if(n->infected->next!=NULL){/* if node is infected, with prob.
		//						p_recover, it recovers from and becomes immune to a strain*/

		s_inf=n->infected;

		while(s_inf->next!=NULL){

			if(params->p_recover>ran2(&seed)){
				s_imm=n->immune;
				while(s_imm->next!=NULL)s_imm=s_imm->next; /* go to end of list */

				s_imm->next=s_inf->next;       /* put infected strain into immune pile*/
				s_inf->next=s_inf->next->next; /* remove from infected list */
				s_imm->next->next=NULL;        /* end immune list */
//				n->immune->string=1;
			}

			if(s_inf->next!=NULL)s_inf=s_inf->next;
		}

		if(n->infected->next==NULL)n->infected->string=0;
	}

  }

}


double MeasureDiversity(int inf1,int inf2,int inf3,int inf4){

  double diversity;
  int num_strains=NUM_STRAINS;
  double *p;
  int totalinf=inf1+inf2+inf3+inf4;
  double sum;
  
  int i;
  
  p=(double *)calloc(num_strains,sizeof(double));

  p[0]=(double)inf1/totalinf;
  p[1]=(double)inf2/totalinf;
  p[2]=(double)inf3/totalinf;
  p[3]=(double)inf4/totalinf;
  sum=0;
  for(i=0;i<num_strains;i++){
    if(p[i]!=0)
      sum+=-p[i]*log10(p[i]);
  }
  diversity=sum/((double)NUM_LOCI*log10(2.0));

  return diversity;
}


double MeasureDiscordance(int inf1,int inf2,int inf3,int inf4){
  int i,j;
  int num_strains=(int)pow((float)2,(float)NUM_LOCI);
  double *p;
  int totalinf=inf1+inf2+inf3+inf4;
  double sum_top,sum_bottom;
  double discordance;

  p=(double *)calloc(num_strains,sizeof(double));
  p[0]=(double)inf1/totalinf;
  p[1]=(double)inf2/totalinf;
  p[2]=(double)inf3/totalinf;
  p[3]=(double)inf4/totalinf;
  sum_top=sum_bottom=0;

  for(i=0;i<num_strains;i++){
    for(j=i+1;j<num_strains;j++){
      sum_top+=((double)HammingDist(i,j)*p[i]*p[j]);
      sum_bottom+=((double)p[i]*p[j]);
    }
  }
  if(sum_bottom!=0)discordance=sum_top/(double)(NUM_LOCI*sum_bottom);
  else discordance=0;
  return discordance;
}

int HammingDist(int i,int j){
  int ii;
  int hamming;
  hamming=0;
  for(ii=0;ii<NUM_LOCI;ii++){
    if(BitTst(i,ii)!=BitTst(j,ii))hamming++;
  }
  return hamming;
}
strain *ClearStrains(strain *str){

	strain *temp=str->next;
	strain *tmp;
	if(temp==NULL)return NULL;
	while(temp->next!=NULL){
		tmp=temp;
		temp=temp->next;
		free(tmp);
	}
	free(temp);
	return NULL;
}

void ClearInfection(struct node_gra *net){
  
	struct node_gra *n=net;

	while(n->next!=NULL){
		n=n->next;
		ClearStrains(n->infected);
		free(n->infected);
		ClearStrains(n->immune);
		free(n->immune);
		ClearStrains(n->queue);
		free(n->queue);
	}
}


void PrintNodeStatus(FILE *trackinff,FILE *trackimmf,
		     int timestep,struct node_gra *n){
  strain *s;
  int i,j;

  fprintf(trackinff,"%d ",timestep);
  fprintf(trackimmf,"%d ",timestep);

  s=n->infected;
  
  i=0;
  while(s->next!=NULL){
    s=s->next;
    fprintf(trackinff,"%d ",s->string);
    i++;
  }
  for(j=i;j<4;j++)fprintf(trackinff,"%d ",-1);
  fprintf(trackinff,"\n");
  
  s=n->immune;
  i=0;
  while(s->next!=NULL){
    s=s->next;
    fprintf(trackimmf,"%d ",s->string);
    i++;
  }
  for(j=i;j<4;j++)fprintf(trackimmf,"%d ",-1);
  fprintf(trackimmf,"\n");
}

void PrintNodeStatus2(FILE *trackinf1f,FILE *trackinf2f,FILE *trackinf3f,
		      FILE *trackinf4f,
		      FILE *trackimm1f,FILE *trackimm2f,FILE *trackimm3f,
		      FILE *trackimm4f,
		      int timestep,struct node_gra *n){
  strain *s;
  int i,j;
  int ok1,ok2,ok3,ok4;

  fprintf(trackinf1f,"%d ",timestep);
  fprintf(trackinf2f,"%d ",timestep);
  fprintf(trackinf3f,"%d ",timestep);
  fprintf(trackinf4f,"%d ",timestep);
 
  fprintf(trackimm1f,"%d ",timestep);
  fprintf(trackimm2f,"%d ",timestep);  
  fprintf(trackimm3f,"%d ",timestep);  
  fprintf(trackimm4f,"%d ",timestep);
  
  s=n->infected;
  
  ok1=ok2=ok3=ok4=0;
  while(s->next!=NULL){
    s=s->next;
    switch(s->string){
    case 0:
      fprintf(trackinf1f,"%d\n",s->string);
      ok1=1;
      break;
    case 1:
      fprintf(trackinf2f,"%d\n",s->string);
      ok2=1;
      break;
    case 2:
      fprintf(trackinf3f,"%d\n",s->string);
      ok3=1;
      break;
    case 3:
      fprintf(trackinf4f,"%d\n",s->string);
      ok4=1;
      break;
    }
    
  }
  if(ok1==0)fprintf(trackinf1f,"-1\n");
  if(ok2==0)fprintf(trackinf2f,"-1\n");
  if(ok3==0)fprintf(trackinf3f,"-1\n");
  if(ok4==0)fprintf(trackinf4f,"-1\n");

  s=n->immune;
  ok1=ok2=ok3=ok4=0;
  while(s->next!=NULL){
    s=s->next;
    switch(s->string){
    case 0:
      fprintf(trackimm1f,"%d\n",s->string);
      ok1=1;
      break;
    case 1:
      fprintf(trackimm2f,"%d\n",s->string);
      ok2=1;
      break;
    case 2:
      fprintf(trackimm3f,"%d\n",s->string);
      ok3=1;
      break;
    case 3:
      fprintf(trackimm4f,"%d\n",s->string);
      ok4=1;
      break;
    }
    
  }
  
  if(ok1==0)fprintf(trackimm1f,"-1\n");
  if(ok2==0)fprintf(trackimm2f,"-1\n");
  if(ok3==0)fprintf(trackimm3f,"-1\n");
  if(ok4==0)fprintf(trackimm4f,"-1\n");
  
}


void PrintStrainClusterFile(struct node_gra *net){
 
  FILE *cluinff=fopen("infected1.clu","w");

  FILE *cluimmf=fopen("immune1.clu","w");
  
  struct node_gra *n;
  strain *s;

  n=net;
  while(n->next!=NULL){
    n=n->next;
    s=n->infected;
    if(s->next==NULL)fprintf(cluinff,"%d\n",0);
    while(s->next!=NULL){
      switch(s->next->string){
      case 0:
	fprintf(cluinff,"%d\n",1);
	break;
      case 1:
	fprintf(cluinff,"%d\n",2);
	break;
      case 2:
	fprintf(cluinff,"%d\n",3);
	break;
      case 3:
	fprintf(cluinff,"%d\n",4);
	break;
      }
      s=s->next;
      break;
    }  
  
    s=n->immune;
    if(s->next==NULL)fprintf(cluinff,"%d\n",0);
    while(s->next!=NULL){
      switch(s->next->string){
      case 0:
	fprintf(cluimmf,"%d\n",1);
	break;
      case 1:
	fprintf(cluimmf,"%d\n",2);
	break;
      case 2:
	fprintf(cluimmf,"%d\n",3);
	break;
      case 3:
	fprintf(cluimmf,"%d\n",4);
	break;
      }
      
      s=s->next;
      break;
    }
  }
  system("unix2dos immune1.clu");
  system("unix2dos infected1.clu");
}
