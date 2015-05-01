#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>

#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_statistics_double.h>   // gnu scientific library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>


#include "rand_gen.h"
#include "network_lib.h"
#include "quicksort.h"
#include "pathogens.h"
#include "globals.h"
#include "bitops.cpp"

void MeasureTimeOfExtinction(int extinctions[], int inf_data[NUM_STRAINS+1][ITERATIONS]);
void OutputRunToFile(int inf_data[NUM_STRAINS+1][ITERATIONS],const char *fname,int format);

int main (int argc, const char *argv[])
{

	// run program like this: str network name (double)R_o  (double)gamma 

	int i,r;
	struct node_gra *net,*n,*nodes_array[N_NODE];


	parameters *params;
	double p_init,p_infect,p_recover,p_mutate,p_loss_imm;
	double gamma;
	double p_recom;
	double R_o;
	double k,k2,theta;
	double enhancement;

	int inf_data[NUM_STRAINS+1][ITERATIONS];
	int imm_data[NUM_STRAINS+1][ITERATIONS];
	
	int extinct[NUM_STRAINS+1];

	FILE *outFile;


	//	int nnod;
	//	float Z_tot;
	//	float Z_out,Z_out_start;

	//	double beta;

	//double sum,dtmp;
	//int itmp;
	//int cutoff=512; //when calculating average, 
	//the number of steps to discard at beginning of run

	seed=-1*time(NULL);
	outFile=fopen("k_v_extinctions.dat","w");


	//	R_o=2;
	// A time step represents one day.

	p_init=(double)100/(N_NODE);			/* initial infection prob, how many people you seed infection with*/
	p_infect=0.0;	                /* transmissibility beta*/
	p_recover=0.1;			/* recovery rate sigma*/
	p_loss_imm=(double) 1/(1*100);	/* loss of immunity = 1/average length span = 1/50years*/
	p_mutate=0.001;		/* mutation rate*/
	p_recom=0.001;			/* not set here */
	gamma=1;				/* not set here I HAVE TO CHANGE THE WAY THIS WORKS*/
	k=1000;					/* mean of the gamma function */
	theta=0.5;				/* sharpness of the gamma function */
	k2=400;
	enhancement=2;			/* level of Immune Dependent Enhancement */

	/* R_o=beta / p_recover>1 */
	/* beta=p_infect*<k> */
	/* =>R_o=p_infect*<k>/p_recover */
	/* =>p_infect=R_o*p_recover/<k> */



	//printf("argc %d\n",argc);
	if(argc>1){
		switch(argc){
			case 1:
				break;
			case 2:
				R_o=atof(argv[1]);
				printf("You better be specifying R0 %lf\n",R_o);
				//				getchar();
				break;
			case 3:
				R_o=atof(argv[1]);
				p_recover=atof(argv[2]);
				printf("You better be specifying R0 %lf \n and recovery rate %lf\n",R_o,p_recover);
				//				getchar();

				break;
			case 4:
				R_o=atof(argv[1]);
				p_recover=atof(argv[2]);
				k=atoi(argv[3]);

				printf("You better be specifying R0 %lf \n and recovery rate %lf \n and length of cross immunity %lf \n",R_o,p_recover,k);
				//				getchar();

				break;
			case 5:
				R_o=atof(argv[1]);
				p_recover=atof(argv[2]);
				k=atoi(argv[3]);
				gamma=atof(argv[4]);

				printf("You better be specifying R0 %lf \n and recovery rate  %lf \n and length of cross immunity %lf \n  and gamma %lf\n",R_o,p_recover,k,gamma);
				//				getchar();

				break;
			case 6:
				R_o=atof(argv[1]);
				p_recover=atof(argv[2]);
				k=atoi(argv[3]);
				gamma=atof(argv[4]);
				enhancement=atof(argv[5]);

				printf("You better be specifying R0 %lf \n and recovery rate   %lf \n and length of cross immunity %lf \n and gamma %lf \n and enhancement %lf\n",R_o,p_recover,k,gamma,enhancement);
				//				getchar();
				break;
		}
	}
	//for(k=0;k<400;k+=10){
	for(k=199;k<201;k++){
		params=InitialiseParameters(p_init,p_infect,p_recover,
			p_mutate,p_loss_imm,p_recom,gamma,k,theta,k2);


		params->p_infect  =  R_o * (params->p_recover + params->p_loss_imm)/(double)N_NODE;

		if(params->p_infect >= 1.0){
			printf("Make Derek Happy :) params::p_infect = %lf\n",params->p_infect);
			getchar();
			return 0;
		}

		params->enhancement=enhancement;

		printf("p_infect %lf recovery rate %lf\n", params->p_infect,params->p_recover);

		net=BuildNetNoLinksMetapops(nodes_array);
		
		n=net;
		
		printf("k %lf k2 %lf\n",params->k,params->k2);
		params->k2=params->k;
		for(i=0;i<=NUM_STRAINS;i++)extinct[i]=0;

		for(r=0;r<REALS;r++){
			printf("REALISATION %d\n",r+1);

			//printf("running %d\n",ITERATIONS);
			RunStrainsOnMetapopDengue(net,params,ITERATIONS,nodes_array,inf_data,imm_data);
			OutputRunToFile(inf_data,"Infecteds.dat",2);
			printf("done\n");
			MeasureTimeOfExtinction(extinct,inf_data);
			for(i=1;i<=NUM_STRAINS;i++)printf(" %lf",(double)extinct[i]);
		}
		fprintf(outFile,"%lf ",k);
		for(i=1;i<=NUM_STRAINS;i++)fprintf(outFile," %lf",(double)extinct[i]/REALS);
		//printf("\n");	
		fprintf(outFile,"\n");
		fflush(outFile);
		
	}	
	fclose(outFile);
	RemoveGraph(net);
	//_CrtDumpMemoryLeaks();
	getchar();

	return 0;
}
void OutputRunToFile(int inf_data[NUM_STRAINS+1][ITERATIONS],const char *fname,int format){
	// format = 1 means xmgrace output
	// format = 2 mean matlab output
	
	FILE *outF=fopen(fname,"w");

	int i,j;
	if(format==1){
		for(i=1;i<NUM_STRAINS+1;i++){

			for(j=0;j<ITERATIONS;j++){

				fprintf(outF,"%d %d\n",j,inf_data[i][j]);

			}

			fprintf(outF,"\n");

		}
	}
	else {

		for(j=0;j<ITERATIONS;j++){
		//	fprintf(outF,"%d",j);
			for(i=1;i<NUM_STRAINS+1;i++){

				fprintf(outF," %d",inf_data[i][j]);
			
			}
			fprintf(outF,"\n");
	
		}	
	}

	fclose(outF);

}
void MeasureTimeOfExtinction(int extinctions[], int inf_data[NUM_STRAINS+1][ITERATIONS]){
	

	int i;
	int i1,t1=0;


	for(i=0;i<=NUM_STRAINS;i++){
		
		i1=1;
		while(inf_data[i][t1]>0){
			t1++;
			if(i1==0){
				extinctions[i]=t1;
			}
		}
	}	
	return;
}
