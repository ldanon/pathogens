#ifndef _PATHOGENS_H
#define _PATHOGENS_H
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "globals.h"

typedef struct strain{
  int string;
  int timer;
  int length;
  struct strain *next;
  struct node_gra *host;
}strain;

typedef struct parameters{
  double p_infect;
  double p_recover;
  double p_loss_imm;
  double p_init;
  double p_mutate;
  double p_recom;
  double gamma;
  double vulnerability[10];
  int susc_times[NUM_STRAINS][NUM_STRAINS]; // the time for the node to be susceptible again is a function of 
						  // history and what's trying to infect now. 
  double p_out; /* probability of having a link with other com*/
  double p_in; /* probability of having a link inside the same com*/
  double k,k2,theta; // parameters of gamma function, used for setting length of cross protection
  double enhancement;
}parameters;

				/* Functions */
void RunSIRSOnNetwork(struct node_gra *net,parameters *params,int iterations);

parameters *InitialiseParameters(double p_init,double p_infect, 
				 double p_recover,double p_mutate,
				 double p_loss_imm,double p_recom,
				 double gamma,
				 double k, double theta,double k2);

double CalculateAverage(char infname[50],int cutoff);

void UpdateVulnerabilityVector(parameters *par);
strain *CreateStrain(int num);
strain *CreateStrainTimestamp(int num, int T);
void RemoveStrain(strain *str);

void StartRandomInfectionSIRS(struct node_gra *net,double p_init);

void StartInfectionInCom1SIRS(struct node_gra *net,double p_init);

void InfectNeighboursSIRS(struct node_gra *node, double prob);

void InfectStrainsNoStructureDengue(struct node_gra *n, parameters *params, 
									struct node_gra *nodes_array[N_NODE],int T,gsl_rng *rngen);

void InfectStrainsNoStructureDengueStrainSpecificCrossImmunity(struct node_gra *n, parameters *params, 
									struct node_gra *nodes_array[N_NODE],int T,gsl_rng *rngen);
									
void InfectStrainsNoStructureDengueStrainSpecificCrossImmunityAndTransmissionEnhancement(struct node_gra *n, parameters *params, 
									struct node_gra *nodes_array[N_NODE],int T,gsl_rng *rngen);									

void RecoverInfectedLoseImmunitySIRS(struct node_gra *net,parameters *params);

void UpdateInfectionSIRS(struct node_gra *net,parameters *params);

void StartInfectionStrains(struct node_gra *net,
			   double p_init,int start_strain);

void StartInfection2Strains(struct node_gra *net,
			    double p_init,int start_strain1,int start_strain2);

void StartInfectionDengue(struct node_gra *net,
							parameters *params, gsl_rng *rngen);

void MutateStrain(strain *s,parameters *params);

int RecombineStrains(struct node_gra *n,parameters *params);

int InfectedWithStrain(struct node_gra *n,int num);
int InfectedWithAnything(struct node_gra *n);
int IsSomethingInQueue(struct node_gra *n);
int ImmuneToStrain(struct node_gra *n,int num);

void UpdateInfectionStrains(struct node_gra *net, parameters *params);
void UpdateInfectionStrainsInstantImmunity(struct node_gra *net, 
					   parameters *params);
void UpdateInfectionStrainsDengue(struct node_gra *net,parameters *params);


void InfectNeighboursStrains(struct node_gra *node,parameters *params);

double CalculateVulnerabilityOfNode(struct node_gra *n,
				    int s1,parameters *params);
int CalculateVulnerabilityOfNodeDengue(struct node_gra *n,
										  int s1,parameters *params, int T);
double CalculateVulnerabilityOfNodeTest(struct node_gra *n,
				    int s1s,parameters *params);

void RunStrainsOnNetwork(struct node_gra *net,parameters *params,
			 int iterations);
void RunStrainsOnComsNetwork(struct node_gra *net,parameters *params,
			 int iterations);
void RunStrainsOnMetapop(struct node_gra *net,parameters *params,
						int iterations,struct node_gra *node_array[N_NODE]);

void RunStrainsOnMetapopDengue(struct node_gra *net,parameters *params,
						  int iterations, struct node_gra *nodes_array[N_NODE],
						  int inf_data[NUM_STRAINS+1][ITERATIONS],int imm_data[NUM_STRAINS+1][ITERATIONS]);

void RecoverInfectedLoseImmunityStrains(struct node_gra *net,
					parameters *params);

void RecoverInfectedDengue(struct node_gra *net,
					parameters *params);
double MeasureDiversity(int inf1,int inf2,int inf3,int inf4);

double MeasureDiscordance(int inf1,int inf2,int inf3,int inf4);

int HammingDist(int i,int j);

void ClearInfection(struct node_gra *net);
strain *ClearStrains(strain *str);

void PrintNodeStatus(FILE *trackinff,FILE *trackimm,
		     int timestep,struct node_gra *n);
void PrintNodeStatus2(FILE *trackinf1f,FILE *trackinf2f,FILE *trackinf3f,
		      FILE *trackinf4f,
		      FILE *trackimm1f,FILE *trackimm2f,FILE *trackimm3f,
		      FILE *trackimm4f,
		      int timestep,struct node_gra *n);





#endif