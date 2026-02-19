#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_coupling.h>
#include "mpi.h"

#define MAX(a,b)	((a)>(b))?(a):(b)
#define MIN(a,b)	((a)<(b))?(a):(b)
#define cube(x)	((x)*(x)*(x))
#define sqr(x)		((x)*(x))
#define LINESIZE 1024

#include "genFunctions/genStructures.h"
#include "genFunctions/genUtils.h"
#include "genFunctions/genSurfUtils.h"
#include "genFunctions/genEnCalcs.h"
#include "genFunctions/genConstraintCalcs.h"
#include "genFunctions/genMCMove.h"
#include "genFunctions/genParallelTempering.h"
#include "genFunctions/genFileInput.h"



	int main (int argc, char *argv[])
	{
		/*************************************************** MPI VARIABLES START **************************************************/

	int rank, size, domrank, received, notreceived,checknTot;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	MPI_Request request;
	MPI_Group groupall, group;

	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_group( MPI_COMM_WORLD, &groupall);
	printf("created MPI conditions\n");

/*************************************************** MPI VARIABLES END **************************************************/

	
		FILE *fpPara,*fpInteract,*fpChromProp,*fpIPD,*fpWrite,*fpTmp;

	char buffer[800],buffer2[800],buffer3[800],spambuff[100];
	double W,Wold=0.0,elastEn=0.0,oldelastEn=0.0;

	double random=0.0,trackerspam=0.0;
	int i,j,ii,jj,indi,indj;
	int ii2,jj2,ii3,jj3,ii4,jj4;



	printf("creating structs\n");

	struct SYSTEM structbox;
	struct SYSTEM tmpbox;
	struct PARAINPUT p;

		printf("args %s %s %s %s\n",argv[1],argv[2],argv[3],argv[4]);

	sprintf(buffer,"%s",argv[1]); // parafile
	printf("%s\n",buffer);
	fpPara=fopen(buffer,"r");
	if(fpPara==NULL)
	{
		printf("no parafile\n");
		exit(0);
	}
	readPara(fpPara,&p,&structbox);
	fclose(fpPara);
	p.nDChrom = 2*p.nChrom; // diploidy
	p.nAto = 2*(p.nChrom + p.nPatchTot) + p.nDChrom; // 2xnchrom + 2x(npatches-per-chrom * nchrom) + 2x laminpatches-per-chrom*nchrom
	printf("read para, simulating %d particles (well here only %d chrom centres)\n",p.nAto,p.nDChrom);
		char **fileList,**IADfileList,**ChromPropsfileList,**IPDfileList;
		double **source;
		// Read source rows: T eps seed outpath IAD ChromProps IPD.
		source=calloc(size,sizeof(double *));
		fileList=calloc(size,sizeof(char *));
		IADfileList=calloc(size,sizeof(char *));
		ChromPropsfileList=calloc(size,sizeof(char *));
		IPDfileList=calloc(size,sizeof(char *));
		for(i=0;i<size;i++)
		{
			source[i]=calloc(3,sizeof(double *)); // T eps_hertz seed
			fileList[i]=calloc(500,sizeof(char));
			IADfileList[i]=calloc(500,sizeof(char));
			ChromPropsfileList[i]=calloc(500,sizeof(char));
			IPDfileList[i]=calloc(500,sizeof(char));

		}
	//printf("Start reading source file %s...\n",argv[2]);
	fpPara=fopen(argv[2],"r");		// Source file (source.dat)
	if(fpPara==NULL)
	{
		printf("Error: %s not found!\n",argv[2]);
		exit(1);
	}
		readSource(fpPara,source,fileList,IADfileList,ChromPropsfileList,IPDfileList,size);
	fclose(fpPara);
	printf("Done reading source file %s!\n",argv[2]);

	MPI_Barrier(MPI_COMM_WORLD);

	p.mctemp = source[rank][0];
	p.eps_hertz = source[rank][1];
	p.seed=(int)(source[rank][2]);
	// p.IPDfile etc assign
	sprintf(p.mcfile,"%s",fileList[rank]);
		sprintf(p.IADfile,"%s",IADfileList[rank]);
		sprintf(p.ChromPropsfile,"%s",ChromPropsfileList[rank]);
		sprintf(p.IPDfile,"%s",IPDfileList[rank]);

		// Allocate and initialise arrays for this rank.
		initAll(&p,&structbox,&tmpbox);
		countList(&p,source,size,rank);

	printf("read para and initialised arrays\n");

	MPI_Barrier(MPI_COMM_WORLD);
		double spam = 0.33333*(pow(p.ella*p.ellb,1.6) + pow(p.ella*p.ellc,1.6) + pow(p.ellb*p.ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p.surfPhi*4.0*surf/(M_PI*p.nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	double countacc, countrej=0.0,rand2=0.0;

	int l,t,mcstep;
	double mspam = 1.0; //pow(p.massscalefact/8.0,0.3333);


	// setting simulations settings

	srand48(p.seed);



		// Read chromosome properties (npatches/mass/type/fractions).
		fpChromProp=fopen(p.ChromPropsfile,"r");
		if(fpChromProp==NULL)
		{
			printf("chrom prop null pointer %s\n",p.ChromPropsfile);
			exit(0);
		}
		readChromInfo(fpChromProp,&p,&structbox);
	fclose(fpChromProp);


	printf("read info for %d chrom\n",p.nChrom);
	// read the interaction matrix
	fpInteract=fopen(p.IADfile,"r");
	if(fpInteract==NULL)
	{
		printf("interaction file null pointer %s\n",argv[2]);
		exit(0);
	}
	rewind(fpInteract);
	printf("reading interactions from %s\n",argv[2]);
	int npairs = readIntMatrix(fpInteract,&p,&structbox);
	printf("%d pairwise interactions doubled up on %d chrom\n",npairs,p.nChrom);
	fclose(fpInteract);

		// Read inter-chromosome IPD matrix.
		fpIPD=fopen(p.IPDfile,"r");
	if(fpIPD==NULL)
	{
		printf("IPD file null pointer %s\n",p.IPDfile);
		exit(0);
	}
		readIPDInfo(fpIPD,&p,&structbox);
		printf("read inter chrom centroid distances\n");


		// constraint mat, and constraint calc
	double spamconst=0.0;
	printf("placed surface\n");
	if(p.restartPt==0)
	{
		printf("restarting from %d means starting anew\n",p.restartPt);
		placeChromMCInit(&p,&structbox);
	}
	else
	{
		printf("restarting from conf at %d\n",p.restartPt);
		printf("NOTE!!!! Restarting does not preserve the histogram of l1 info\n");
		sprintf(buffer,"%s-dump_%d.dat",p.mcfile,p.restartPt);
		structbox.fpRestart=fopen(buffer,"r");
		if(structbox.fpRestart==NULL)
		{
			printf("did not find restart file %s\n",buffer);
			exit(0);
		}
		readChromPos(&p,&structbox);
		fclose(structbox.fpRestart);
	}


	for(ii=0;ii<p.nChrom;ii++)
	{
		for(jj=0;jj<p.nChrom;jj++)
		{
				if(ii!=jj)
				{
					structbox.expectedMat[ii][jj]=structbox.ipd[ii][jj];
					structbox.expectedMat[jj][ii]=structbox.ipd[ii][jj];

				if(p.fipd == 1)
				{
					structbox.toCalc[ii][jj] = 1;
					structbox.toCalc[jj][ii] = 1;
					if(ii==p.nChrom-1 | jj == p.nChrom-1)
					{
					structbox.toCalc[ii][jj] = 0;
					structbox.toCalc[jj][ii] = 0;
					}
					}

				}
				else
				{
					structbox.expectedMat[ii][jj]=0.0;
				}
				printf("exp mat %d %d %d %f\n",ii,jj,structbox.toCalc[ii][jj],structbox.expectedMat[ii][jj]);

			}

			// Keep centroid slot present for file-shape compatibility but inactive.
			structbox.expectedMat[ii][p.nChrom] = 0.0;
			structbox.expectedMat[p.nChrom][ii] = 0.0;
			structbox.toCalc[ii][p.nChrom] = 0;
			structbox.toCalc[p.nChrom][ii] = 0;

		}

	structbox.expectedMat[p.nChrom][p.nChrom]=0;
// 	exit(0);

	genSurf(&p,&structbox);




	// generate thermo, constraint and histo files

	if(p.restartPt ==0)
	{
		sprintf(buffer,"%s-thermo.dat",p.mcfile);
		printf("%s thermo file\n",buffer);
		structbox.fpmcthermo = fopen(buffer,"w");

		sprintf(buffer,"%s-constraints.dat",p.mcfile);
		printf("%s constraint file\n",buffer);
		structbox.fpconstraint = fopen(buffer,"w");
	}
	else
	{
		// open the thermo and constraint files -- read and write until restartpt
		// make a tmp dir -- write to tmp dir
		// once written, copy back
		sprintf(buffer,"%s-thermo.dat",p.mcfile);
		printf("%s thermo file\n",buffer);
// 		structbox.fpmcthermo = fopen(buffer,"r");
// 		rewind(structbox.fpmcthermo);
		sprintf(buffer2,"mkdir -p %s-spam/",p.mcfile);
		system(buffer2);

		sprintf(buffer2,"%s-spam/tmpThermo.dat",p.mcfile);
		sprintf(buffer3,"cp %s-thermo.dat %s",p.mcfile,buffer2);
		system(buffer3);
// 		fpTmp = fopen(buffer2,"w");
// 		while(!(feof(structbox.fpmcthermo)))
// 		{
// 			// read a line, write a line
// 			fscanf(structbox.fpmcthermo,"%d %lf %lf %lf %lf %lf %d %d %d %lf %lf\n",&mcstep,&structbox.currW,&structbox.curren,&structbox.cosine_sim,&structbox.l_infty,&structbox.l1,&structbox.maxdev,&structbox.maxi,&structbox.maxj,&structbox.devdir,&structbox.proxLinf);
//
// 			if(mcstep < p.restartPt)
// 			{
// 			fprintf(fpTmp,"%d %f %f %f %f %f %d %d %d %f %f\n",mcstep,structbox.currW,structbox.curren,structbox.cosine_sim,structbox.l_infty,structbox.l1,structbox.maxdev,structbox.maxi,structbox.maxj,structbox.devdir,structbox.proxLinf);
// 			fflush(structbox.fpmcthermo);
// 			}
// 		}
// 		fclose(structbox.fpmcthermo);
// 		fclose(fpTmp);

		// now over-write the old thermo with new one
		sprintf(buffer3,"cp %s %s-thermo.dat",buffer2,p.mcfile);
		system(buffer3);
		sprintf(buffer,"%s-thermo.dat",p.mcfile);
		printf("%s reopened thermo file\n",buffer);
		structbox.fpmcthermo = fopen(buffer,"a");



// 		sprintf(buffer,"%s-constraints.dat",p.mcfile);
// 		printf("%s constraint file\n",buffer);
// 		structbox.fpconstraint = fopen(buffer,"r");
// 		rewind(structbox.fpconstraint);
		sprintf(buffer2,"%s-spam/tmpConstraint.dat",p.mcfile);
// 		fpTmp = fopen(buffer2,"w");
//
// 		while(!(feof(structbox.fpconstraint)))
// 		{
// 			// read a line, write a line
//
// 			fscanf(structbox.fpconstraint,"%lf ",&spamconst);
// 			for(i=0;i<(p.nChrom+1)*(p.nChrom+1);i++)
// 			{
// 				indi = (int) i/(p.nChrom+1);
// 				indj = i % (p.nChrom+1);
// 				if(indj >= indi)
// 				{
// 					fscanf(structbox.fpconstraint,"%lf ",&structbox.constraintMat[indi][indj]);
//
// 				}
// 			}
// 			fscanf(structbox.fpconstraint,"%lf ",&structbox.cosine_sim);
//
// 			if(spamconst < p.restartPt)
// 			{
// 				fprintf(fpTmp,"%g ",spamconst);
// 				for(i=0;i<(p.nChrom+1)*(p.nChrom+1);i++)
// 				{
// 					indi = (int) i/(p.nChrom+1);
// 					indj = i % (p.nChrom+1);
// 					if(indj >= indi)
// 					{
// 						fprintf(fpTmp,"%g ",structbox.constraintMat[indi][indj]);
//
// 					}
// 				}
// 				fprintf(fpTmp,"%g\n",structbox.cosine_sim);
// 				fflush(fpTmp);
// 			}
// 		}
// 		fclose(structbox.fpconstraint);
// 		fclose(fpTmp);


		// now over-write the old thermo with new one
		sprintf(buffer3,"cp %s %s-constraints.dat",buffer2,p.mcfile);
		system(buffer3);
		sprintf(buffer,"%s-constraints.dat",p.mcfile);
		printf("%s reopened constraints file\n",buffer);
		structbox.fpconstraint = fopen(buffer,"a");


	}


	MPI_Barrier(MPI_COMM_WORLD);
	structbox.curren = fullEnergy(&p,&structbox);
	printf("init with old en %f\n",structbox.curren);
	structbox.currW = 0.0;
	structbox.cumuldev = 0.0;
	structbox.cumuldenom = 0.0;
// 	structbox.currW += IPDBias(&p,&structbox,0);

	/// set the initial W and en
	// below we only care about function calls and updating the constraint matrices so we lose the W info almost immediately

	// calcBias makes the function calls -- updates the l_infty and calls surfEnergy
	calcBias(&p,&structbox);
// 	structbox.currW = IPDBias(&p,&structbox,1);
// 	calcDeviations(&p,&structbox,0);
//
// 	structbox.currW = p.k_ipd*structbox.l_infty;

// 	for(ii=0;ii<p.nDChrom;ii++)
// 	{
// 		structbox.currW += surfEnergy(&p,&structbox,ii);
// 	}


	MPI_Barrier(MPI_COMM_WORLD);
// 	sprintf(buffer,"%s-histo.dat",p.mcfile);
// 	printf("%s thermo file\n",buffer);
// 	structbox.fpmcthermo = fopen(buffer,"w");


	// check annealing rate consistency

// 	if(p.mctemp - (1.0*p.MCIter/p.anneal_freq)*p.dtemp <=0.0)
// 	{
// 		printf("annealing gap too low or annealing delta temp too high -- EXITING\n");
// 		exit(0);
// 	}


// 	double mspam = pow(p.massscalefact/4.0,0.3333);
	double inittemp = p.mctemp;
	double deltemp = pow((p.dtemp/inittemp),(1.0*p.anneal_freq)/(1.0*p.MCIter));
// 	if(p.restartPt==0)
// 	{
// 	for(mcstep=p.restartPt;mcstep<p.restartPt+1000;mcstep++)
// 	{
// 		ChromMCStepElast(&p,&structbox);
// 	}
// 	}
	for(mcstep=p.restartPt;mcstep<p.MCIter;mcstep++)
	{

		if(p.gradMove==0)
		{
// 			printf("MC moving as usual -- %d %d\n",rank,mcstep);
			ChromMCStepLinfConstraint(&p,&structbox); // movement based on constraints

			// recompute elastic energy to track

			random=drand48();
			if(random>p.constraint_freq) // do an elastic relaxation every 10th step
			{

				ChromMCStepElast(&p,&structbox);
// 				if(random>0.99)
// 					ChromSwapStep(&p,&structbox);
			}
			calcBias(&p,&structbox);
		}

		if(p.gradMove==1) // update all chrom according to gradient
		{

			ChromGradStep(&p,&structbox,0);
			// configuration updated, recalculated constraint matrix

			calcBias(&p,&structbox);
// 			structbox.currW = IPDBias(&p,&structbox,1);
// 			calcDeviations(&p,&structbox,mcstep);
//
// 			structbox.currW = p.k_ipd*structbox.l_infty;
// // 			printf("All grad -- rank %d step %d pre surf W %f -- %f\n",rank,mcstep,structbox.currW,structbox.l_infty);
//
// 			for(i=0;i<p.nDChrom;i++)
// 			{
// 				structbox.currW += surfEnergy(&p,&structbox,i);
// 			}

			// elastic energy makes sure we don't exit bounds etc
			random=drand48();
			if(random>p.constraint_freq) // do an elastic relaxation every 10th step
			{

				ChromMCStepElast(&p,&structbox);
			}
		}
		if(p.gradMove==2) // update chromosome positions only according to largest gradient
		{
			random=drand48();
			if(random>p.constraint_freq) // && mcstep > 0.1*p.MCIter) // do an elastic relaxation every 10th 
			{

				ChromMCStepElast(&p,&structbox);

// 				if(random>0.99)
// 					ChromSwapStep(&p,&structbox);
//
				calcBias(&p,&structbox);
			}

			ChromGradStep(&p,&structbox,1);
			// configuration updated, recalculated constraint matrix
			calcBias(&p,&structbox);
// 			structbox.currW = IPDBias(&p,&structbox,1);
// 			calcDeviations(&p,&structbox,mcstep);
//
// 			structbox.currW = p.k_ipd*structbox.l_infty;
// // 			printf("Max grad -- rank %d step %d pre surf W %f -- %f\n",rank,mcstep,structbox.currW,structbox.l_infty);
//
// 			for(i=0;i<p.nDChrom;i++)
// 			{
// 				structbox.currW += surfEnergy(&p,&structbox,i);
// 			}
			// elastic energy makes sure we don't exit bounds etc

		}


		// parallel tempering swaps (T and elastEn)

		if(mcstep%p.Tswap==0 && mcstep>100000)
		{
// 			printf("entering T swap routine, %d %d on %d\n",mcstep,p.Tswap,rank);
			if(rank==0)		// Only master node does 'coin flip'
			{
				if(drand48()<0.5)
				{
					structbox.whoSwap=0;
				}
				else
				{
					structbox.whoSwap=1;
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&structbox.whoSwap,1,MPI_INT,0,MPI_COMM_WORLD);	// Master tells the value of whoSwap
			MPI_Barrier(MPI_COMM_WORLD);

			if(structbox.whoSwap==0)				// Attempt swap 0-1, 2-3, ...
			{
// 				printf("T swap whoswap 0\n");
				for(l=0;l<p.neps;l++)
				{
				for(t=0;t<p.nT-1;t+=2)
				{
					if(rank==p.procList[l][t] || rank==p.procList[l][t+1])		// Go into with 0-1, 2-3, ...
					{
						if(rank==p.procList[l][t])		// Proc 0,2,4,... send
						{
							TswapSend(&structbox,p.procList[l][t+1]);
						}
						else							// Proc 1,3,5,... receive and send response
						{
							TswapReceive(&structbox,p.procList[l][t],source,&p);
						}
						if(structbox.yestoSwap==1)	// Swap accepted
						{
							if(rank==p.procList[l][t])		// Proc 0,2,4,... send tempbox
							{
								SendTmpBox(&structbox,&tmpbox,&p,p.procList[l][t+1]);
							}
							else							// Proc 1,3,5,... received StructBox
							{
								ReceiveTmpBox(&structbox,&tmpbox,&p,p.procList[l][t]);
							}

							structbox.curren = fullEnergy(&p,&structbox);
							printf("new en (T swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.curren);
							structbox.currW = 0.0;
							structbox.cumuldev = 0.0;
							structbox.cumuldenom = 0.0;

							// configuration updated, recalculated constraint matrix
							calcBias(&p,&structbox);
// 							structbox.currW = IPDBias(&p,&structbox,1);
// 							calcDeviations(&p,&structbox,mcstep);
//
// 							structbox.currW = p.k_ipd*structbox.l_infty;
//
// 							for(ii=0;ii<p.nDChrom;ii++)
// 							{
// 								structbox.currW += surfEnergy(&p,&structbox,ii);
// 							}


							printf("new en (T swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.currW);

							if(rank==p.procList[l][t])
							{
								printf("whoswap %d eps_hertz: %g - %g: ",structbox.whoSwap,source[ p.procList[l][t] ][1],source[ p.procList[l][t+1] ][1]);
								printf("MC %d: Swap in T completed between %g e %g -- ranks are %d %d\n",t,source[ p.procList[l][t] ][0],source[ p.procList[l][t+1] ][0],p.procList[l][t],p.procList[l][t+1]);
							}
						}

					}
				}
				}
			}
			else if(structbox.whoSwap==1 && p.nT>2)			// Attempt swap 1-2, 3-4, ...
			{

				for(l=0;l<p.neps;l++)
				{
				for(t=1;t<p.nT-1;t+=2)
				{
					if(rank==p.procList[l][t] || rank==p.procList[l][t+1])		// Go into with 1-2, 3-4, ...
					{
						if(rank==p.procList[l][t])		// Proc 1,3,5,... send
						{
							TswapSend(&structbox,p.procList[l][t+1]);
						}
						else							// Proc 2,4,6,... receive and send response
						{
							TswapReceive(&structbox,p.procList[l][t],source,&p);
						}
						if(structbox.yestoSwap==1)	// Swap accepted
						{
							if(rank==p.procList[l][t])		// Proc 1,3,5,... send tempBox
							{
								SendTmpBox(&structbox,&tmpbox,&p,p.procList[l][t+1]);
							}
							else							// Proc 2,4,6,... receive StructBox
							{
								ReceiveTmpBox(&structbox,&tmpbox,&p,p.procList[l][t]);
							}
							structbox.curren = fullEnergy(&p,&structbox);
							printf("new en (T swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.curren);
							structbox.currW = 0.0;
							structbox.cumuldev = 0.0;
							structbox.cumuldenom = 0.0;
							// configuration updated, recalculated constraint matrix

							calcBias(&p,&structbox);
// 							structbox.currW = IPDBias(&p,&structbox,1);
// 							calcDeviations(&p,&structbox,mcstep);
//
// 							structbox.currW = p.k_ipd*structbox.l_infty;
//
// 							for(ii=0;ii<p.nDChrom;ii++)
// 							{
// 								structbox.currW += surfEnergy(&p,&structbox,ii);
// 							}


							printf("new en (T swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.currW);

							if(rank==p.procList[l][t])
							{
								printf("whoswap %deps: %g - %g: ",structbox.whoSwap,source[ p.procList[l][t] ][1],source[ p.procList[l][t+1] ][1]);
								printf("MC %d: other Swap in T completed between %g e %g -- ranks are %d %d\n",t,source[ p.procList[l][t] ][0],source[ p.procList[l][t+1] ][0],p.procList[l][t],p.procList[l][t+1]);
							}

						}
					}
				}
				}
			}
		}

		if(mcstep%p.epsSwap==0 && mcstep>100000)
		{
			if(rank==0)
			{
				if(drand48()<0.5)		// Only master node does 'coin flip'
				{
					structbox.whoSwap=0;
				}
				else
				{
					structbox.whoSwap=1;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&structbox.whoSwap,1,MPI_INT,0,MPI_COMM_WORLD);	// Master node tells the value of whoSwap
			MPI_Barrier(MPI_COMM_WORLD);

			if(structbox.whoSwap==0)				// Attempt swap 0-1, 2-3, ...
			{
				for(t=0;t<p.nT;t++)
				{
				for(l=0;l<p.neps-1;l+=2)
				{
					if(rank==p.procList[l][t] || rank==p.procList[l+1][t])		// Go into with 0-1, 2-3, ...
					{
						if(rank==p.procList[l][t])		// Proc 0,2,4,... send
						{
							epsSwapSend(&structbox,&p,p.procList[l+1][t]);
						}
						else							// Proc 1,3,5,... receive and send response
						{
							epsSwapRiceive(&structbox,p.procList[l][t],source,&p);
						}
						if(structbox.yestoSwap==1)		// Swap accepted
						{
							if(rank==p.procList[l][t])		// Proc 0,2,4,... send tempBox
							{
								SendTmpBox(&structbox,&tmpbox,&p,p.procList[l+1][t]);
							}
							else							// Proc 1,3,5,... receive StructBox
							{
								ReceiveTmpBox(&structbox,&tmpbox,&p,p.procList[l][t]);
							}

							structbox.curren = fullEnergy(&p,&structbox);
							printf("new en (eps swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.curren);
							structbox.currW = 0.0;
							structbox.cumuldev = 0.0;
							structbox.cumuldenom = 0.0;
							// configuration updated, recalculated constraint matrix
							calcBias(&p,&structbox);

// 							structbox.currW = IPDBias(&p,&structbox,1);
// 							calcDeviations(&p,&structbox,mcstep);
//
// 							structbox.currW = p.k_ipd*structbox.l_infty;
//
// 							for(ii=0;ii<p.nDChrom;ii++)
// 							{
// 								structbox.currW += surfEnergy(&p,&structbox,ii);
// 							}


							printf("new en (eps swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.currW);

							if(rank==p.procList[l][t])
							{
								printf("whoswap %d T: %g - %g ",structbox.whoSwap,source[ p.procList[l][t] ][0],source[ p.procList[l+1][t] ][0]);
								printf("MC %d: Swap in eps completed between %g e %g - ranks are %d %d\n",t,source[ p.procList[l][t] ][1],source[ p.procList[l+1][t] ][1],p.procList[l][t],p.procList[l+1][t]);
							}


						}

					}
				}
				}
			}

			else if(structbox.whoSwap ==1 && p.neps>2)			// Attempt swap 1-2, 3-4, ...
			{
				for(t=0;t<p.nT;t++)
				{
				for(l=1;l<p.neps-1;l+=2)
				{
					if(rank==p.procList[l][t] || rank==p.procList[l+1][t])		// Go into with 1-2, 3-4, ...
					{
						if(rank==p.procList[l][t])		// Proc 0,2,4,... send
						{
							epsSwapSend(&structbox,&p,p.procList[l+1][t]);
						}
						else							// Proc 1,3,5,... receive and send response
						{
							epsSwapRiceive(&structbox,p.procList[l][t],source,&p);
						}
						if(structbox.yestoSwap==1)	// Swap accepted
						{
							if(rank==p.procList[l][t])		// Proc 0,2,4,... send tempBox
							{
								SendTmpBox(&structbox,&tmpbox,&p,p.procList[l+1][t]);
							}
							else							// Proc 1,3,5,... receive StructBox
							{
								ReceiveTmpBox(&structbox,&tmpbox,&p,p.procList[l][t]);
							}

							structbox.curren = fullEnergy(&p,&structbox);
							printf("new en (eps swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.curren);
							structbox.currW = 0.0;
							structbox.cumuldev = 0.0;
							structbox.cumuldenom = 0.0;
							// configuration updated, recalculated constraint matrix
							calcBias(&p,&structbox);

// 							structbox.currW = IPDBias(&p,&structbox,1);
// 							calcDeviations(&p,&structbox,mcstep);
//
// 							structbox.currW = p.k_ipd*structbox.l_infty;
//
// 							for(ii=0;ii<p.nDChrom;ii++)
// 							{
// 								structbox.currW += surfEnergy(&p,&structbox,ii);
// 							}


							printf("new en (eps swap %d) on rank %d post swap %f\n",structbox.whoSwap,rank,structbox.currW);

							if(rank==p.procList[l][t])
							{
								printf("whoswap %dT: %g - %g",structbox.whoSwap,source[ p.procList[l][t] ][0],source[ p.procList[l+1][t] ][0]);
								printf("MC %d: Swap in eps completed between %g e %g - ranks are %d %d\n",t,source[ p.procList[l][t] ][1],source[ p.procList[l+1][t] ][1],p.procList[l][t],p.procList[l+1][t]);
							}
						}
					}
				}
				}
			}
		}

		if(mcstep%p.mcthermo==0 || mcstep == 0)
		{

			for(i=0;i<p.nDChrom;i++)
			{
				for(j=0;j<p.nDChrom;j++)
				{
					if (j>i)
					{
						trackerspam = monitorIPD(&p,&structbox,i,j);

					}

				}
				trackerspam = monitorCentroid(&p,&structbox,i);
				trackerspam = monitorHomologues(&p,&structbox,i);
			}
			structbox.cumuldev = 0.0;
			structbox.cumuldenom = 0.0;


			structbox.curren=fullEnergy(&p,&structbox);
// 			printf("thermo print\n");
			//Include code to check cosine similarity and L_inf max error also


			// configuration updated, recalculated constraint matrix
			calcBias(&p,&structbox);

// 			structbox.currW = IPDBias(&p,&structbox,1);
// 			calcDeviations(&p,&structbox,mcstep);
//
// 			structbox.currW = p.k_ipd*structbox.l_infty;
// // 			printf("rank %d step %d pre surf W %f\n",rank,mcstep,structbox.currW);
//
// 			for(i=0;i<p.nDChrom;i++)
// 			{
// 				structbox.currW += surfEnergy(&p,&structbox,i);
// 			}




			fprintf(structbox.fpmcthermo,"%d %f %f %f %f %f %d %d %d %f %f\n",mcstep,structbox.currW,structbox.curren,structbox.cosine_sim,structbox.l_infty,structbox.l1,structbox.maxdev,structbox.maxi,structbox.maxj,structbox.devdir,structbox.proxLinf);
			fflush(structbox.fpmcthermo);
			printf("rank %d step %d T %f W %g En %f cos_sim %f, L_inf %f, L1 %f, maxdev %d, ij %d %d,devdir %f, ex %f n constraint %d, proxLinf %f\n",
				   rank,mcstep,p.mctemp,structbox.currW,structbox.curren,structbox.cosine_sim,structbox.l_infty,structbox.l1,
		  structbox.maxdev,structbox.maxi,structbox.maxj,structbox.devdir,structbox.expectedMat[structbox.maxi][structbox.maxj],(int)((p.nChrom+1)*(p.nChrom+2)/2),structbox.proxLinf);
			// reset trackers for next step
			structbox.cumuldev = 0.0;
			structbox.cumuldenom = 0.0;
			structbox.maxavgIPDdist = 0.0;
			structbox.maxavgcentroiddist = 0.0;
			structbox.maxavghomologuedist = 0.0;


			// print out the expected mat and a spam 1 for the cosine similarity
			if(mcstep == 0)
			{
				fprintf(structbox.fpconstraint,"%f ",1.1);
				for(i=0;i<(p.nChrom+1)*(p.nChrom+1);i++)
				{
					indi = (int) i/(p.nChrom+1);
					indj = i % (p.nChrom+1);
					if(indj >= indi)
					{
						fprintf(structbox.fpconstraint,"%f ",structbox.expectedMat[indi][indj]);

					}
				}
				fprintf(structbox.fpconstraint,"%f\n",1.0);
				// print the toCalcMat also
				fprintf(structbox.fpconstraint,"%f ",1.2);
				for(i=0;i<(p.nChrom+1)*(p.nChrom+1);i++)
				{
					indi = (int) i/(p.nChrom+1);
					indj = i % (p.nChrom+1);
					if(indj >= indi)
					{
						fprintf(structbox.fpconstraint,"%d ",structbox.toCalc[indi][indj]);

					}
				}
				fprintf(structbox.fpconstraint,"%f\n",1.0);
			}

			// now print out the list of all constraints and the current cosine similarity in the last column
			fprintf(structbox.fpconstraint,"%d ",mcstep);
			for(i=0;i<(p.nChrom+1)*(p.nChrom+1);i++)
			{
				indi = (int) i/(p.nChrom+1);
				indj = i % (p.nChrom+1);
				if(indj >= indi)
				{
					fprintf(structbox.fpconstraint,"%f ",structbox.constraintMat[indi][indj]);

				}
			}
			fprintf(structbox.fpconstraint,"%f\n",structbox.cosine_sim);

			// print histogram
			if(mcstep > 200000)
			{
				sprintf(buffer2,"%s-L1Histo.dat",p.mcfile);
				structbox.fpHisto=fopen(buffer2,"w");
// 				printf("norm %f\n",structbox.norml1histo);
				for(i=0;i<p.nbinl1histo;i++)
				{
// 					deltacd = (1.0)/(1.0*p.nbinl1histo);
// 					printf("%d %d %f %f %f\n",mcstep,i,structbox.l1histo[i],structbox.norml1histo,structbox.l1histo[i]/structbox.norml1histo);
					fprintf(structbox.fpHisto,"%f %f\n",(i+0.5)/(1.0*p.nbinl1histo),structbox.l1histo[i]/structbox.norml1histo);
				}
				fclose(structbox.fpHisto);
			}




		}
// 			printf("printed thermo data\n");
		if(mcstep%p.mcdump==0)
		{
// 			printf("iter %d curr en %f -- curr temp %f -- others %d %d\n",mcstep,structbox.curren,p.mctemp,p.Tswap,p.epsSwap);
			sprintf(buffer2,"%s-dump_%d.dat",p.mcfile,mcstep);
// 			printf("dump file %s\n",buffer2);
			structbox.fpmcdump=fopen(buffer2,"w");
			fprintf(structbox.fpmcdump,"%d\n",p.nDChrom+p.nSurfPoints);
			fprintf(structbox.fpmcdump,"spam\n");

			for(ii=0;ii<p.nSurfPoints;ii++)
			{
				fprintf(structbox.fpmcdump,"%d %d %f %f %f %f %f\n",ii,20,0.8,0.5*avesep,structbox.surfPointx[ii],structbox.surfPointy[ii],structbox.surfPointz[ii]);
			}
			for(ii=0;ii<p.nDChrom;ii++)
			{
				i = structbox.chromList[ii];
				fprintf(structbox.fpmcdump,"%d %d %f %f %f %f %f\n",ii,structbox.chromType[ii],0.1,mspam*structbox.radius[i],structbox.X[i],structbox.Y[i],structbox.Z[i]);
			}
			fclose(structbox.fpmcdump);

			printf("writing separations\n");
			sprintf(buffer,"%s-separations.dat",p.mcfile);
			fpWrite=fopen(buffer,"w");
			writeSeparations(fpWrite,&p,&structbox);
			fclose(fpWrite);
			printf("wrote separations -- %s\n",buffer);
		}

		if(mcstep%(20*p.mcdump)==0)
		{
// 			checknTot = placePatchesVectorised(&p,&structbox);
// 			MPI_Barrier(MPI_COMM_WORLD);
// 			printf("checking %d %d %d\n",checknTot,p.nAto,p.nPatchTot);
//
// 			checknTot = placePatchesRandom(&p,&structbox);
// 			MPI_Barrier(MPI_COMM_WORLD);
// 			printf("checking again (add yourself) %d %d %d\n",checknTot,p.nAto,p.nPatchTot);
//
// 			checknTot = placeLaminPatches(&p,&structbox);
// 			MPI_Barrier(MPI_COMM_WORLD);
// 			printf("checking lamin (add yourself) %d %d %d\n",checknTot,p.nAto,p.nPatchTot);


			// Periodically write SimStart snapshot consumed by downstream mesh/simulation stages.
			printf("writing init conf\n");
			sprintf(buffer,"%s-%s.dat",p.mcfile,argv[3]);
			fpWrite=fopen(buffer,"w");
			dropSimStart(fpWrite,&p,&structbox);
			fclose(fpWrite);
			printf("wrote init conf -- %s\n",buffer);


		}

		if(mcstep> 0 && mcstep%p.anneal_freq==0)
		{

// 			p.mctemp += p.anneal_freq*deltemp;
			printf("%f %f %f -- %d %d\n",p.mctemp,deltemp,(1.0*p.anneal_freq)/(1.0*p.MCIter),p.anneal_freq,p.MCIter);
			p.mctemp *=deltemp;
		}
	}
	fclose(structbox.fpmcthermo);
	fclose(structbox.fpconstraint);

	// place patches
// 	int checknTot = placePatchesRandom(&p,&structbox);

	checknTot = placePatchesVectorised(&p,&structbox);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("checking %d %d %d\n",checknTot,p.nAto,p.nPatchTot);

	checknTot = placePatchesRandom(&p,&structbox);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("checking again (add yourself) %d %d %d\n",checknTot,p.nAto,p.nPatchTot);

	checknTot = placeLaminPatches(&p,&structbox);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("checking lamin (add yourself) %d %d %d\n",checknTot,p.nAto,p.nPatchTot);

	printf("writing init conf\n");
	sprintf(buffer,"%s-%s.dat",p.mcfile,argv[3]);
	fpWrite=fopen(buffer,"w");
	// Final SimStart handoff file.
	dropSimStart(fpWrite,&p,&structbox);
	fclose(fpWrite);
	printf("wrote init conf -- %s\n",buffer);

	int globCount=0;


	sprintf(buffer,"%s-visuglob%s",p.mcfile,argv[3]);
	printf("writing glob visu init conf -- %s\n",buffer);
	fpWrite=fopen(buffer,"w");
	fprintf(fpWrite,"%d\n",p.nAto+p.nSurfPoints);
	fprintf(fpWrite,"%lf\n",8.0*p.ella);

	for(ii=0;ii<p.nSurfPoints;ii++)
	{
		fprintf(fpWrite,"%d %d %f %f %f %f %f %f %f %f %d\n",globCount,200,0.8,0.5*avesep,structbox.surfPointx[ii],structbox.surfPointy[ii],structbox.surfPointz[ii],0.0,0.0,0.0,0);
		globCount++;
	}
	for(ii=0;ii<p.nDChrom;ii++)
	{
		// first write the info for this chrom
		//fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&chromindex,&type,&spamx,&spamy,&spamz,&orix,&oriy,&oriz,&rad,&mass,&k_spring,&act,&lamin);
		i=structbox.chromList[ii];

// 		printf("wrote %f, writing %f %d %d\n",p.ella,structbox.X[i],ii,i);
		fprintf(fpWrite,"%d %d %f %f %f %f %f %f %f %f %d\n",globCount,structbox.chromType[ii],0.1,mspam*structbox.radius[i],structbox.X[i],structbox.Y[i],structbox.Z[i],structbox.ox[i],structbox.oy[i],structbox.oz[i],0);
		globCount++;
// 		printf("printing chrom info no issue %d -- %d\n",ii,structbox.patchList[ii][0]);
		for(j=1;j<=structbox.patchList[ii][0];j++)
		{
			jj=structbox.patchList[ii][j];
			//printf("chrom %d raw %d patch %d raw %d\n",ii,i,j,jj);
			fprintf(fpWrite,"%d %d %f %f %f %f %f %f %f %f %d\n",globCount,structbox.patchPointer[ii][j],0.3,structbox.radius[jj],structbox.X[jj],structbox.Y[jj],structbox.Z[jj],structbox.actox[ii][j],structbox.actoy[ii][j],structbox.actoz[ii][j],structbox.lamin[ii][j]);


			globCount++;
		}
	}
	fclose(fpWrite);
	printf("wrote glob visu init conf -- %s\n",buffer);
	MPI_Barrier(MPI_COMM_WORLD);
	// Build IntMat entries consumed by prepSurf and simulation engines.
	sprintf(buffer,"%s-%s.dat",p.mcfile,argv[4]);
	fpWrite=fopen(buffer,"w");

	for(i=0;i<p.nAto-1;i++)
	{
		ii=structbox.whichChrom[i];
		ii2=structbox.chromType[ii];
		ii3=structbox.whichPatch[i];
		ii4=structbox.patchPointer[ii][ii3];

		// is there an entry eps_ij[ii3][ii4]?? if so, print it
		printf("iatom %d is on chrom %d (type %d) as patch %d -- points to %d\n",i,ii,ii2,ii3,ii4);
		if(ii3==-1)
			continue;


		// Scan candidate partner patches and emit only valid cross-type interactions.
		for(j=i+1;j<p.nAto;j++)
		{
			jj=structbox.whichChrom[j];
			jj2=structbox.chromType[jj];
			jj3=structbox.whichPatch[j];
			jj4=structbox.patchPointer[jj][jj3];
			printf("jatom %d is on chrom %d as patch %d (type %d)-- points to %d\n",j,jj,jj2,jj3,jj4);
			if(jj3==-1 || ii==jj || structbox.chromType[ii] == structbox.chromType[jj])
				continue;


			if(structbox.eps_ij[ii2][ii4]>0.01 && jj2==ii4 && jj4==ii2) // check if non-zero entry and if the chromtype for j is ii4, pointed to by i, the patch pointed to by j is of chromtype ii2
			{
				printf("%d %d  -- i %d i %d -- j %d j %d -- interaction %2.12f %2.12f %2.12f \n",i,j,ii2,jj4,jj2,ii4,structbox.eps_ij[ii2][ii4],structbox.eps_ij[ii2][jj2],structbox.eps_ij[jj4][jj2]);
				fprintf(fpWrite,"%d %d %d %d %d %d %f\n",i,j,ii,ii3,jj,jj3,structbox.eps_ij[ii2][jj2]);
			}
		}
	}


	fclose(fpWrite);
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(buffer2,"%s-constraintMat.dat",p.mcfile);
	printf("constraint mat file %s\n",buffer2);
	fpWrite=fopen(buffer2,"w");
	// Persist current and expected constraints for test/debug handoff.
	for(ii=0;ii<p.nChrom+1;ii++)
	{
		for(jj=0;jj<p.nChrom+1;jj++)
		{

			fprintf(fpWrite,"%d %d %f %f\n",ii,jj,structbox.constraintMat[ii][jj],structbox.expectedMat[ii][jj]);
		}
// 		fprintf(fpWrite,"\n");
	}
	fclose(fpWrite);



	MPI_Finalize();



}
