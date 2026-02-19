#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gsl/gsl_rng.h>

#define MAX(a,b)	((a)>(b))?(a):(b)
#define MIN(a,b)	((a)<(b))?(a):(b)
#define cube(x)	((x)*(x)*(x))
#define sqr(x)		((x)*(x))

#define SKIN 		1.2
#define minGrow	0.95
#define N 			624
#define minDist2	0.8
#define minPatchDist2	0.5
#define LINESIZE 1024
#define alat 1.675


#include "functions/Structures.h"
#include "functions/RigidBody.h"
#include "functions/FileInput.h"
#include "functions/AuxFileInput.h"
#include "functions/Surface_constrained.h"
#include "functions/Random.h"

static FILE *open_or_die(const char *path, const char *mode)
{
	FILE *fp = fopen(path, mode);
	if(fp == NULL)
	{
		printf("Error: could not open %s\n", path);
		exit(1);
	}
	return fp;
}

static void run_cmd_or_die(const char *cmd)
{
	int rc = system(cmd);
	if(rc != 0)
	{
		printf("Error: command failed (%d): %s\n", rc, cmd);
		exit(1);
	}
}

int main (int argc, char *argv[])
{
	if(argc!=5 && argc!=6)
	{
		printf("Usage: %s <Para File> <Topo> <Int_Mat> <Output> [ConstraintsVector]\n", argv[0]);
		exit(1);
	}
	const char *expected_mat_path = (argc == 6) ? argv[5] : "ConstraintsVector.dat";
	
	struct SYSTEM structbox;
	struct SYSTEM tmpbox;
	struct PARAINPUT p;
	
	
	FILE *fpLog;
	char filepart[200],filename[500],sffilename[500],cmdstr[1000];
	
	int spampairs=0;
	int mdstart=0;




	
	fpLog=open_or_die(argv[1],"r");
	readPara(fpLog,&p,&structbox);
	printf("read para ref rounds - %d\n",p.refineRounds);
	fclose(fpLog);
	inizialized(&p,&structbox,&tmpbox);
	
	printf("Done input file %s!\n",argv[1]);
	snprintf(filepart,sizeof(filepart),"%s/%s",argv[4],p.label);

	// Load the expected constraints vector used to initialize constraint bookkeeping.
	structbox.fpIPD=open_or_die(expected_mat_path,"r");
	readExpectedMat(&p,&structbox);
	fclose(structbox.fpIPD);
	printf("read expected constraint vector\n");

	structbox.currsig0x=p.sig0x;
	structbox.currsig0y=p.sig0y;
	structbox.currsig0z=p.sig0z;

	printf("going to genrand %d\n",p.seed);
	gsl_rng_set(p.gsl_r,p.seed);
	printf("first rand num %ld -- %f\n",gsl_rng_get(p.gsl_r),gsl_rng_uniform(p.gsl_r));

	fpLog=open_or_die(argv[2],"r");
	printf("initialisation\n");
	readTopo(fpLog,&p,&structbox);
	printf("read topo\n");
	fclose(fpLog);
	p.Rw = pow(3.0*p.nChrom/(4.0*M_PI*p.rho),0.3333);
	printf("done topo file %s!\n",argv[2]);

	printf("reading interactions --%s\n",argv[3]);
	fpLog=open_or_die(argv[3],"r");
	spampairs=readIntMatrix(fpLog,&p,&structbox);
	printf("read interaction matrix! %d pairwise interactions specified\n",spampairs);
	fclose(fpLog);
	snprintf(filename,sizeof(filename),"%s-surfConn.dat",filepart);
	fpLog=open_or_die(filename,"w");
	printf("entering surf gen -- %d\n",p.refineRounds);

	// initRead=0: generate base surface, refine it, and save final refined mesh.
	if(p.initRead==0)
	{
		printf("generating new surf -- %d rounds\n",p.refineRounds);
		genSurf(fpLog,&p,&structbox);
		fclose(fpLog);

		MCRelaxSurf(&p,&structbox,50000);

		for(int rd =1;rd<=p.refineRounds;rd++)
		{

			// Dump current mesh, run refinement, then load the refined result.
			snprintf(sffilename,sizeof(sffilename),"%s-SurfFile_%d.dat",filepart,rd);
			fpLog=open_or_die(sffilename,"w");
			dropSurf(fpLog,&p,&structbox);
			fclose(fpLog);
			snprintf(cmdstr,sizeof(cmdstr),"python3 fixSurf_randSel.py %s %d %s-fixedSurfFile_ea%g_eb%g_ec%g_%d.dat %f %f %f",
				sffilename,p.nSurfBasicPoints,filepart,p.ellai,p.ellbi,p.ellci,rd,p.ellai,p.ellbi,p.ellci);
			printf("%s\n",cmdstr);
			run_cmd_or_die(cmdstr);

			snprintf(sffilename,sizeof(sffilename),"%s-fixedSurfFile_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,rd);
			fpLog=open_or_die(sffilename,"r");
			readSurf(fpLog,&p,&structbox);
			fclose(fpLog);

			// Smooth after each refinement round before the next drop/read cycle.
			MCRelaxSurf(&p,&structbox,100000);
		}
		snprintf(sffilename,sizeof(sffilename),"%s-finalfixedSurfFile_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds);
		printf("%s",sffilename);
		fpLog=open_or_die(sffilename,"w");
		dropSurf(fpLog,&p,&structbox);
		fclose(fpLog);
	}
	// initRead=1: relax an already-created surface and rewrite final mesh.
	else if(p.initRead==1)
	{
		fclose(fpLog);
		MCRelaxSurf(&p,&structbox,100000);

		snprintf(sffilename,sizeof(sffilename),"%s-finalfixedSurfFile_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds);
		printf("%s",sffilename);
		fpLog=open_or_die(sffilename,"w");
		dropSurf(fpLog,&p,&structbox);
		fclose(fpLog);
	}
	// initRead=2: load an existing final surface and continue downstream.
	else if(p.initRead==2)
	{
		fclose(fpLog);
		snprintf(sffilename,sizeof(sffilename),"%s-finalfixedSurfFile_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds);
		fpLog=open_or_die(sffilename,"r");
		readSurf(fpLog,&p,&structbox);
		fclose(fpLog);
	}
	else
	{
		fclose(fpLog);
		printf("Error: unsupported initRead value %d\n", p.initRead);
		exit(1);
	}

	printf("done initial prep ---v%d\n",p.initRead);

	// Generate cycle points from the final refined surface.
	printf("generating cycle points %d\n",p.nSurfCycPoints);
	snprintf(sffilename,sizeof(sffilename),"%s-SurfCyc_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds);
	snprintf(cmdstr,sizeof(cmdstr),"python3 getNHole_perturb.py %s-finalfixedSurfFile_ea%g_eb%g_ec%g_%d.dat %d %d %d %f %f %f %f %s",
		filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds,p.nSurfBasicPoints,p.nSurfCycPoints,p.nGrid,p.ellai,p.ellbi,p.ellci,p.surfNeighfact,sffilename);
	printf("%s",cmdstr);
	run_cmd_or_die(cmdstr);

	snprintf(filename,sizeof(filename),"%s-surfPoints_%d.dat",filepart,mdstart);
	structbox.fpSurfPoints=open_or_die(filename,"w");
	writeSurfPoints(&p,&structbox,mdstart);
	fclose(structbox.fpSurfPoints);


	printf("generated prepped surf and surface cycles\n");
	gsl_rng_free(p.gsl_r);
	return 0;
}
