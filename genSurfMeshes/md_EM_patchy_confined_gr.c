#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_coupling.h>
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
// #include "functions/Patchy_hertz.h"
#include "functions/Patchy_IntHertzSurfNorm.h"
// #include "functions/Integrator_EMActive_wTorque.h"
// #include "functions/Integrator_OABP_wTorque.h"
#include "functions/Integrator_OABP_wTorque_gr.h"
#include "functions/Analysis.h"
//#include "functions/NeighList.h"
// #include "functions/create_patchy_spheres.h"


int main (int argc, char *argv[])
{
	if(argc!=5)
	{
		printf("Error: Launch %s <Para File> <Topo> <Int_Mat> <Output> %d \n", argv[0],argc);
		exit(1);
	}
	
	struct SYSTEM structbox;
	struct SYSTEM tmpbox;
	struct APPOGGIO app;
	struct PARAINPUT p;
	
	
	FILE *fpLog,*fpThermo,*fpRestart;
	char filepart[200],filename[500],thermofilename[500],sffilename[500],cmdstr[1000];
	
	int t=0,spampairs=0;
	int mdstart;
	double rgx,rgtheta,rgz,sigtheta;

	
	fpLog=fopen(argv[1],"r");					// Input file (paraUS.dat)
	if(fpLog==NULL) 
	{
		printf("Error: %s not found!\n",argv[1]);
		exit(1);
	}
	readPara(fpLog,&p,&structbox);
	fclose(fpLog);
// 	inizialized(&p,&structbox,&tmpbox,&app);
	inizialized(&p,&structbox,&tmpbox);

// 	exit(0);
	
	printf("Done input file %s!\n",argv[1]);
	sprintf(filepart,"%s/%s",argv[4],p.label);

	// read the expected constraints vector
	structbox.fpIPD=fopen("ConstraintsVector.dat","r");
	readExpectedMat(&p,&structbox);
	printf("read expected constraint vector\n");

	structbox.currsig0x=p.sig0x;
	structbox.currsig0y=p.sig0y;
	structbox.currsig0z=p.sig0z;
	if(p.initRead==0)					// initread 0 = Create atoms
	{
//     structbox.fpInit=fopen(argv[4],"r");
		printf("going to genrand %d\n",p.seed);
// 		init_genrand(p.seed,&app);
		gsl_rng_set(p.gsl_r,p.seed);
		printf("first rand num %ld -- %f\n",gsl_rng_get(p.gsl_r),gsl_rng_uniform(p.gsl_r));

		fpLog=fopen(argv[2],"r");
		printf("initialisation -- %s\n",argv[2]);
		readTopo(fpLog,&p,&structbox);
		printf("read topo\n");
		fclose(fpLog);
		p.Rw = pow(3.0*p.nChrom/(4.0*M_PI*p.rho),0.3333);
		printf("done topo file %s!\n",argv[2]);
// 		init_vel(&p,&structbox);
		mdstart=0.0;
		
		printf("reading interactions --%s\n",argv[3]);
		fpLog=fopen(argv[3],"r");
		spampairs=readIntMatrix(fpLog,&p,&structbox);
		printf("read interaction matrix! %d pairwise interactions specified\n",spampairs);
		fclose(fpLog);
		


		// read doubly refined
		sprintf(sffilename,"%s-finalfixedSurfFile_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds);
		printf("opening %s\n",sffilename);
		fpLog=fopen(sffilename,"r");
		readSurf(fpLog,&p,&structbox);
		fclose(fpLog);

		// read cycles
		sprintf(sffilename,"%s-SurfCyc_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds);
		printf("opening %s\n",sffilename);
		fpLog=fopen(sffilename,"r");
		readCyc(fpLog,&p,&structbox);
		fclose(fpLog);

		findGrid(&p,&structbox,0); // also places surfpoints outside surfFundx
		vnlist(&p,&structbox);


		sprintf(filename,"%s-surfPoints_%d.dat",filepart,mdstart);
		structbox.fpSurfPoints=fopen(filename,"w");
		writeSurfPoints(&p,&structbox,mdstart);
		fclose(structbox.fpSurfPoints);

		structbox.prevArea = structbox.stdArea;
		structbox.currArea = structbox.stdArea;




		printf("generated surface and identified grid!\n");

		sprintf(thermofilename,"%s-thermo.dat",filepart);
		fpThermo = fopen(thermofilename,"w");

		sprintf(filename,"%s-Constraints.dat",filepart);
		structbox.fpIPD = fopen(filename,"w");

		structbox.sigtheta = M_PI+atan(structbox.currsig0z/structbox.currsig0x) + 2.0*gsl_rng_uniform(p.gsl_r) - 1.0;
		printf("initial angle %f -- %f %f\n",structbox.sigtheta,structbox.currsig0z,structbox.currsig0x);
		
	}
	
	if(p.initRead ==1)
	{
		printf("we ready to deal with all this fancy stuff\n");


		printf("going to genrand RESTART  %d\n",p.seed);
// 		init_genrand(p.seed,&app);
		gsl_rng_set(p.gsl_r,p.seed);
// 		printf("first rand num %ld -- %f\n",gsl_rng_get(p.gsl_r),gsl_rng_uniform(p.gsl_r));

		// NOTE: read RNG state here!!!
		sprintf(filename,"%s-RNGstart.dat",filepart);

		printf("rng restart %s\n",filename);
		structbox.fpRNG=fopen(filename,"rb");
		gsl_rng_fread(structbox.fpRNG,p.gsl_r);
		fclose(structbox.fpRNG);


		/// get config from restart file

		sprintf(filename,"%s-Restart.dat",filepart);
		printf("RESTART initialisation %s\n",filename);

		fpRestart=fopen(filename,"r");
		mdstart=readSimRestart(fpRestart,&p,&structbox);
		fclose(fpRestart);

		printf("RESTART read config -- %d\n",mdstart);
		p.Rw = pow(3.0*p.nChrom/(4.0*M_PI*p.rho),0.3333);
		printf("RESTART done file %s!\n",filename);
// 		init_vel(&p,&structbox);
// 		mdstart=0.0;

		// read interactions
		printf("RESTART reading interactions --%s\n",argv[3]);
		fpLog=fopen(argv[3],"r");
		spampairs=readIntMatrix(fpLog,&p,&structbox);
		printf("RESTART read interaction matrix! %d pairwise interactions specified\n",spampairs);
		fclose(fpLog);

		// read surface cycles -- need to get triangle indices
		printf("reading cycle first -- fund and point will be reassigned by read surfPoints\n");
		sprintf(sffilename,"%s-SurfCyc_ea%g_eb%g_ec%g_%d.dat",filepart,p.ellai,p.ellbi,p.ellci,p.refineRounds);
		fpLog=fopen(sffilename,"r");
		readCyc(fpLog,&p,&structbox);
		fclose(fpLog);

		printf("read cycle -- reading points\n");

		/// read surface points -- updates the fund positions
		sprintf(filename,"%s-surfPoints_%d.dat",filepart,mdstart);
		printf("reading surf from %s\n",filename);
		structbox.fpSurfPoints=fopen(filename,"r");
		readSurfPoints(&p,&structbox,mdstart);
		fclose(structbox.fpSurfPoints);

		printf("finding grid\n");
		findGrid(&p,&structbox,1); // make no modification to surfPoints -- just assign grids

		vnlist(&p,&structbox); // have to explicitly call neigh list because no call to neighlist

		printf("reading hic\n");
		/// read HiC
		sprintf(filename,"%s-HiCMat.dat",filepart);
		structbox.fpHiC=fopen(filename,"r");
		readCurrHiC(&p,&structbox);
		fclose(structbox.fpHiC);



		sprintf(thermofilename,"%s-thermo.dat",filepart);
		fpThermo = fopen(thermofilename,"a");

		sprintf(filename,"%s-Constraints.dat",filepart);
		structbox.fpIPD = fopen(filename,"a");

// 		structbox.sigtheta = M_PI+atan(structbox.currsig0z/structbox.currsig0x) + 2.0*gsl_rng_uniform(p.gsl_r) - 1.0;
// 		printf("initial angle %f -- %f %f\n",structbox.sigtheta,structbox.currsig0z,structbox.currsig0x);
	}
	
// 	lcount=vnlist(&p,&structbox);							// Build neighbor list
// 	printf("built first neighbourlist\n");

	// NOTE: calling energyForce resets forces and then calculates based on position
	printf("doing first energy calculation\n");
	structbox.potEne=energyForce(&structbox,&p);			// Compute potential energy
    kinEnergy(&structbox,&p);
    structbox.totEne=structbox.potEne+structbox.KE;

	// create neighbourlist (careful with radii etc -- leave it for later)




	
	// later interactions are patch-wall and a hertzian/viscoelastic chromosome-chromosome interaction which can have COM and torque components
	// update positions based on force and orientation based on torque if any
	int annealcount=0;
// 	sigtheta = M_PI+atan(currsig0z/currsig0x) + 2.0*drand48() - 1.0;

		structbox.osmoflag=0;
	printf("starting dynamics -- %d %d steps\n",p.eqRun,p.totRun);
	for(t=mdstart;t<p.totRun;t++)
    {
// 		if(t>=p.eqRun)
// 			structbox.osmoflag=1;
		if(t%p.thermo==0 || t==0)
		{
			printf("Step %d -- IMe %f PE %f chPE %f sfPE %f KE %f geom %f (%f) %f (%f) %f (%f) (angle) %f  \
			-- str %f %f restore %f %f act %f %f -- %f %f\n",t,structbox.IMEne/p.epsij_scale,structbox.potEne,structbox.chromPot, \
			structbox.surfPot,structbox.KE*2.0/(3.0*p.nAto),p.ella,p.ellai,p.ellb,p.ellbi,p.ellc,p.ellci, \
			structbox.sigtheta - atan(p.sig0z/p.sig0x),structbox.currsig0x,structbox.currsig0z,p.tau_theta*p.dt*(p.sig0x - structbox.currsig0x), \
			p.tau_theta*p.dt*(p.sig0z - structbox.currsig0z),cos(structbox.sigtheta - atan(p.sig0z/p.sig0x)), \
			sin(structbox.sigtheta - atan(p.sig0z/p.sig0x)),sqrt(2.0*p.dt*p.sig0x*p.sigdel),sqrt(2.0*p.dt*p.sig0x*p.sigactdel));

			fprintf(fpThermo,"%d %f %f %f %f %f %f %f %f %f\n",t,structbox.IMEne/p.epsij_scale, \
			structbox.potEne,structbox.KE,p.ella,p.ellb,p.ellc,structbox.currsig0x,structbox.currsig0z,structbox.sigtheta);

			fflush(fpThermo);
		}
		if(t%p.dump==0 || t==0 || (t<p.dump && t%100 ==0))
		{
// 			printf("printing output %d\n",t);
// 			sprintf(filename,"%s-Chromdump_%d.dat",filepart,t);
// 			printf("output file is %s\n",filename);
// 			structbox.fpLAMMPSdump=fopen(filename,"w");
// 			writeChromdump(&p,&structbox,&app,t);
// 			fclose(structbox.fpLAMMPSdump);

// 			printf("printing output %d\n",t);
			sprintf(filename,"%s-Fulldump_%d.dat",filepart,t);
			printf("output file is %s\n",filename);
			structbox.fpLAMMPSdump=fopen(filename,"w");
			writeFulldump(&p,&structbox,&app,t);
			fclose(structbox.fpLAMMPSdump);

// 			sprintf(filename,"%s-Transpdump_%d.dat",filepart,t);
// 			printf("output file is %s\n",filename);
// 			structbox.fpLAMMPSdump=fopen(filename,"w");
// 			writeTranspdump(&p,&structbox,&app,t);
// 			fclose(structbox.fpLAMMPSdump);

			sprintf(filename,"%s-Anadump_%d.dat",filepart,t);
			printf("output file is %s\n",filename);
			structbox.fpLAMMPSdump=fopen(filename,"w");
			writeAnadump(&p,&structbox,&app,t);
			fclose(structbox.fpLAMMPSdump);
			
			
			sprintf(filename,"%s-surfPoints_%d.dat",filepart,t);
			structbox.fpSurfPoints=fopen(filename,"w");
			writeSurfPoints(&p,&structbox,t);
			fclose(structbox.fpSurfPoints);

// 			sprintf(filename,"%s-HalfOpen_%d.dat",filepart,t);
// 			structbox.fpHalfPrint=fopen(filename,"w");
// // 			writeHalfOpen(&p,&structbox,t);
// 			fclose(structbox.fpHalfPrint);

			// print analysis files (HiC -- IPD)
			getProxyHiC(&p,&structbox);
			sprintf(filename,"%s-HiCMat.dat",filepart);
			structbox.fpHiC=fopen(filename,"w");
			writeCurrHiC(&p,&structbox);
			fclose(structbox.fpHiC);


			getIPD(&p,&structbox);
			writeConstraints(&p,&structbox,t);
// 			fflush(structbox.fpIPD);

		}
		if (t%(10*p.dump)==0)
		{
			// NOTE: Write restart config, restart surf, save rng state

			sprintf(filename,"%s-Restart.dat",filepart);
			fpRestart=fopen(filename,"w");
			printf("RESTART initialisation\n");
			writeSimRestart(fpRestart,&p,&structbox,t);
			printf("RESTART read config\n");
			fclose(fpRestart);

			sprintf(filename,"%s-surfPoints_%d.dat",filepart,t);
			structbox.fpSurfPoints=fopen(filename,"w");
			writeSurfPoints(&p,&structbox,t);
			fclose(structbox.fpSurfPoints);

			// print analysis files (HiC -- IPD)
			getProxyHiC(&p,&structbox);
			sprintf(filename,"%s-HiCMat.dat",filepart);
			structbox.fpHiC=fopen(filename,"w");
			writeCurrHiC(&p,&structbox);
			fclose(structbox.fpHiC);

			sprintf(filename,"%s-RNGstart.dat",filepart);
			structbox.fpRNG=fopen(filename,"wb");
			gsl_rng_fwrite(structbox.fpRNG,p.gsl_r);
			fclose(structbox.fpRNG);
		}
		
        // update -- adds random force and/or torque (NOTE: later add active force) and then updates positions
        integrateLangevinEM(&structbox,&p);
		integrateLangevinEMSurf(&structbox,&p);
		updateGeom(&p,&structbox,structbox.currsig0x,structbox.currsig0z);

		// calc interaction -- resets forces and torques
        structbox.potEne=energyForce(&structbox,&p);
		 // calc kin energy
        kinEnergy(&structbox,&p);
        structbox.totEne=structbox.potEne+structbox.KE;

		// active stress fluctuation direction
		rgtheta = gauss_var_AT(0,1,&p);
		structbox.sigtheta = structbox.sigtheta + sqrt(2*p.dt/p.tau_theta)*rgtheta;




		// update stressFlucFreq
		rgx = gauss_var_AT(0,1,&p);
		rgz = gauss_var_AT(0,1,&p);


// 		currsig0x = currsig0x + p.dt*(p.sig0x - currsig0x + p.sigactdel*cos(sigtheta)) + sqrt(2.0*p.dt*p.sigdel)*rgx;
// 		currsig0z = currsig0z + p.dt*(p.sig0z - currsig0z + p.sigactdel*sin(sigtheta)) + sqrt(2.0*p.dt*p.sigdel)*rgz ;

// 		currsig0x = currsig0x + p.tau_theta*p.dt*(p.sig0x - currsig0x) + sqrt(2.0*p.dt*p.sig0x*p.sigactdel)*cos(sigtheta - atan(p.sig0z/p.sig0x)) + sqrt(2.0*p.dt*p.sig0x*p.sigdel)*rgx;
// 		currsig0z = currsig0z + p.tau_theta*p.dt*(p.sig0z - currsig0z) + sqrt(2.0*p.dt*p.sig0z*p.sigactdel)*sin(sigtheta - atan(p.sig0z/p.sig0x)) + sqrt(2.0*p.dt*p.sig0z*p.sigdel)*rgz ;

		structbox.currsig0x = structbox.currsig0x + p.surfrestore*p.dt*(p.sig0x - structbox.currsig0x) + sqrt(2.0*p.dt*p.sig0x*p.sigactdel)*cos(structbox.sigtheta) + sqrt(2.0*p.dt*p.sig0x*p.sigdel)*rgx;

		structbox.currsig0z = structbox.currsig0z + p.surfrestore*p.dt*(p.sig0z - structbox.currsig0z) + sqrt(2.0*p.dt*p.sig0z*p.sigactdel)*sin(structbox.sigtheta) + sqrt(2.0*p.dt*p.sig0z*p.sigdel)*rgz ;


        if(t>=p.eqRun && t < 0.75*p.totRun && t%p.dewin==0)
		{
			annealcount+=1;
			p.tem = p.temi + p.dtem*annealcount;
			p.surfSpringCross = p.surfSpringCrossi + p.dssc*annealcount;
			p.sigdel = p.sigdeli + p.dsigdel*annealcount;
			p.lamin_scale = p.lamin_scalei + p.dls*annealcount;


// 			p.sigactdel = p.sigactdeli + annealcount*(p.sigactdel - 0.001)*p.dewin/(1.0*(0.75*p.totRun - p.eqRun));
			p.sigactdel = p.sigactdeli + p.dsigactdel*annealcount;


			p.sig0x = (p.ellai + p.dea*annealcount)*(p.surfSpringSelf/p.surfEta)/(p.surftauinv/p.surfEta);
// 			p.sig0y = (p.ellbi + p.deb*annealcount)*(p.surfSpringSelf/p.surfEta)/(p.surftauinv/p.surfEta);
			p.sig0z = (p.ellci + p.dec*annealcount)*(p.surfSpringSelf/p.surfEta)/(p.surftauinv/p.surfEta);

			if(p.ellai != p.ellaf && p.ellci != p.ellcf)
			{
				p.deb = (p.ellbf - p.ellb)*p.dewin/(1.0*(0.75*p.totRun - t));
				printf("current eb %f -- will change to %f -- %f %f\n",p.ellb,p.ellbi + p.deb*annealcount,p.deb,0.75*p.totRun-t);
				p.ellb = p.ellbi + p.deb; //*annealcount;
				updateGeom(&p,&structbox,structbox.currsig0x,structbox.currsig0z);
// 			p.sig0x = (p.ellai + p.dea*annealcount)*p.surfSpringSelf;
// // 			p.sig0y = (p.ellbi + p.deb*annealcount)*p.surfSpringSelf;
// 			p.sig0z = (p.ellci + p.dec*annealcount)*p.surfSpringSelf;

				printf("annealing geometry e %f sig %f -- e %f sig %f -- e %f sig %f -- tem %f spc %f str fluc %f\n",p.ella,p.sig0x,p.ellb,p.sig0y,p.ellc,p.sig0z,p.tem,p.surfSpringCross,p.sigdel);
			}
        // output configuration and thermo data
		}
        
    }   
	fclose(fpThermo);
	fclose(structbox.fpIPD);
	gsl_rng_free(p.gsl_r);
	
}
