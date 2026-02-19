struct SYSTEM
{

	
	// Thermo variables
	double T_meas,press_meas,potEne,chromPot,surfPot,KE,totEne,IMEne;
	double lx,ly,lz,modRw;

	//Neighbour list related//
	double maxdispsq,probrij,probr0;
	int probi;

    // attributes
    double *radius,*mass,*invmass,*isBond,*eps_hertz,*eps_blob,*k_spring,*activity,*lamin,*oldactivity;
	double **eps_ij;
    int *type,*printflag;
    
	// Coordinates 
	
	double *X,*Y,*Z;
    double *trialX, *trialY, *trialZ;
    double *oldX,*oldY,*oldZ;
    
	// orientation
    double *qx,*qy,*qz,*qw;
    double *trialqx,*trialqy,*trialqz,*trialqw;
	double *oldqx,*oldqy,*oldqz,*oldqw; // chromosome orientation
	double *axisx,*axisy,*axisz;
	double *cosRotAngle;
	double *ox,*oy,*oz;
	double *actx,*acty,*actz;
	double *oldactx,*oldacty,*oldactz;
    
	double **patchposx,**patchposy,**patchposz,**Gpatchposx,**Gpatchposy,**Gpatchposz;
    double **trialGpatchposx,**trialGpatchposy,**trialGpatchposz; // patch positions chrom reference, patch active orientation global reference.
	
	// pointers organisers
	int *chromList;
    int **patchList,**patchPointer;
	int *whichChrom, *whichPatch;
	int *closesurf;
	
	// velocities
	double *vx,*vy,*vz;
	
	
    // forces
    double *fx,*fy,*fz;
	double *fRx,*fRy,*fRz;
	double *rgx,*rgy,*rgz;
	double *argx,*argy,*argz;


    // torques
    double *torquex,*torquey,*torquez;
    
	// auxiliary and dynamic
	double *en,*isIntermingling,*checkIntermingling;
// 	//FILES
 	FILE *fpthermo,*fprestart,*fpLAMMPSdump,*fpHisto,*fpHisto2,*fpTraj,*fpInit,*fpSurfPoints,*fpHalfPrint;
	FILE *fpHiC, *fpIPD, *fpRNG;
	
	
	
	// surface points
	
	double *surfPointx,*surfPointy,*surfPointz;
	double *surfFundx,*surfFundy,*surfFundz;
	double *oldsurfFundx,*oldsurfFundy,*oldsurfFundz;
	int *igridx,*igridy,*igridz;
	int *surfActive;
	int **surfcycpt;
	double *centSurfx,*centSurfy,*centSurfz;
	double *surffx,*surffy,*surffz;
	double *surfdispX,*surfdispY,*surfdispZ;
	double maxsurfdispsq;
	int **surfNeigh;
	int *surfneigpar;
	double *isMembIntermingling;
	double minsep,maxOffEll;
	int osmoflag;
	double sigtheta,currsig0x,currsig0z,currsig0y;
	double prevArea,currArea,stdArea;
	int activeAccum;


	// Analysis
	double **HiCMat;
	double *expectedMat,*constraintMat;


    char ens[3];

};

struct PARAINPUT
{

	int nChrom,nPatchTot,nAto,chromAct,seed;

	
	double Rw,rho,tem,gamma,rotgamma,pressure,thermostatFact,dt,origdt,taup,lamin_scale,epsij_scale,spring_scale,actscale;
	
	// file related
	int thermo,dump,startGeom,initRead,eqRun,totRun;
		char label[100];
	
		
	// repulsion related
		
	double eps_hertz,eps_blob;
	
	
	
	// surface related
	double ella,ellb,ellc,ellai,ellbi,ellci,ellaf,ellbf,ellcf,dea,deb,dec;
	int nSurfPoints,nSurfBasicPoints,nSurfCycPoints,nGrid,dewin,refineRounds;
	double surfPhi,surfNeighfact,SurfArea,delSurfArea,d_ar_thresh;



	// active stress fluctuations
	double tau_theta,sigactdel,sigdel;
	double sig0x,sig0y,sig0z;

	// geometry response
	double surftauinv,surfSpringSelf,surfDiffu,surfEta,surfrestore;

	// local geometry
	double surfSpringCross;

	// changing values

	double temi,temf,surfSpringCrossi,surfSpringCrossf,sigdeli,sigdelf,lamin_scalei,lamin_scalef,dtem,dssc,dsigdel,dls;
	double sigactdeli,sigactdelf,dsigactdel,actscalei,actscalef,dactscale;

	// analysis related

	double hic_cut,hic_decay;

	const gsl_rng_type * gsl_T;
	gsl_rng * gsl_r;


};

struct APPOGGIO
{
	//Random number related//            
            

 
	unsigned long *mt;
	int mti;
	int **NEIGLISTold;
	float *buff_float;
	double *Xold,*Yold,*Zold,*dispXold,*dispYold,*dispZold;
	double lxold, lyold, lzold;

	double vold,vnew,lnvnew;
	double arg, resfactX, resfactY, resfactZ;
	double penew, pediff;
	
	double deltaBeta;
};


void inizialized(struct PARAINPUT *p,struct SYSTEM *structbox,struct SYSTEM *tmpbox)
{
	int i;


	p->gsl_T = gsl_rng_default;
	p->gsl_r = gsl_rng_alloc (p->gsl_T);

    
    // attributes
    structbox->type=calloc(p->nAto,sizeof(int));
    structbox->printflag=calloc(p->nAto,sizeof(int));
    
	structbox->mass=calloc(p->nAto,sizeof(double));
	structbox->invmass=calloc(p->nAto,sizeof(double));
    structbox->radius=calloc(p->nAto,sizeof(double));
    structbox->eps_hertz=calloc(p->nAto,sizeof(double));
    structbox->eps_blob=calloc(p->nAto,sizeof(double));
    structbox->k_spring=calloc(p->nAto,sizeof(double));
    structbox->isBond=calloc(p->nAto,sizeof(double));
	structbox->activity=calloc(p->nAto,sizeof(double));
	structbox->oldactivity=calloc(p->nAto,sizeof(double));
	structbox->lamin=calloc(p->nAto,sizeof(double));
	
	structbox->eps_ij=calloc(p->nAto,sizeof(double *));
	for(i=0;i<p->nAto;i++)
	{
		structbox->eps_ij[i] = calloc(p->nAto,sizeof(double));
	}
    
    //positions
	structbox->X=calloc(p->nAto,sizeof(double));
	structbox->Y=calloc(p->nAto,sizeof(double));
	structbox->Z=calloc(p->nAto,sizeof(double));
    structbox->trialX=calloc(p->nAto,sizeof(double));
	structbox->trialY=calloc(p->nAto,sizeof(double));
	structbox->trialZ=calloc(p->nAto,sizeof(double));
	structbox->oldX=calloc(p->nAto,sizeof(double));
	structbox->oldY=calloc(p->nAto,sizeof(double));
	structbox->oldZ=calloc(p->nAto,sizeof(double));
	
// 	structbox->lx = structbox->ly = structbox->lz = 8.0*p->Rw;
	//velocity
	structbox->vx=calloc(p->nAto,sizeof(double));
 	structbox->vy=calloc(p->nAto,sizeof(double));
 	structbox->vz=calloc(p->nAto,sizeof(double));
    
    // forces
    
    structbox->fx=calloc(p->nAto,sizeof(double));
 	structbox->fy=calloc(p->nAto,sizeof(double));
 	structbox->fz=calloc(p->nAto,sizeof(double));
    
    structbox->fRx=calloc(p->nAto,sizeof(double));
 	structbox->fRy=calloc(p->nAto,sizeof(double));
 	structbox->fRz=calloc(p->nAto,sizeof(double));
    
    structbox->rgx=calloc(p->nAto,sizeof(double));
 	structbox->rgy=calloc(p->nAto,sizeof(double));
 	structbox->rgz=calloc(p->nAto,sizeof(double));

	structbox->argx=calloc(p->nAto,sizeof(double));
 	structbox->argy=calloc(p->nAto,sizeof(double));
 	structbox->argz=calloc(p->nAto,sizeof(double));


    
	
	// torques 
	
	structbox->torquex=calloc(p->nAto,sizeof(double));
	structbox->torquey=calloc(p->nAto,sizeof(double));
	structbox->torquez=calloc(p->nAto,sizeof(double));

    //orientation
	
	structbox->axisx=calloc(p->nAto,sizeof(double));
	structbox->axisy=calloc(p->nAto,sizeof(double));
	structbox->axisz=calloc(p->nAto,sizeof(double));
	
	structbox->cosRotAngle=calloc(p->nAto,sizeof(double));

    structbox->qx=calloc(p->nAto,sizeof(double));
	structbox->qy=calloc(p->nAto,sizeof(double));
	structbox->qz=calloc(p->nAto,sizeof(double));
    structbox->qw=calloc(p->nAto,sizeof(double));
    structbox->trialqx=calloc(p->nAto,sizeof(double));
	structbox->trialqy=calloc(p->nAto,sizeof(double));
	structbox->trialqz=calloc(p->nAto,sizeof(double));
    structbox->trialqw=calloc(p->nAto,sizeof(double));
    structbox->oldqx=calloc(p->nAto,sizeof(double));
	structbox->oldqy=calloc(p->nAto,sizeof(double));
	structbox->oldqz=calloc(p->nAto,sizeof(double));
    structbox->oldqw=calloc(p->nAto,sizeof(double));
	
    structbox->ox=calloc(p->nAto,sizeof(double));
	structbox->oy=calloc(p->nAto,sizeof(double));
	structbox->oz=calloc(p->nAto,sizeof(double));
    
	structbox->actx=calloc(p->nAto,sizeof(double));
 	structbox->acty=calloc(p->nAto,sizeof(double));
 	structbox->actz=calloc(p->nAto,sizeof(double));

	structbox->oldactx=calloc(p->nAto,sizeof(double));
 	structbox->oldacty=calloc(p->nAto,sizeof(double));
 	structbox->oldactz=calloc(p->nAto,sizeof(double));


    // chrom and chrom patchs
	
	structbox->whichChrom=calloc(p->nAto,sizeof(int));	
	structbox->whichPatch=calloc(p->nAto,sizeof(int));
	structbox->closesurf=calloc(p->nAto,sizeof(int));

	
    structbox->chromList=calloc(p->nChrom,sizeof(int));
	structbox->patchposx=calloc(p->nChrom,sizeof(double *));
	structbox->patchposy=calloc(p->nChrom,sizeof(double *));
	structbox->patchposz=calloc(p->nChrom,sizeof(double *));
    structbox->Gpatchposx=calloc(p->nChrom,sizeof(double *));
	structbox->Gpatchposy=calloc(p->nChrom,sizeof(double *));
	structbox->Gpatchposz=calloc(p->nChrom,sizeof(double *));
    structbox->trialGpatchposx=calloc(p->nChrom,sizeof(double *));
	structbox->trialGpatchposy=calloc(p->nChrom,sizeof(double *));
	structbox->trialGpatchposz=calloc(p->nChrom,sizeof(double *));
	structbox->patchList=calloc(p->nChrom,sizeof(int *));
	structbox->patchPointer=calloc(p->nChrom,sizeof(int *));
	for(i=0;i<p->nChrom;i++)
	{
		structbox->patchposx[i]=calloc(p->nPatchTot+1,sizeof(double));
		structbox->patchposy[i]=calloc(p->nPatchTot+1,sizeof(double));
		structbox->patchposz[i]=calloc(p->nPatchTot+1,sizeof(double));
        structbox->Gpatchposx[i]=calloc(p->nPatchTot+1,sizeof(double));
		structbox->Gpatchposy[i]=calloc(p->nPatchTot+1,sizeof(double));
		structbox->Gpatchposz[i]=calloc(p->nPatchTot+1,sizeof(double));
        structbox->trialGpatchposx[i]=calloc(p->nPatchTot+1,sizeof(double));
		structbox->trialGpatchposy[i]=calloc(p->nPatchTot+1,sizeof(double));
		structbox->trialGpatchposz[i]=calloc(p->nPatchTot+1,sizeof(double));
		structbox->patchList[i]=calloc(p->nPatchTot+1,sizeof(int));
		structbox->patchPointer[i]=calloc(p->nPatchTot+1,sizeof(int));

	}
	
	
	// auxiliary
	
	structbox->en=calloc(p->nAto,sizeof(double));
	structbox->isIntermingling=calloc(p->nAto,sizeof(double));
	structbox->checkIntermingling=calloc(p->nAto,sizeof(double));
	structbox->isMembIntermingling=calloc(p->nSurfPoints,sizeof(double));


	// Analysis
	int nhap = (int)(0.5*p->nChrom);
	structbox->HiCMat=calloc(nhap,sizeof(double *));

	for(i=0;i<nhap;i++)
	{
		structbox->HiCMat[i] = calloc(nhap,sizeof(double));
	}
// 	structbox->IPD = calloc(nhap+1,sizeof(double *));
// 	for(i=0;i<nhap+1;i++)
// 	{
// 		structbox->IPD[i] = calloc(nhap+1,sizeof(double));
// 	}

	int nconstraint = (nhap + 1)*nhap - (int) (nhap*(nhap-1)/2) + nhap+1 - nhap;
	structbox->expectedMat = calloc(nconstraint,sizeof(double));
	structbox->constraintMat = calloc(nconstraint,sizeof(double));


	
	
	// surface points
	structbox->surfPointx=calloc(p->nSurfPoints,sizeof(double));
	structbox->surfPointy=calloc(p->nSurfPoints,sizeof(double));
	structbox->surfPointz=calloc(p->nSurfPoints,sizeof(double));

	structbox->surfFundx=calloc(p->nSurfPoints,sizeof(double));
	structbox->surfFundy=calloc(p->nSurfPoints,sizeof(double));
	structbox->surfFundz=calloc(p->nSurfPoints,sizeof(double));

	structbox->oldsurfFundx=calloc(p->nSurfPoints,sizeof(double));
	structbox->oldsurfFundy=calloc(p->nSurfPoints,sizeof(double));
	structbox->oldsurfFundz=calloc(p->nSurfPoints,sizeof(double));

	structbox->surffx=calloc(p->nSurfPoints,sizeof(double));
	structbox->surffy=calloc(p->nSurfPoints,sizeof(double));
	structbox->surffz=calloc(p->nSurfPoints,sizeof(double));

	structbox->igridx=calloc(p->nSurfPoints,sizeof(int));
	structbox->igridy=calloc(p->nSurfPoints,sizeof(int));
	structbox->igridz=calloc(p->nSurfPoints,sizeof(int));

	structbox->surfActive=calloc(p->nSurfPoints,sizeof(int));
	for(i=0;i<p->nSurfBasicPoints;i++)
		structbox->surfActive[i] = 1;

	for(i=p->nSurfBasicPoints;i<p->nSurfPoints;i++)
		structbox->surfActive[i]=0;

	structbox->surfcycpt=calloc(p->nSurfCycPoints,sizeof(int *));

	for(i=0;i<p->nSurfCycPoints;i++)
		structbox->surfcycpt[i]=calloc(3,sizeof(int));

	structbox->centSurfx=calloc(p->nSurfPoints,sizeof(double));
	structbox->centSurfy=calloc(p->nSurfPoints,sizeof(double));
	structbox->centSurfz=calloc(p->nSurfPoints,sizeof(double));


	structbox->surfdispX=calloc(p->nSurfPoints,sizeof(double));
	structbox->surfdispY=calloc(p->nSurfPoints,sizeof(double));
	structbox->surfdispZ=calloc(p->nSurfPoints,sizeof(double));

	structbox->surfNeigh=calloc(p->nSurfPoints+1,sizeof(int *));
	structbox->surfneigpar=calloc(p->nSurfPoints+1,sizeof(int));
	for(i=0;i<=p->nSurfPoints;i++)
	{
		structbox->surfNeigh[i]=calloc(p->nSurfPoints+1,sizeof(int));
	}
	
// 	a->mti=N+1;
// 	a->mt=calloc(N,sizeof(unsigned long));
// 	a->buff_float=calloc(p->nChrom,sizeof(float *));

  
}



