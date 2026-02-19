struct SYSTEM
{
    double *X,*Y,*Z;
    double *ox,*oy,*oz;
    double *fx,*fy,*fz;
    double **patchposx,**patchposy,**patchposz;

    double **actox,**actoy,**actoz;

    int **patchList,**patchPointer,**patchDone,**lamin;
    int *chromList, *chromType;
    int *whichChrom,*whichPatch;
    double *radius,*mass;

    double **ipd,**ipddev,**ipdhomo;
    double **eps_ij;
    double **constraintMat,**expectedMat;
    int **toCalc;

    double *chromDist,*chromDistDev,*chrompatchmass,*chromcoremass,*chromlaminmass,*chrompatchrad,*chromcorerad,*chromlaminrad;
    double *l1histo;

    double *surfPointx,*surfPointy,*surfPointz;

    double maxavgIPDdist,maxavghomologuedist,maxavgcentroiddist,l_infty,l1,cosine_sim,norml1histo,proxLinf,devdir;
    double olden,curren,cumuldev,cumuldenom,countConstraints,currW;
    int maxdev;
    int maxi,maxj,chi1,chi2,chj1,chj2;



    /// replica and bias related
    double *Xpre,*Ypre,*Zpre;
    double deltaBeta;

    int whoSwap,yestoSwap;

    FILE *fpmcthermo,*fpmcdump,*fpconstraint,*fpHisto,*fpRestart;
};

struct PARAINPUT
{
    int nChrom,nDChrom,nPatchTot,nAto,anneal_freq;
    int gradMove;
    double ella,ellb,ellc,epsij_rescale,lamin_scale,massscalefact,k_spring_level,f_act_level;
    double k_ipd,k_homologue,k_centroid,eps_hertz;
    double spamipddev;

    int MCIter,mcthermo,mcdump,seed;
    double mcdisp,mctemp,dtemp,constraint_freq,IPD_scale;
    int fipd,fhomo,fcent;

    int nbinl1histo;

    int nSurfPoints;
    double surfPhi;
    char IPDfile[500];
    char IPDHomofile[500];
    char ChromPropsfile[500];
    char IADfile[500];
    char mcfile[500];

    int nT,neps,Tswap,epsSwap,restartPt;
    int **procList;

};

void readPara(FILE *fp, struct PARAINPUT *p, struct SYSTEM *box)
{
	char string[100];
	double spam;
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->nChrom=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->nPatchTot=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ella=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ellb=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ellc=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->epsij_rescale=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->MCIter=atoi(string);
//     if(fscanf(fp, "%s",string)>0); 				// thermostatFact
// 	if(fscanf(fp, "%s",string)>0);
// 	p->mctemp=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->mcdisp=atof(string);
//     if(fscanf(fp, "%s",string)>0); 				// thermostatFact
// 	if(fscanf(fp, "%s",string)>0);
//     sprintf(p->mcfile,"%s",string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->mcthermo=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->mcdump=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->nSurfPoints=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfPhi=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->lamin_scale=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->massscalefact=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->k_ipd=atof(string);
    p->k_homologue=0.0;
    p->k_centroid=0.0;
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->spamipddev=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->k_spring_level=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->f_act_level=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->dtemp=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->anneal_freq=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->constraint_freq=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->IPD_scale=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->Tswap=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->epsSwap=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->restartPt=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->gradMove=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->fipd=atoi(string);
    p->fhomo=0;
    p->fcent=0;
}


void initAll(struct PARAINPUT *p,struct SYSTEM *structbox,struct SYSTEM *tmpbox)
{
    int i,j;

    structbox->maxavgcentroiddist = 0.0;
    structbox->maxavgIPDdist = 0.0;
    structbox->maxavghomologuedist = 0.0;
    structbox->cumuldev = 0.0;
    structbox->cumuldenom = 0.0;
    structbox->mass=calloc(p->nAto,sizeof(double));
    structbox->radius=calloc(p->nAto,sizeof(double));

    structbox->X=calloc(p->nAto,sizeof(double));
	structbox->Y=calloc(p->nAto,sizeof(double));
	structbox->Z=calloc(p->nAto,sizeof(double));

    structbox->fx=calloc(p->nAto,sizeof(double));
	structbox->fy=calloc(p->nAto,sizeof(double));
	structbox->fz=calloc(p->nAto,sizeof(double));

    structbox->ox=calloc(p->nAto,sizeof(double));
	structbox->oy=calloc(p->nAto,sizeof(double));
	structbox->oz=calloc(p->nAto,sizeof(double));

	structbox->actox=calloc(p->nDChrom,sizeof(double *));
 	structbox->actoy=calloc(p->nDChrom,sizeof(double *));
 	structbox->actoz=calloc(p->nDChrom,sizeof(double *));

    structbox->patchposx=calloc(p->nDChrom,sizeof(double *));
    structbox->patchposy=calloc(p->nDChrom,sizeof(double *));
    structbox->patchposz=calloc(p->nDChrom,sizeof(double *));

    structbox->patchList=calloc(p->nDChrom,sizeof(int *));
    structbox->patchPointer=calloc(p->nDChrom,sizeof(int *));
    structbox->patchDone=calloc(p->nDChrom,sizeof(int *));
    structbox->lamin=calloc(p->nDChrom,sizeof(int *));

    for(i=0;i<p->nDChrom;i++) // plus two for lamin patch + buffer
    {
        structbox->actox[i]=calloc(p->nPatchTot+2,sizeof(double));
        structbox->actoy[i]=calloc(p->nPatchTot+2,sizeof(double));
        structbox->actoz[i]=calloc(p->nPatchTot+2,sizeof(double));

        structbox->patchposx[i]=calloc(p->nPatchTot+2,sizeof(double));
        structbox->patchposy[i]=calloc(p->nPatchTot+2,sizeof(double));
        structbox->patchposz[i]=calloc(p->nPatchTot+2,sizeof(double));

        structbox->patchList[i]=calloc(p->nPatchTot+2,sizeof(int));
        structbox->patchDone[i]=calloc(p->nPatchTot+2,sizeof(int));
        structbox->patchPointer[i]=calloc(p->nPatchTot+2,sizeof(int));

        structbox->lamin[i]=calloc(p->nPatchTot+2,sizeof(int));

    }

    structbox->whichChrom=calloc(p->nAto,sizeof(int));
	structbox->whichPatch=calloc(p->nAto,sizeof(int));

    structbox->chromList=calloc(p->nDChrom,sizeof(int));
    structbox->chromType=calloc(p->nDChrom,sizeof(int));
    structbox->chromDist=calloc(p->nDChrom,sizeof(double));
    structbox->chromDistDev=calloc(p->nDChrom,sizeof(double));

    structbox->chrompatchmass=calloc(p->nDChrom,sizeof(double));
    structbox->chromcoremass=calloc(p->nDChrom,sizeof(double));
    structbox->chromlaminmass=calloc(p->nDChrom,sizeof(double));

    structbox->chrompatchrad=calloc(p->nDChrom,sizeof(double));
    structbox->chromcorerad=calloc(p->nDChrom,sizeof(double));
    structbox->chromlaminrad=calloc(p->nDChrom,sizeof(double));

    structbox->eps_ij=calloc(p->nChrom+1,sizeof(double *));
    structbox->constraintMat=calloc(p->nChrom+1,sizeof(double *));
    structbox->expectedMat=calloc(p->nChrom+1,sizeof(double *));
    structbox->toCalc=calloc(p->nChrom+1,sizeof(int *));
    for(i=0;i<p->nChrom+1;i++)
    {

        structbox->eps_ij[i]=calloc(p->nChrom+1,sizeof(double *));
        structbox->constraintMat[i]=calloc(p->nChrom+1,sizeof(double *));
        structbox->expectedMat[i]=calloc(p->nChrom+1,sizeof(double *));
        structbox->toCalc[i]=calloc(p->nChrom+1,sizeof(int *));
    }


    structbox->ipd=calloc(p->nDChrom,sizeof(double *));
    structbox->ipddev=calloc(p->nDChrom,sizeof(double *));
    structbox->ipdhomo=calloc(p->nDChrom,sizeof(double *));

    for(i=0;i<p->nDChrom;i++)
    {
        structbox->ipd[i]=calloc(p->nDChrom,sizeof(double *));
        structbox->ipddev[i]=calloc(p->nDChrom,sizeof(double *));
        structbox->ipdhomo[i]=calloc(2,sizeof(double *));

    }
    structbox->countConstraints=0;
    structbox->surfPointx=calloc(p->nSurfPoints,sizeof(double));
    structbox->surfPointy=calloc(p->nSurfPoints,sizeof(double));
    structbox->surfPointz=calloc(p->nSurfPoints,sizeof(double));

    p->nbinl1histo = 20;
    structbox->l1histo = calloc(p->nbinl1histo,sizeof(double));




    structbox->Xpre=calloc(p->nAto,sizeof(double));
	structbox->Ypre=calloc(p->nAto,sizeof(double));
	structbox->Zpre=calloc(p->nAto,sizeof(double));


    tmpbox->X=calloc(p->nAto,sizeof(double));
    tmpbox->Y=calloc(p->nAto,sizeof(double));
    tmpbox->Z=calloc(p->nAto,sizeof(double));

    tmpbox->surfPointx=calloc(p->nSurfPoints,sizeof(double));
    tmpbox->surfPointy=calloc(p->nSurfPoints,sizeof(double));
    tmpbox->surfPointz=calloc(p->nSurfPoints,sizeof(double));
}
