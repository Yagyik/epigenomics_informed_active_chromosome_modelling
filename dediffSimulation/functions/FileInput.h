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
	p->nAto = p->nChrom + p->nPatchTot;
	if(fscanf(fp, "%s",string)>0); 				// tem
	if(fscanf(fp, "%s",string)>0);
	p->temi=atof(string);
	p->tem=p->temi;
	if(fscanf(fp, "%s",string)>0); 				// tem
	if(fscanf(fp, "%s",string)>0);
	p->temf=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// tem
	if(fscanf(fp, "%s",string)>0);
	p->dt=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// gamma (6*pi*eta)
	if(fscanf(fp, "%s",string)>0);
	p->gamma=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// rotation gamma (6*pi*eta)
	if(fscanf(fp, "%s",string)>0);
	p->rotgamma=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// activity same for full chromosome if ==1, diff for each patch if ==0
	if(fscanf(fp, "%s",string)>0);
	p->chromAct=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// activity scaler for patches
	if(fscanf(fp, "%s",string)>0);
	p->actscalei=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// activity scaler for patches
	if(fscanf(fp, "%s",string)>0);
	p->actscalef=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// rho
	if(fscanf(fp, "%s",string)>0);
	p->rho=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// pressure
	if(fscanf(fp, "%s",string)>0);
	p->pressure=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// initRead
	if(fscanf(fp, "%s",string)>0);
	p->initRead=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// startGeom
	if(fscanf(fp, "%s",string)>0);
	p->startGeom=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// seed
	if(fscanf(fp, "%s",string)>0);
	p->seed=atoi(string);
	printf("seed is %d\n",p->seed);
	if(fscanf(fp, "%s",string)>0);				// Label
	if(fscanf(fp, "%s",p->label)>0);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->thermostatFact=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->eqRun=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->totRun=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->thermo=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->dump=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ellai=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ellbi=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ellci=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ellaf=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->ellcf=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->dewin=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->eps_blob=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->spring_scale=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->epsij_scale=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->eps_hertz=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->nSurfBasicPoints=atoi(string);
	printf("surf points -- %d\n",p->nSurfBasicPoints);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfNeighfact=atof(string);
	printf("read surf %f\n",p->surfNeighfact);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->nSurfCycPoints=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->refineRounds=atoi(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->d_ar_thresh=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfPhi=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfSpringSelf=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfSpringCrossi=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfSpringCrossf=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfDiffu=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->lamin_scalei=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->lamin_scalef=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->taup=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->sigdeli=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->sigdelf=atof(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->nGrid=atoi(string);
    if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->tau_theta=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->sigactdeli=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->sigactdelf=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surftauinv=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->surfrestore=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->hic_cut=atof(string);
	if(fscanf(fp, "%s",string)>0); 				// thermostatFact
	if(fscanf(fp, "%s",string)>0);
	p->hic_decay=atof(string);


	p->nSurfPoints = p->nSurfBasicPoints + p->nSurfCycPoints;
    p->origdt=p->dt;

	p->surfEta = 1.0/(p->surfDiffu); // ultimately to be rescaled with the number of surf fpSurfPoints

	p->tem=p->temi;
	p->surfSpringCross=p->surfSpringCrossi;
	p->sigdel=p->sigdeli;
	p->lamin_scale=p->lamin_scalei;
	p->sigactdel=p->sigactdeli;
	p->actscale=p->actscalei;

	spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	box->stdArea = 4.0*M_PI*pow(spam,1.0/1.6);


// 	p->ellcf = p->ellai*p->ellbi*p->ellci/(p->ellaf*p->ellbf);

	// volume preserving
// 	p->ellbf = p->ellai*p->ellbi*p->ellci/(p->ellaf*p->ellcf);

	// non-volume preserving
	p->ellbf = p->ellbi - 0.5*(p->ellai - p->ellaf + p->ellci - p->ellcf);


	p->ella=p->ellai;
	p->ellb=p->ellbi;
	p->ellc=p->ellci;


	p->sig0x = p->ellai*(p->surfSpringSelf/p->surfEta)/(p->surftauinv/p->surfEta);
	p->sig0y = p->ellbi*(p->surfSpringSelf/p->surfEta)/(p->surftauinv/p->surfEta);
	p->sig0z = p->ellci*(p->surfSpringSelf/p->surfEta)/(p->surftauinv/p->surfEta);



	p->dea=(p->ellaf - p->ellai)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));
	p->deb=(p->ellbf - p->ellbi)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));
	p->dec=(p->ellcf - p->ellci)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));






	p->dtem = (p->temf - p->temi)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));
	p->dssc = (p->surfSpringCrossf - p->surfSpringCrossi)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));
	p->dsigdel = (p->sigdelf - p->sigdeli)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));
	p->dls = (p->lamin_scalef - p->lamin_scalei)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));
	p->dsigactdel = (p->sigactdelf - p->sigactdeli)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));
	p->dactscale = (p->actscalef - p->actscalei)*p->dewin/(1.0*(0.5*p->totRun - p->eqRun));

	printf("CHECKING !!!!!! %f %f %f %f %f %f\n",p->dtem,p->dssc,p->dsigdel,p->dls,p->dsigactdel,p->dactscale);


}

void readTopo(FILE *fp, struct PARAINPUT *p, struct SYSTEM *box)
{
	int i,j,ii,jj,spamindex;
	int nLines;
	char string[100];
	double rotpatch[3];

	double spam;
	int chromindex, type;
	double spamx,spamy,spamz,orix,oriy,oriz,angle,pot_factor,rad,mass,k_spring,act,lamin;
	
	double norm2,rotpatchnorm,quatnorm;
	printf("reading topo\n");
	fscanf(fp,"%s",string);
	box->lx = atof(string); // extent of square box, replace later
    box->ly = box->lx;
    box->lz = box->lx;
	printf("chrom, patch, atoms, box -- %d %d %d %f\n",p->nChrom,p->nPatchTot,p->nAto,box->lx);
	for(i=0;i<p->nAto;i++) 
	{
		// each line has type -chrom/patch, npatches (next n are patches for this particle and are placed on body coords wrt COM of chrom), orientation vector, radius, mass, spring constant
		fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&chromindex,&type,&spamx,&spamy,&spamz,&orix,&oriy,&oriz,&rad,&mass,&k_spring,&act,&lamin);
// 		printf("checking reads %d %d %f %f %f %f %f %f %f %f %f %f %f\n",chromindex,type,spamx,spamy,spamz,orix,oriy,oriz,rad,mass,k_spring,act,lamin);
		// if particle type chromosome, assign particle type, pos, orientation,potential factor, radius,mass
		if(type > 0) // a chromosome, therefore store appropriate info
		{
			box->type[i]=type;
			box->X[i] = spamx;
			box->Y[i] = spamy;
			box->Z[i] = spamz;
			box->oldX[i] = spamx;
			box->oldY[i] = spamy;
			box->oldZ[i] = spamz;
			
			box->ox[i] = orix;
			box->oy[i] = oriy;
			box->oz[i] = oriz;
			
// 			norm2 = sqrt(box->ox[i]*box->ox[i] + box->oy[i]*box->oy[i] + box->oz[i]*box->oz[i]);
// 			printf("norm of this vector is %f\n",norm2);
//			printf("chrom %d pos %f %f %f\n",i,box->X[i],box->Y[i],box->Z[i]);
// 			box->ox[i] /= norm2;
// 			box->oy[i] /= norm2;
// 			box->oz[i] /= norm2;
			
            // set attributes
// 			box->eps_hertz[i] = pot_factor;
			box->radius[i] = rad;			
			box->mass[i] = mass;
			box->invmass[i] = 1.0/mass;
			box->activity[i] = p->chromAct*act*p->actscale; //NOTE::: remove this hard-coded obscenity
			box->lamin[i] = p->lamin_scale*lamin; // NOTE: lamin_scale modifier is multiplied to maintain consistency with energy scales
			box->chromList[chromindex-1] = i; //the core of the chromindexth chromosome
			box->whichChrom[i] = chromindex - 1;
			box->whichPatch[i] = -1;
			j=0; // reset j for patch positions
			
		}
		// if particle type patch, assign type 0 global position based on body position of last chrom and it's ori vector -- then put potential factor, radius, mass and spring constant
		if(type == 0) // not a chromosome, a patch
		{
			box->type[i] = type;
// 			box->X[i] = spamx; // body coordinate positions -- will have to be modified
// 			box->Y[i] = spamy;
// 			box->Z[i] = spamz;
// 			box->oldX[i] = spamx;
// 			box->oldY[i] = spamy;
// 			box->oldZ[i] = spamz;
			j++; //next patch
			box->whichChrom[i] = chromindex - 1;
			box->whichPatch[i] = j;
			box->patchList[chromindex-1][0]++; // one more patch
			box->patchList[chromindex-1][j] = i; // particle i is connected to the jth patch on this chromosome
			
			box->actx[i] = orix;
			box->acty[i] = oriy;
			box->actz[i] = oriz;

			norm2 = sqrt(box->actx[i]*box->actx[i] + box->acty[i]*box->acty[i] + box->actz[i]*box->actz[i]);
// 			printf("norm of this vector is %f\n",norm2);

			box->actx[i] /= norm2;
			box->acty[i] /= norm2;
			box->actz[i] /= norm2;

			box->oldactx[i] = box->actx[i];
			box->oldacty[i] = box->acty[i];
			box->oldactz[i] = box->actz[i];

			// body-frame patch coordinates -- radius agnostic and immutable
			box->patchposx[chromindex-1][0]++; 
			box->patchposy[chromindex-1][0]++;
			box->patchposz[chromindex-1][0]++;
			box->patchposx[chromindex-1][j] = spamx;
			box->patchposy[chromindex-1][j] = spamy;
			box->patchposz[chromindex-1][j] = spamz;
			
			box->radius[i] = rad;
// 			box->radius[i] = 5.0;
			box->mass[i] = mass;
			box->invmass[i] = 1.0/mass;
// 			box->eps_blob[i] = pot_factor;
			box->k_spring[i] = p->spring_scale*k_spring;
			if(lamin > 0.0)
			box->activity[i] = 0.0;
			else
			box->activity[i] = act*p->actscale;
			box->lamin[i] = p->lamin_scale*lamin;
			box->isBond[i] = 1.0; // NOTE: very bad idea but placeheld to later be replaced with a 2D array
			
		}
		
		
		// also assign tether position as body coord position on surface of 
	}
	printf("finding axes for init ori\n");
	findAxesInitOri(p,box);
	// now rotate to find global coordinates of each patch
	printf("found axes, rotating patches\n");
	for(ii=0;ii<p->nChrom;ii++)
	{
		// set quarternion for ith chromosome based on axis and angle
		i = box->chromList[ii];
//		printf("chrom %d par %d\n",ii,i);
		box->qx[i] = box->axisx[i]*sin(0.5*acos(box->cosRotAngle[i]));
		box->qy[i] = box->axisy[i]*sin(0.5*acos(box->cosRotAngle[i]));
		box->qz[i] = box->axisz[i]*sin(0.5*acos(box->cosRotAngle[i]));
		box->qw[i] = cos(0.5*acos(box->cosRotAngle[i]));
		
		box->oldqx[i] = box->qx[i];
		box->oldqy[i] = box->qy[i];
		box->oldqz[i] = box->qz[i];
		box->oldqw[i] = box->qw[i];
		
		quatnorm = sqrt(box->oldqx[i]*box->oldqx[i] + box->oldqy[i]*box->oldqy[i] + box->oldqz[i]*box->oldqz[i] + box->oldqw[i]*box->oldqw[i]);
		
		box->oldqx[i] /= quatnorm;
		box->oldqy[i] /= quatnorm;
		box->oldqz[i] /= quatnorm;
		box->oldqw[i] /= quatnorm;
// 		printf("set quat from rotations chrom %d par %d\n",ii,i);
		for(j=1;j<=box->patchList[ii][0];j++)
		{
			rotate(box->patchposx[ii][j],box->patchposy[ii][j],box->patchposz[ii][j],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],rotpatch);
			rotpatchnorm=sqrt(rotpatch[0]*rotpatch[0] + rotpatch[1]*rotpatch[1] + rotpatch[2]*rotpatch[2]);
			// now assign the global positions of the patch positions -- COM +
			box->Gpatchposx[ii][j] = box->X[i] + box->radius[i]*rotpatch[0]/rotpatchnorm; // position of patches
			box->Gpatchposy[ii][j] = box->Y[i] + box->radius[i]*rotpatch[1]/rotpatchnorm;
			box->Gpatchposz[ii][j] = box->Z[i] + box->radius[i]*rotpatch[2]/rotpatchnorm;
//
			// for now the real positions are the patch positions + some small displacement
			spamindex = box->patchList[ii][j];
			box->X[spamindex] = box->Gpatchposx[ii][j] + 0.5*box->radius[i]*(2.0*drand48()-1.0); // position of connected bead
			box->Y[spamindex] = box->Gpatchposy[ii][j] + 0.5*box->radius[i]*(2.0*drand48()-1.0);
			box->Z[spamindex] = box->Gpatchposz[ii][j] + 0.5*box->radius[i]*(2.0*drand48()-1.0);
			
			box->oldX[spamindex] = box->X[spamindex];
			box->oldY[spamindex] = box->Y[spamindex];
			box->oldZ[spamindex] = box->Z[spamindex];
			
// 			printf("ch %d patch %d (real %d %d) Init patch pos %f %f %f -- bead pos %f %f %f\n",ii,j,i,spamindex,box->Gpatchposx[ii][j],box->Gpatchposy[ii][j],box->Gpatchposz[ii][j],box->oldX[spamindex],box->oldY[spamindex],box->oldZ[spamindex]);
			
		}
	}

// 	exit(0);

}

int readIntMatrix(FILE *fp, struct PARAINPUT *p, struct SYSTEM *box)
{
	// subroutine to read out the NxM matrix of interactions
	// first read chromi1 patchj1 and chromi2 patchj2 then read the epsilon value
	// interaction matrix is 0 for all others
	int j,j1,j2,cp1,cp2,pp1,pp2;
	int scp1,spp1,scp2,spp2;
	double spameps;
	int npairs=0;
	rewind(fp);
	int ii,jj;
	printf("whut whut\n");
	for(ii=0;ii<p->nChrom;ii++)
	{
//		printf("chrom %d has %d patches \n",ii,box->patchList[ii][0]);
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
			j = box->patchList[ii][jj];
// 			printf("chrom %d patch %d -- glob %d lamin %f\n",ii,jj,j,box->lamin[j]);
			if(box->lamin[j]==0) // not lamin
			box->patchPointer[ii][jj] = p->nAto+2;
			else
			box->patchPointer[ii][jj] = p->nAto+3;
		}
	}

	while(!feof(fp))
	{
		fscanf(fp,"%d %d %d %d %d %d %lf\n",&j1,&j2,&scp1,&spp1,&scp2,&spp2,&spameps);
		box->eps_ij[j1][j2] = p->epsij_scale*spameps;
		box->eps_ij[j2][j1] = p->epsij_scale*spameps;
// 		printf("read int mat for %d %d as %f\n",j1,j2,box->eps_ij[j1][j2]);
		pp1 = box->whichPatch[j1];
		pp2 = box->whichPatch[j2];
		cp1 = box->whichChrom[j1];
		cp2 = box->whichChrom[j2];
// 		printf("assigning to c1 %d p1 %d -- c2 %d p2 %d -- compare %d %d %d %d\n",cp1,pp1,cp2,pp2,scp1,spp1,scp2,spp2);

		box->patchPointer[cp1][pp1] = box->type[box->chromList[cp2]]-1;
		box->patchPointer[cp2][pp2] = box->type[box->chromList[cp1]]-1;

		npairs++;
	}
// 	exit(0);
	return npairs;
}

// NOTE: Subroutine below deprecated
int readIntMatrixOld(FILE *fp, struct PARAINPUT *p, struct SYSTEM *box)
{
	// subroutine to read out the NxM matrix of interactions
	// first read chromi1 patchj1 and chromi2 patchj2 then read the epsilon value
	// interaction matrix is 0 for all others
	int spamci1,spamci2,spampj1,spampj2,j1,j2;
	double spameps;
	int npairs=0;
	rewind(fp);
	while(!feof(fp))
	{
		fscanf(fp,"%d %d %d %d %lf\n",&spamci1,&spampj1,&spamci2,&spampj2,&spameps);
		if(spampj1 <=0 || spampj2 <=0)
		{
			printf("patch pos must be 1 or higher check interactions file\n");
			exit(0);
		}
		if(spamci1 <=0 || spamci2 <=0)
		{
			printf("chrom indices must be 1 or, will be mod in code to index correctly, higher check interactions file\n");
			exit(0);
		}
		j1 = box->patchList[spamci1-1][spampj1];
		j2 = box->patchList[spamci2-1][spampj2];
		box->eps_ij[j1][j2] = spameps;
		box->eps_ij[j2][j1] = spameps;
// 		printf("read int mat for %d %d as %f\n",j1,j2,box->eps_ij[j1][j2]);
		npairs++;
	}
// 	exit(0);
	return npairs;
}
//



int readSimRestart(FILE *fp,struct PARAINPUT *p, struct SYSTEM *box)
{
	int i,ii,j,jj,chromindex,spamindex,type,mdstep;
	double spamx,spamy,spamz,orix,oriy,oriz,oriw,pot_factor,rad,mass,k_spring,act,lamin,spampx,spampy,spampz,spamgx,spamgy,spamgz;

	int nLines;
	double ex,ey,ez,sigx,sigz,sigt,prevArea;
	double spam,angle,norm2,rotpatchnorm,quatnorm,spamdiff;
	char string[100];
	double rotpatch[3];
	// first read the MD step
	fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",&mdstep,&spam,&ex,&ey,&ez,&sigx,&sigz,&sigt,&prevArea);
	box->lx = spam;
	box->ly = box->lx;
    box->lz = box->lx;
	p->ella = ex;
	p->ellb = ey;
	p->ellc = ez;
	spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	box->currArea = 4.0*M_PI*pow(spam,1.0/1.6);
	box->prevArea = prevArea;
	box->currsig0x = sigx;
	box->currsig0z = sigz;
	box->sigtheta = sigt;

	for(i=0;i<p->nAto;i++)
	{
// 		fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", \
// 		&chromindex,&type,&spamx,&spamy,&spamz,&orix,&oriy,&oriz,&spampx,&spampy,&spampz,&spamgx,&spamgy,&spamgz, \
// 		&rad,&mass,&k_spring,&act,&lamin);

		fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", \
		&chromindex,&type,&spamx,&spamy,&spamz,&orix,&oriy,&oriz,&oriw,&spampx,&spampy,&spampz,&spamgx,&spamgy,&spamgz, \
		&rad,&mass,&k_spring,&act,&lamin);

		if(type > 0) // a chromosome, therefore store appropriate info
		{
			box->type[i]=type;
			box->X[i] = spamx;
			box->Y[i] = spamy;
			box->Z[i] = spamz;
			box->oldX[i] = spamx;
			box->oldY[i] = spamy;
			box->oldZ[i] = spamz;

// 			box->ox[i] = orix;
// 			box->oy[i] = oriy;
// 			box->oz[i] = oriz;

			box->axisx[i] = orix;
			box->axisy[i] = oriy;
			box->axisz[i] = oriz;
			box->cosRotAngle[i] = oriw;

            // set attributes
// 			box->eps_hertz[i] = pot_factor;
			box->radius[i] = rad;
			box->mass[i] = mass;
			box->invmass[i] = 1.0/mass;
			box->activity[i] = p->chromAct*act*p->actscale; //NOTE::: remove this hard-coded obscenity
			box->lamin[i] = p->lamin_scale*lamin; // NOTE: lamin_scale modifier is multiplied to maintain consistency with energy scales
			box->chromList[chromindex-1] = i; //the core of the chromindexth chromosome
			box->whichChrom[i] = chromindex - 1;
			box->whichPatch[i] = - 1;
			j=0; // reset j for patch positions

		}
		// if particle type patch, assign type 0 global position based on body position of last chrom and it's ori vector -- then put potential factor, radius, mass and spring constant
		if(type == 0) // not a chromosome, a patch
		{
			box->type[i] = type;
			box->X[i] = spamx; // body coordinate positions -- will have to be modified
			box->Y[i] = spamy;
			box->Z[i] = spamz;
			box->oldX[i] = spamx;
			box->oldY[i] = spamy;
			box->oldZ[i] = spamz;
			j++; //next patch
			box->whichChrom[i] = chromindex - 1;
			box->whichPatch[i] = j;
			box->patchList[chromindex-1][0]++; // one more patch
			box->patchList[chromindex-1][j] = i; // particle i is connected to the jth patch on this chromosome

			box->actx[i] = orix;
			box->acty[i] = oriy;
			box->actz[i] = oriz;

// 			norm2 = sqrt(box->actx[i]*box->actx[i] + box->acty[i]*box->acty[i] + box->actz[i]*box->actz[i]);
// // 			printf("norm of this vector is %f\n",norm2);
//
// 			box->actx[i] /= norm2;
// 			box->acty[i] /= norm2;
// 			box->actz[i] /= norm2;

			box->oldactx[i] = box->actx[i];
			box->oldacty[i] = box->acty[i];
			box->oldactz[i] = box->actz[i];

			// body-frame patch coordinates -- radius agnostic and immutable
			box->patchposx[chromindex-1][0]++;
			box->patchposy[chromindex-1][0]++;
			box->patchposz[chromindex-1][0]++;
			box->patchposx[chromindex-1][j] = spampx;
			box->patchposy[chromindex-1][j] = spampy;
			box->patchposz[chromindex-1][j] = spampz;

			box->Gpatchposx[chromindex-1][j] = spamgx;
			box->Gpatchposy[chromindex-1][j] = spamgy;
			box->Gpatchposz[chromindex-1][j] = spamgz;



			box->radius[i] = rad;
			box->mass[i] = mass;
			box->invmass[i] = 1.0/mass;
			box->k_spring[i] = p->spring_scale*k_spring;
			box->activity[i] = act*p->actscale;
			box->lamin[i] = p->lamin_scale*lamin;
			box->isBond[i] = 1.0; // NOTE: very bad idea but placeheld to later be replaced with a 2D array

		}
	}

// 	printf("finding axes for init ori\n");
// 	findAxesInitOri(p,box);
	// now rotate to find global coordinates of each patch
	printf("found axes, rotating patches\n");
	for(ii=0;ii<p->nChrom;ii++)
	{
		// set quarternion for ith chromosome based on axis and angle
		i = box->chromList[ii];
//		printf("chrom %d par %d\n",ii,i);
		box->qx[i] = box->axisx[i]*sin(0.5*acos(box->cosRotAngle[i]));
		box->qy[i] = box->axisy[i]*sin(0.5*acos(box->cosRotAngle[i]));
		box->qz[i] = box->axisz[i]*sin(0.5*acos(box->cosRotAngle[i]));
		box->qw[i] = cos(0.5*acos(box->cosRotAngle[i]));

		box->oldqx[i] = box->qx[i];
		box->oldqy[i] = box->qy[i];
		box->oldqz[i] = box->qz[i];
		box->oldqw[i] = box->qw[i];

		quatnorm = sqrt(box->oldqx[i]*box->oldqx[i] + box->oldqy[i]*box->oldqy[i] + box->oldqz[i]*box->oldqz[i] + box->oldqw[i]*box->oldqw[i]);

		box->oldqx[i] /= quatnorm;
		box->oldqy[i] /= quatnorm;
		box->oldqz[i] /= quatnorm;
		box->oldqw[i] /= quatnorm;
// 		printf("set quat from rotations chrom %d par %d\n",ii,i);
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
			j = box->patchList[ii][jj];
			rotate(box->patchposx[ii][jj],box->patchposy[ii][jj],box->patchposz[ii][jj],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],rotpatch);
			rotpatchnorm=sqrt(rotpatch[0]*rotpatch[0] + rotpatch[1]*rotpatch[1] + rotpatch[2]*rotpatch[2]);

			// test global position vs local position

			spamdiff = box->Gpatchposx[ii][jj] - (box->X[i] + box->radius[i]*rotpatch[0]/rotpatchnorm);
			spamdiff += box->Gpatchposy[ii][jj] - (box->Y[i] + box->radius[i]*rotpatch[1]/rotpatchnorm);
			spamdiff += box->Gpatchposz[ii][jj] - (box->Z[i] + box->radius[i]*rotpatch[2]/rotpatchnorm);

			if (fabs(spamdiff) > 0.01)
			{
				printf("chr %d patch %d raw %d not rotated properly %4.15f %4.15f %4.15f -- %4.15f %4.15f %4.15f\n", \
				ii,jj,j, box->Gpatchposx[ii][jj],box->Gpatchposy[ii][jj],box->Gpatchposz[ii][jj], \
				box->X[i] + box->radius[i]*rotpatch[0]/rotpatchnorm, box->Y[i] + box->radius[i]*rotpatch[1]/rotpatchnorm,\
				box->Z[i] + box->radius[i]*rotpatch[2]/rotpatchnorm);

				box->Gpatchposx[ii][jj] = box->X[i] + box->radius[i]*rotpatch[0]/rotpatchnorm; // position of patches
				box->Gpatchposy[ii][jj] = box->Y[i] + box->radius[i]*rotpatch[1]/rotpatchnorm;
				box->Gpatchposz[ii][jj] = box->Z[i] + box->radius[i]*rotpatch[2]/rotpatchnorm;
// 				exit(0);
			}


// 			// now assign the global positions of the patch positions -- COM +
// 			box->Gpatchposx[ii][jj] = box->X[i] + box->radius[i]*rotpatch[0]/rotpatchnorm; // position of patches
// 			box->Gpatchposy[ii][jj] = box->Y[i] + box->radius[i]*rotpatch[1]/rotpatchnorm;
// 			box->Gpatchposz[ii][jj] = box->Z[i] + box->radius[i]*rotpatch[2]/rotpatchnorm;
//
			// for now the real positions are the patch positions + some small displacement
// 			spamindex = box->patchList[ii][j];
// 			box->X[spamindex] = box->Gpatchposx[ii][jj] + 0.5*box->radius[i]*(2.0*drand48()-1.0); // position of connected bead
// 			box->Y[spamindex] = box->Gpatchposy[ii][jj] + 0.5*box->radius[i]*(2.0*drand48()-1.0);
// 			box->Z[spamindex] = box->Gpatchposz[ii][jj] + 0.5*box->radius[i]*(2.0*drand48()-1.0);

			box->oldX[spamindex] = box->X[spamindex];
			box->oldY[spamindex] = box->Y[spamindex];
			box->oldZ[spamindex] = box->Z[spamindex];

// 			printf("ch %d patch %d (real %d %d) Init patch pos %f %f %f -- bead pos %f %f %f\n",ii,j,i,spamindex,box->Gpatchposx[ii][j],box->Gpatchposy[ii][j],box->Gpatchposz[ii][j],box->oldX[spamindex],box->oldY[spamindex],box->oldZ[spamindex]);

		}
	}
	return mdstep;

}



void writeSimRestart(FILE *fp,struct PARAINPUT *p, struct SYSTEM *box,int spamstep)
{


	// FILE to write down the restart info needed to respawn a simulation --
	// note current geometry, stress and a other sim state are stored in a sim state file -- different from this
	// likewise, the RNG state is stored in a third file

	int i,j,ii,jj,spamindex;
	double spamaxis[3];
	double spamact=0.0,spamlam=0.0;


	// write the box dimensions and stress state
	fprintf(fp,"%d %f %f %f %f %f %f %f %f\n",spamstep,box->lx,p->ella,p->ellb,p->ellc,box->currsig0x,box->currsig0z,box->sigtheta,box->prevArea);

	for(ii=0;ii<p->nChrom;ii++)
	{
		i=box->chromList[ii];

		box->cosRotAngle[i] = findAxesQuat(box->qx[i],box->qy[i],box->qz[i],box->qw[i],spamaxis);

        box->axisx[i] = spamaxis[0];
        box->axisy[i] = spamaxis[1];
        box->axisz[i] = spamaxis[2];



        fprintf(fp,"%d %d %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %f %f %f %f %f\n", \
        ii+1,box->type[i],box->X[i],box->Y[i],box->Z[i],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],0.0,0.0,0.0,0.0,0.0,0.0, \
        box->radius[i],box->mass[i],0.0,box->activity[i],box->lamin[i]);
//              printf("printing chrom info no issue %d -- %d\n",ii,structbox.patchList[ii][0]);
        for(jj=1;jj<=box->patchList[ii][0];jj++)
        {
			j=box->patchList[ii][jj];

			// write patchposxyz and Gpatchposxyz -- when reading, test compare (anchor pos) then write XYZ (extruded patch pos)
			spamact=0;
			if(p->actscale > 0)
			spamact = box->activity[j]/p->actscale;
			else
			spamact = 0.0;

			spamlam=0;
			if(p->lamin_scale > 0)
			spamlam =box->lamin[j]/p->lamin_scale;
			else
			spamlam = 0.0;

			fprintf(fp,"%d %d %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %4.12f %f %f %f %f %f\n", \
			ii+1,0,box->X[j],box->Y[j],box->Z[j],box->actx[j],box->acty[j],box->actz[j],0.0,box->patchposx[ii][jj],box->patchposy[ii][jj], \
			box->patchposz[ii][jj],box->Gpatchposx[ii][jj],box->Gpatchposy[ii][jj],box->Gpatchposz[ii][jj],box->radius[j],box->mass[j], \
			box->k_spring[j]/p->spring_scale,spamact,spamlam);

        }
	}

}



void writeAnadump(struct PARAINPUT *p, struct SYSTEM *box, struct APPOGGIO *a, int MCstep)
{
	// subroutine to write output in format that mimics LAMMPS
	// this is so that AnaCode can be used on configs printed out of this code as well
	int i,j,ii,jj;
	double spamx,spamy,spamz;
	double rotpatch[3];

	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)




	printf("writing to dump file\n");
	fprintf(box->fpLAMMPSdump,"ITEM: TIMESTEP\n");
	fprintf(box->fpLAMMPSdump,"%d\n",MCstep);		// TIMESTEP
	fprintf(box->fpLAMMPSdump,"ITEM: NUMBER OF ATOMS\n");
	fprintf(box->fpLAMMPSdump,"%d\n",2*p->nPatchTot+p->nChrom + p->nSurfPoints);
	fprintf(box->fpLAMMPSdump,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(box->fpLAMMPSdump,"%f %f\n",-p->ella,p->ella);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-p->ellb,p->ellb);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-p->ellc,p->ellc);
	fprintf(box->fpLAMMPSdump,"ITEM: ATOMS id type radius Transparency xu yu zu vx vy vz ox oy oz\n");


	for(i=0;i<p->nSurfPoints;i++)
	{
		if(box->surfActive[i]==1)
		fprintf(box->fpLAMMPSdump,"%d %d %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,p->nAto+4,0.5*avesep,0.85*(1.0 - box->isMembIntermingling[i]),box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],0.0,0.0,0.0,0.0,0.0,0.0);
		if(box->surfActive[i]==0)
			fprintf(box->fpLAMMPSdump,"%d %d %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,p->nAto+5,0.5*avesep,0.85*(1.0 - box->isMembIntermingling[i]),box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],0.0,0.0,0.0,0.0,0.0,0.0);
// 		fprintf(box->fpSurfPoints,"%d %d %f %lf %lf %lf\n",i,200,0.5*avesep,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}

	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		fprintf(box->fpLAMMPSdump,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,ii+1,box->patchList[ii][0],box->type[i],box->radius[i],box->X[i],box->Y[i],box->Z[i],box->vx[i],box->vy[i],box->vz[i],box->axisx[i],box->axisy[i],box->axisz[i]);
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
// 			rotate(box->patchposx[i][j],box->patchposy[i][j],box->patchposz[i][j],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],rotpatch);
// // 			spamx = box->X[i] + box->patchposx[i][j];
// // 			spamy = box->Y[i] + box->patchposy[i][j];
// // 			spamz = box->Z[i] + box->patchposz[i][j];
//
// 			spamx = box->X[i] + rotpatch[0];
// 			spamy = box->Y[i] + rotpatch[1];
// 			spamz = box->Z[i] + rotpatch[2];

            spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
            // patches only have positions

			if(box->patchPointer[ii][jj]==p->nAto+2)
			{
			fprintf(box->fpLAMMPSdump,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,jj,box->patchPointer[ii][jj]+1,0.5*box->radius[j],spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
            // beads have positions and velocities

            fprintf(box->fpLAMMPSdump,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,jj,box->patchPointer[ii][jj]+1,box->radius[j],box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
			else
			{
			fprintf(box->fpLAMMPSdump,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,jj,box->patchPointer[ii][jj]+1,0.5*box->radius[j],spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
            // beads have positions and velocities

            fprintf(box->fpLAMMPSdump,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,jj,box->patchPointer[ii][jj]+1,box->radius[j],box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
		}

	}
}

void writeFulldump(struct PARAINPUT *p, struct SYSTEM *box, struct APPOGGIO *a, int MCstep)
{
	// subroutine to write output in format that mimics LAMMPS
	// this is so that AnaCode can be used on configs printed out of this code as well
	int i,j,ii,jj;
	double spamx,spamy,spamz;
	double rotpatch[3];

	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)




	printf("writing to dump file\n");
	fprintf(box->fpLAMMPSdump,"ITEM: TIMESTEP\n");
	fprintf(box->fpLAMMPSdump,"%d\n",MCstep);		// TIMESTEP
	fprintf(box->fpLAMMPSdump,"ITEM: NUMBER OF ATOMS\n");
	fprintf(box->fpLAMMPSdump,"%d\n",2*p->nPatchTot+p->nChrom + p->nSurfPoints);
	fprintf(box->fpLAMMPSdump,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->lx,0.5*box->lx);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->ly,0.5*box->ly);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->lz,0.5*box->lz);
	fprintf(box->fpLAMMPSdump,"ITEM: ATOMS id type pointer radius Transparency xu yu zu vx vy vz ox oy oz\n");


	for(i=0;i<p->nSurfPoints;i++)
	{
		if(box->surfActive[i]==1)
		fprintf(box->fpLAMMPSdump,"%d %d %d %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,p->nAto+4,p->nAto+4,0.5*avesep,0.85*(1.0 - box->isMembIntermingling[i]),box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],box->surffx[i],box->surffy[i],box->surffz[i],0.0,0.0,0.0);
		if(box->surfActive[i]==0)
		fprintf(box->fpLAMMPSdump,"%d %d %d %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,p->nAto+5,p->nAto+5,0.5*avesep,0.85,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i],box->surffx[i],box->surffy[i],box->surffz[i],0.0,0.0,0.0);
// 		fprintf(box->fpSurfPoints,"%d %d %f %lf %lf %lf\n",i,200,0.5*avesep,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}

	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
// 			rotate(box->patchposx[i][j],box->patchposy[i][j],box->patchposz[i][j],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],rotpatch);
// // 			spamx = box->X[i] + box->patchposx[i][j];
// // 			spamy = box->Y[i] + box->patchposy[i][j];
// // 			spamz = box->Z[i] + box->patchposz[i][j];
//
// 			spamx = box->X[i] + rotpatch[0];
// 			spamy = box->Y[i] + rotpatch[1];
// 			spamz = box->Z[i] + rotpatch[2];

            spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
            // patches only have positions

			if(box->patchPointer[ii][jj]==p->nAto+2)
			{
			fprintf(box->fpLAMMPSdump,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,box->patchPointer[ii][jj]+1,0.5*box->radius[j],0.95,spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
            // beads have positions and velocities

            fprintf(box->fpLAMMPSdump,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,box->patchPointer[ii][jj]+1,box->radius[j],0.95,box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
			else
			{
			fprintf(box->fpLAMMPSdump,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,box->patchPointer[ii][jj]+1,0.5*box->radius[j],0.8*(1.0-box->isIntermingling[j]),spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
            // beads have positions and velocities

            fprintf(box->fpLAMMPSdump,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",j,ii+1,box->patchPointer[ii][jj]+1,box->radius[j],0.8*(1.0- box->isIntermingling[j]),box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
		}
		fprintf(box->fpLAMMPSdump,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,ii+1,box->type[i],box->radius[i],0.65,box->X[i],box->Y[i],box->Z[i],box->vx[i],box->vy[i],box->vz[i],box->axisx[i],box->axisy[i],box->axisz[i]);
	}
}

void readSurfPoints(struct PARAINPUT *p, struct SYSTEM *box,int step)
{
	int i,j,n,spamt2;
	double spamt1,surfrad,spamtransp,spamx,spamy,spamz;
	char datoSTR[400];


	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)



	fscanf(box->fpSurfPoints,"%s", datoSTR); 		// ITEM:
	fscanf(box->fpSurfPoints, "%s", datoSTR);		// TIMESTEP
	fscanf(box->fpSurfPoints, "%s", datoSTR);
	n=atoi(datoSTR);
	//printf("time actually, not number of atoms %d\n",n);
	if(n!=step)
	{
		printf("Error!\nI am reading surf file %d instead of %d\n",n,step);
		exit(0);
	}

	for(n=0;n<5;n++)				// ITEM: NUMBER OF ATOMS
	{
		fscanf(box->fpSurfPoints, "%s", datoSTR);
	}
	n=atoi(datoSTR);
	printf("is this number of atoms?? %d -%d \n",n,2*p->nSurfPoints);

	if(n!=2*p->nSurfPoints)
	{
		printf("Error!\nI am reading %d atoms instead of %d\n",n,2*p->nSurfPoints);
		exit(0);
	}
	for(n=0;n<6;n++)				// ITEM: BOX BOUNDS pp pp pp
	{
		fscanf(box->fpSurfPoints, "%s", datoSTR);
		printf("%s, wtf\n",datoSTR);
	}
	for(n=0;n<3;n++)
	{
		fscanf(box->fpSurfPoints, "%s", datoSTR);
		spamx=atof(datoSTR);
		fscanf(box->fpSurfPoints, "%s", datoSTR);
		spamx-=atof(datoSTR);
		spamx*=-1;
		printf("%f\n",spamx);
		//printf("the box dimension of %d is %f\n",n,d->lBox[n]);
	}
//
	if(spamx != box->lx)
	{
		printf("%f box size mismatch %f\n",spamx,box->lx);
	}

	printf("entering loop\n");

	for(n=0;n<10;n++) // ITEM: ATOMS id type color radius Transparency xu yu zu
	{
		fscanf(box->fpSurfPoints, "%s",datoSTR);
// 		printf("%s\n",datoSTR);
	}
// 	printf("nAto %d -- vz %s",s->Ato,datoSTR);
	for(n=0;n<p->nSurfPoints;n++) // xu yu and zu
	{
		fscanf(box->fpSurfPoints,"%d %lf %d %lf %lf %lf %lf %lf\n",&i,&spamt1,&spamt2,&surfrad,&spamtransp,&spamx,&spamy,&spamz);
// 		printf("Fund -- %d %d %f %f %f %f %f\n",i,n,surfrad,spamtransp,spamx,spamy,spamz);
		if(i!=n)
		{
			printf("surf Fund index mismatch at %d -- should be %d -- %f\n",i,n,spamtransp);
			exit(0);
		}
		box->surfFundx[n] = spamx;
		box->surfFundy[n] = spamy;
		box->surfFundz[n] = spamz;
		box->surfActive[n] = spamt2;
		if(fabs(surfrad - 0.5*avesep) > 0.001)
		{
			printf("%d surf point radius mismatch %f should be %f\n",n,surfrad,0.5*avesep);
			exit(0);
		}

	}

	for(n=0;n<p->nSurfPoints;n++) // xu yu and zu
	{
		fscanf(box->fpSurfPoints,"%d %lf %d %lf %lf %lf %lf %lf\n",&i,&spamt1,&spamt2,&surfrad,&spamtransp,&spamx,&spamy,&spamz);
// 		printf("Point -- %d %d %f %f %f %f %f\n",i,n,surfrad,spamtransp,spamx,spamy,spamz);
		if(i!=n)
		{
			printf("surf point index mismatch at %d -- should be %d -- %f\n",i,n,spamtransp);
			exit(0);
		}
		if(fabs(surfrad - 0.5*avesep) > 0.001)
		{
			printf("%d surf point radius mismatch %f should be %f\n",n,surfrad,0.5*avesep);
			exit(0);
		}
		box->surfPointx[n] = spamx;
		box->surfPointy[n] = spamy;
		box->surfPointz[n] = spamz;

	}

}


void writeSurfPoints(struct PARAINPUT *p, struct SYSTEM *box,int MCstep)
{
	int i,j;
	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	fprintf(box->fpSurfPoints,"ITEM: TIMESTEP\n");
	fprintf(box->fpSurfPoints,"%d\n",MCstep);		// TIMESTEP
	fprintf(box->fpSurfPoints,"ITEM: NUMBER OF ATOMS\n");
	fprintf(box->fpSurfPoints,"%d\n",2*p->nSurfPoints);
	fprintf(box->fpSurfPoints,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(box->fpSurfPoints,"%f %f\n",-0.5*box->lx,0.5*box->lx);
	fprintf(box->fpSurfPoints,"%f %f\n",-0.5*box->ly,0.5*box->ly);
	fprintf(box->fpSurfPoints,"%f %f\n",-0.5*box->lz,0.5*box->lz);
	fprintf(box->fpSurfPoints,"ITEM: ATOMS id type color radius Transparency xu yu zu\n");
	for(i=0;i<p->nSurfPoints;i++)
	{
// 		fprintf(box->fpSurfPoints,"%d %d %f %lf %lf %lf\n",i,200,0.5*avesep,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i]);
		if(box->surfActive[i]==1)
		fprintf(box->fpSurfPoints,"%d %d %d %f %f %lf %lf %lf\n",i,100,box->surfActive[i],0.5*avesep,0.1,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
		if(box->surfActive[i]==0)
		fprintf(box->fpSurfPoints,"%d %d %d %f %f %lf %lf %lf\n",i,300,box->surfActive[i],0.5*avesep,0.1,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}
	for(i=0;i<p->nSurfPoints;i++)
	{
		if(box->surfActive[i]==1)
		fprintf(box->fpSurfPoints,"%d %d %d %f %f %lf %lf %lf\n",i,200,box->surfActive[i],0.5*avesep,0.85,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i]);
		if(box->surfActive[i]==0)
		fprintf(box->fpSurfPoints,"%d %d %d %f %f %lf %lf %lf\n",i,400,box->surfActive[i],0.5*avesep,0.85,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i]);
// 		fprintf(box->fpSurfPoints,"%d %d %f %lf %lf %lf\n",i,200,0.5*avesep,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}


}

void readCurrHiC(struct PARAINPUT *p, struct SYSTEM *box)
{
	int i,j;
	int nhap;
	double spam;
	nhap = (int) (0.5*p->nChrom);

	while(!feof(box->fpHiC))
	{
		fscanf(box->fpHiC,"%d %d %lf\n",&i,&j,&spam);
		box->HiCMat[i][j] = spam;

	}

}


void writeCurrHiC(struct PARAINPUT *p, struct SYSTEM *box)
{

	int i,j;
	int nhap;

	nhap = (int) (0.5*p->nChrom);
	for(i=0;i<nhap;i++)
	{
		for(j=0;j<nhap;j++)
		{
			fprintf(box->fpHiC,"%d %d %f\n",i,j,box->HiCMat[i][j]);
		}
// 		fprintf(box->fpHiC,"\n");
	}
}

void readExpectedMat(struct PARAINPUT *p, struct SYSTEM *box)
{
	int i;
	double spam;
	int nhap = (int)(0.5*p->nChrom);
	int nconstraints = (nhap + 1)*nhap - (int) (nhap*(nhap-1)/2) + nhap+1 - nhap;
	fscanf(box->fpIPD,"%lf ",&spam);
	for(i=0;i<nconstraints;i++)
	{
		fscanf(box->fpIPD,"%lf ",&spam);
		box->expectedMat[i] = spam;
// 		printf("%d %f %f\n",i,box->expectedMat[i],spam);

	}

}

void writeConstraints(struct PARAINPUT *p, struct SYSTEM *box,int startflag)
{
	int i;

	int nhap = (int)(0.5*p->nChrom);
	int nconstraints = (nhap + 1)*nhap - (int) (nhap*(nhap-1)/2) + nhap+1 - nhap;
	printf("writing constraints!!\n");

	if(startflag==0)
	{
		fprintf(box->fpIPD,"%d ",startflag-1);
		for(i=0;i<nconstraints;i++)
		{
			fprintf(box->fpIPD,"%f ",box->expectedMat[i]);
		}
		fprintf(box->fpIPD,"\n");
		fprintf(box->fpIPD,"%d ",startflag);
		for(i=0;i<nconstraints;i++)
		{
			fprintf(box->fpIPD,"%f ",box->constraintMat[i]);
		}
		fprintf(box->fpIPD,"\n");
	}

	else
	{
		fprintf(box->fpIPD,"%d ",startflag);
		for(i=0;i<nconstraints;i++)
		{
			fprintf(box->fpIPD,"%f ",box->constraintMat[i]);
		}
		fprintf(box->fpIPD,"\n");
	}

}

void dropSurf(FILE *fp,struct PARAINPUT *p, struct SYSTEM *box)
{
	int i;
	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	for(i=0;i<p->nSurfBasicPoints;i++)
	{
// 		fprintf(box->fpSurfPoints,"%d %d %f %lf %lf %lf\n",i,200,0.5*avesep,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i]);
		fprintf(fp,"%d %lf %lf %lf %lf\n",i,0.5*avesep,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}
}


void readSurf(FILE *fp,struct PARAINPUT *p, struct SYSTEM *box)
{
	int i,j,n;
	double spamrad,spamx,spamy,spamz;
	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)

	printf("%d surf points to be read -- %d %d, sep by %f\n",p->nSurfBasicPoints,p->nSurfCycPoints,p->nSurfPoints,0.5*avesep);
	for(i=0;i<p->nSurfBasicPoints;i++)
	{
		// surfid, type, transparency, radius, x,y,z
		fscanf(fp,"%d %lf %lf %lf %lf\n",&n,&spamrad,&spamx,&spamy,&spamz);

		box->surfFundx[i] = spamx;
		box->surfFundy[i] = spamy;
		box->surfFundz[i] = spamz;

// 		printf("%d %lf %lf %lf %lf\n",n,spamrad,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);

		if(fabs(0.5*avesep - spamrad) > 0.0001)
		{
			printf("surf point size mismatch -- please fix %2.12f should be %2.12f\n",spamrad,0.5*avesep);
// 			exit(0);
		}

		box->surfActive[i]=1;


	}
	printf("read points checking seps\n");
	/// check minsep again
	double delx,dely,delz;
	double spamdist;
	box->minsep=p->ellai;
	box->maxOffEll = 0.0;
	for(i=0;i<p->nSurfBasicPoints-1;i++)
	{
		for(j=i+1;j<p->nSurfBasicPoints;j++)
		{
			delx = box->surfFundx[j] - box->surfFundx[i];
			dely = box->surfFundy[j] - box->surfFundy[i];
			delz = box->surfFundz[j] - box->surfFundz[i];

			spamdist = sqrt(delx*delx + dely*dely + delz*delz);
// 					rij = sqrt(rijsq);

			if(spamdist==0)
			{
				printf("%d %d ij %f %f %f -- %f %f %f\n",i,j,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i],box->surfFundx[j],box->surfFundy[j],box->surfFundz[j]);
			}

			if(spamdist < box->minsep)
			{
				box->minsep = spamdist;
			}
			if(box->minsep > 0.49*avesep)
			{
				printf("initial min sep greater than  %f -- %f surfrad %f\n",0.49*avesep,box->minsep,avesep);
				box->minsep = 0.49*avesep;
			}
		}


		delx = box->surfFundx[i];
		dely = box->surfFundy[i];
		delz = box->surfFundz[i];
	// 		printf("surf %f %f %f -- fund %f %f %f\n",spamx,spamy,spamz,fundx,fundy,fundz);

		spamdist = fabs(sqr(delx/p->ellai) + sqr(dely/p->ellbi) + sqr(delz/p->ellci) - 1);
		if(spamdist > sqr(box->maxOffEll))
			box->maxOffEll = sqrt(spamdist);


	}

	printf("min sep at read -- %f --- max off ell at read %f\n",box->minsep,box->maxOffEll);

}

void readCyc(FILE *fp,struct PARAINPUT *p, struct SYSTEM *box)
{
	int i,j,n,c1,c2,c3;
	double spamrad,spamx,spamy,spamz;
	char spamstring[500];

	printf("%d surf points to be read -- %d %d\n",p->nSurfBasicPoints,p->nSurfCycPoints,p->nSurfPoints);
	for(i=p->nSurfBasicPoints;i<p->nSurfPoints;i++)
	{
// 		printf("%d -- %d index %f %f %f %d %d %d\n",i,i-p->nSurfBasicPoints,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i], \
			   box->surfcycpt[i-p->nSurfBasicPoints][0],box->surfcycpt[i-p->nSurfBasicPoints][1],box->surfcycpt[i-p->nSurfBasicPoints][2]);
// 		fscanf(fp,"%s\n",spamstring);

		fscanf(fp,"%d %lf %lf %lf %lf %d %d %d\n",&n,&spamrad,&spamx,&spamy,&spamz,&c1,&c2,&c3);

		box->surfFundx[i] = spamx;
		box->surfFundy[i] = spamy;
		box->surfFundz[i] = spamz;

		box->surfcycpt[n][0] = c1;
		box->surfcycpt[n][1] = c2;
		box->surfcycpt[n][2] = c3;

//
// 		printf("%d %d %f %f %f %f %d %d %d\n",i,n,spamrad,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i], \
		box->surfcycpt[n][0],box->surfcycpt[n][1],box->surfcycpt[n][2]);
// 		printf("%s\n",spamstring);
	}

	printf("done reading cyc\n");
	box->activeAccum = 0;
}







