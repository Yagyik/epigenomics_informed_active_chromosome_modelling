double gauss_var_AT(double mu,double sig)
{
    double gauss_r,v1,v2,gauss_var;
    gauss_r=2.0;
    while(gauss_r>=1.0)
    {
//         v1=2.0*genrand_real1(&(*a))-1.0;
//         v2=2.0*genrand_real1(&(*a))-1.0;
		v1=2.0*drand48()-1.0;
		v2=2.0*drand48()-1.0;
//         printf("v1 v2 %f %f\n",v1,v2);
        gauss_r = v1*v1 + v2*v2;
    }
//     printf("gr log gr %f %f\n",gauss_r,log(gauss_r));
    gauss_var=v1*sqrt(-2.0*log(gauss_r)/gauss_r);
    gauss_var=mu+sig*gauss_var;
    
    return gauss_var;
}


int init_vel(struct PARAINPUT *p,struct SYSTEM *box)
{
    int i;
    double vxcom=0.0, vycom=0.0,vzcom=0.0;
    double sig=sqrt(p->tem);
    for(i=0;i<p->nAto;i++)
    {
        box->vx[i]=gauss_var_AT(0,sig); ///box->mass[i]; // init momentum and get vel
        box->vy[i]=gauss_var_AT(0,sig); ///box->mass[i];
        box->vz[i]=gauss_var_AT(0,sig); ///box->mass[i];
        //printf("%f %f %f\n",box->vx[i],box->vy[i],box->vz[i]);
        vxcom+=box->vx[i];
        vycom+=box->vy[i];
        vzcom+=box->vz[i];
    }
    
    for(i=0;i<p->nAto;i++)
    {
        box->vx[i]-=vxcom/p->nAto;
        box->vy[i]-=vycom/p->nAto;
        box->vz[i]-=vzcom/p->nAto;
//         printf("velocities %f %f %f\n",box->vx[i],box->vy[i],box->vz[i]);
    }
    
    box->KE = 0.0; // reset for call on this time step
    for(i=0;i<p->nAto;i++)
    {
        box->KE += 0.5*(box->vx[i]*box->vx[i] + box->vy[i]*box->vy[i] + box->vz[i]*box->vz[i]);
    }
    box->T_meas=box->KE*2.0/(3.0*p->nAto); // actually the measured temp
// 	p->tem=1.0;
	printf("last norm for vel %f %f and %d atoms\n",box->T_meas,p->tem,p->nAto);
    for(i=0;i<p->nAto;i++)
    {
        box->vx[i]/=sqrt(box->T_meas/p->tem);
        box->vy[i]/=sqrt(box->T_meas/p->tem);
		box->vz[i]/=sqrt(box->T_meas/p->tem);
		
// 		box->oldvx[i] = box->vz[i];
// 		box->oldvy[i] = box->vy[i];
// 		box->oldvz[i] = box->vz[i];
    }
    printf("initialised velocities\n");
}


// this subroutine integrates the raw positions of the chromosomes, the orientation of the chromosomes and then positions of patches in global coordinate frame

void updateTrial(struct SYSTEM *box, struct PARAINPUT *p,int chainflag)
{
	int ch,i,n1,ii,j,jj;
	double diffX,diffY,diffZ,rij,rijsq,r0,rgori;
    double qdot[4];
    double spamaxis[3];
    double rotpatch[3];
	double rotpatchnorm;
	double spamwx,spamwy,spamwz;
	double quatnorm;
    double actnorm;

    if(chainflag==1)
        printf("attempting with new dt %d %f\n",chainflag,p->dt);
	for(i=0;i<p->nAto;i++)
	{
       
//         if(box->type[i]==0)  // bead random forces
// 		{
// 			box->fRx[i] = sqrt(2.0*p->thermostatFact*p->tem/(6.0*M_PI*p->dt*p->gamma*box->radius[i]))*box->rgx[i];
// 			box->fRy[i] = sqrt(2.0*p->thermostatFact*p->tem/(6.0*M_PI*p->dt*p->gamma*box->radius[i]))*box->rgy[i];
// 			box->fRz[i] = sqrt(2.0*p->thermostatFact*p->tem/(6.0*M_PI*p->dt*p->gamma*box->radius[i]))*box->rgz[i];
//
// 		}
// 		else // chromosome random forces
// 		{
// 			box->fRx[i] = sqrt(2.0*p->tem/(6.0*M_PI*p->dt*p->gamma*box->radius[i]))*box->rgx[i];
// 			box->fRy[i] = sqrt(2.0*p->tem/(6.0*M_PI*p->dt*p->gamma*box->radius[i]))*box->rgy[i];
// 			box->fRz[i] = sqrt(2.0*p->tem/(6.0*M_PI*p->dt*p->gamma*box->radius[i]))*box->rgz[i];
// 		}
			
        // old vectors are unit vectors (Gaussian noise vectors need not be)
// 		box->actx[i] = box->oldactx[i] + sqrt(2.0/(p->dt*p->taup))*(box->oldacty[i]*box->argz[i] - box->oldactz[i]*box->argy[i]);
//         box->acty[i] = box->oldacty[i] + sqrt(2.0/(p->dt*p->taup))*(box->oldactz[i]*box->argx[i] - box->oldactx[i]*box->argz[i]);
//         box->actz[i] = box->oldactz[i] + sqrt(2.0/(p->dt*p->taup))*(box->oldactx[i]*box->argy[i] - box->oldacty[i]*box->argx[i]);


        if(box->type[i]==0)  // bead random forces
		{
			box->fRx[i] = sqrt(2.0*p->thermostatFact*p->tem*p->dt/(6.0*M_PI*p->gamma*box->radius[i]))*box->rgx[i];
			box->fRy[i] = sqrt(2.0*p->thermostatFact*p->tem*p->dt/(6.0*M_PI*p->gamma*box->radius[i]))*box->rgy[i];
			box->fRz[i] = sqrt(2.0*p->thermostatFact*p->tem*p->dt/(6.0*M_PI*p->gamma*box->radius[i]))*box->rgz[i];

		}
		else // chromosome random forces
		{
			box->fRx[i] = sqrt(2.0*p->tem*p->dt/(6.0*M_PI*p->gamma*box->radius[i]))*box->rgx[i];
			box->fRy[i] = sqrt(2.0*p->tem*p->dt/(6.0*M_PI*p->gamma*box->radius[i]))*box->rgy[i];
			box->fRz[i] = sqrt(2.0*p->tem*p->dt/(6.0*M_PI*p->gamma*box->radius[i]))*box->rgz[i];
		}
        box->actx[i] = box->oldactx[i] + sqrt(2.0*p->dt/p->taup)*(box->oldacty[i]*box->argz[i] - box->oldactz[i]*box->argy[i]);
        box->acty[i] = box->oldacty[i] + sqrt(2.0*p->dt/p->taup)*(box->oldactz[i]*box->argx[i] - box->oldactx[i]*box->argz[i]);
        box->actz[i] = box->oldactz[i] + sqrt(2.0*p->dt/p->taup)*(box->oldactx[i]*box->argy[i] - box->oldacty[i]*box->argx[i]);
//         printf("%f old %f\n",box->actx[i],box->oldactx[i]);

        if(box->actx[i]*box->actx[i] + box->acty[i]*box->acty[i] + box->actz[i]*box->actz[i] > 0.0)
        actnorm = sqrt(box->actx[i]*box->actx[i] + box->acty[i]*box->acty[i] + box->actz[i]*box->actz[i]);
        else
        actnorm = 1.0;

        // make sure these are unit vectors
        box->actx[i] /= actnorm;
        box->acty[i] /= actnorm;
        box->actz[i] /= actnorm;


		// bead positions and chromosome COMs
		box->trialX[i] = box->oldX[i] + (p->dt*box->fx[i])/(6.0*M_PI*p->gamma*box->radius[i])  + box->fRx[i] + p->dt*box->activity[i]*box->actx[i];
		box->trialY[i] = box->oldY[i] + (p->dt*box->fy[i])/(6.0*M_PI*p->gamma*box->radius[i])  + box->fRy[i] + p->dt*box->activity[i]*box->acty[i];
		box->trialZ[i] = box->oldZ[i] + (p->dt*box->fz[i])/(6.0*M_PI*p->gamma*box->radius[i])  + box->fRz[i] + p->dt*box->activity[i]*box->actz[i];
		
// 		printf("%d position update from %f %f %f to %f %f %f\n",i,box->oldX[i],box->oldY[i],box->oldZ[i],box->trialX[i],box->trialY[i],box->trialZ[i]);
// 		printf("%d forces %f %f %f random forces are %f %f %f %f\n",i,box->fx[i],box->fy[i],box->fz[i],box->fRx[i],box->fRy[i],box->fRz[i],p->gamma);
		
		
		
		
		// rotational velocity update -- need to make this array later
// 		rgori=gauss_var_AT(0,1);		
// 		box->vori[i] = sqrt(1-p->em2dttp)*rgori;
// 		// remaining half step update positions		
// 		box->ori[i] = box->ori[i] + 0.5*p->dt*box->vori[i];
        
    }
    if(chainflag==1)
        printf("success disp with new dt %d %f\n",chainflag,p->dt);
     
    for(ii=0;ii<p->nChrom;ii++)
    {
        i=box->chromList[ii];
        // update trial angular velocities
        
        spamwx = box->torquex[i]/(6.0*M_PI*p->rotgamma*box->radius[i]); // plus some noise term for the Brownian rotation case
        spamwy = box->torquey[i]/(6.0*M_PI*p->rotgamma*box->radius[i]); //
        spamwz = box->torquez[i]/(6.0*M_PI*p->rotgamma*box->radius[i]);
//         printf("old q %f %f %f %f and ang vel %2.12f %2.12f %2.12f and torque %f %f %f\n",box->oldqx[i],box->oldqy[i],box->oldqz[i],box->oldqw[i],spamwx,spamwy,spamwz,box->torquex[i],box->torquey[i],box->torquez[i]);
        // get qdot from angular velocities in global reference frame
// 		printf("old q pre qdot %f %f %f %f\n",box->oldqx[i],box->oldqy[i],box->oldqz[i],box->oldqw[i]);
        findQDotG(box->oldqx[i],box->oldqy[i],box->oldqz[i],box->oldqw[i],spamwx,spamwy,spamwz,qdot);
//         printf("%d %d qdot coming out to be %2.12f %2.12f %2.12f %2.12f\n",ii,i,qdot[0],qdot[1],qdot[2],qdot[3]);
// 		if(qdot[0]>1000)
// 			exit(0);
// 		printf("old q post qdot %f %f %f %f\n",box->oldqx[i],box->oldqy[i],box->oldqz[i],box->oldqw[i]);
        // update the quaternions
        box->trialqx[i] = box->oldqx[i] + p->dt*qdot[0];
        box->trialqy[i] = box->oldqy[i] + p->dt*qdot[1];
        box->trialqz[i] = box->oldqz[i] + p->dt*qdot[2];
        box->trialqw[i] = box->oldqw[i] + p->dt*qdot[3];
        quatnorm = sqrt(box->trialqx[i]*box->trialqx[i] + box->trialqy[i]*box->trialqy[i] + box->trialqz[i]*box->trialqz[i] + box->trialqw[i]*box->trialqw[i]);
		
		box->trialqx[i] /= quatnorm;
		box->trialqy[i] /= quatnorm;
		box->trialqz[i] /= quatnorm;
		box->trialqw[i] /= quatnorm;
		
        // generate trial positions of patches by rotation with trial quaternion (remember to add trialCOM and radius)
        
        // first we need the angle and the axes -- these will be overwritten every trial
// 		printf("old angle %2.12f and axis %2.12f %2.12f %2.12f for the old quat %f %f %f %f\n",box->cosRotAngle[i],box->axisx[i],box->axisy[i],box->axisz[i],box->oldqx[i],box->oldqy[i],box->oldqz[i],box->oldqw[i]);
		
        box->cosRotAngle[i] = findAxesQuat(box->trialqx[i],box->trialqy[i],box->trialqz[i],box->trialqw[i],spamaxis);
        
        box->axisx[i] = spamaxis[0];
        box->axisy[i] = spamaxis[1];
        box->axisz[i] = spamaxis[2];
//         printf("new angle %2.12f and axis %2.12f %2.12f %2.12f for the trial quat %f %f %f %f\n",box->cosRotAngle[i],box->axisx[i],box->axisy[i],box->axisz[i],box->trialqx[i],box->trialqy[i],box->trialqz[i],box->trialqw[i]);
        // once trial axis and angle are identified, find each patch pos and rotate to generate trialpatchpos
        for(j=1;j<=box->patchList[ii][0];j++)
        {
			jj = box->patchList[ii][j];
// 			printf("why looping %d %d %d %d\n",ii,i,j,box->patchList[ii][0]);
            rotate(box->patchposx[ii][j],box->patchposy[ii][j],box->patchposz[ii][j],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],rotpatch);
//             printf("patch %d on chrom %d fund pos %f %f %f -- rotated pos %f %f %f\n",j,ii,box->patchposx[ii][j],box->patchposy[ii][j],box->patchposz[ii][j],rotpatch[0],rotpatch[1],rotpatch[2]);
            // now assign the global positions of the patch positions -- COM + rotpatch*radius
			// find rotpatch norm and renormalise to one
			rotpatchnorm=sqrt(rotpatch[0]*rotpatch[0] + rotpatch[1]*rotpatch[1] + rotpatch[2]*rotpatch[2]);
// 			printf("norm of rotated patch %f\n",rotpatchnorm);
			box->trialGpatchposx[ii][j] = box->trialX[i] + box->radius[i]*rotpatch[0]/rotpatchnorm; // position of patches
			box->trialGpatchposy[ii][j] = box->trialY[i] + box->radius[i]*rotpatch[1]/rotpatchnorm;
			box->trialGpatchposz[ii][j] = box->trialZ[i] + box->radius[i]*rotpatch[2]/rotpatchnorm;
			
// 			printf("final rotated translated trial pos for patch %d on c %d %f %f %f\n",j,ii,box->trialGpatchposx[ii][j],box->trialGpatchposy[ii][j],box->trialGpatchposz[ii][j]);
// 			printf("final trial pos for fund bead %d on c %d -- %d %f %f %f\n\n",j,ii,jj,box->X[jj],box->Y[jj],box->Z[jj]);
        }
//         printf("indices chrom %d ato %d patch %d patch fund %d -- %d\n",ii,i,j,box->patchList[ii][j],box->patchList[ii][0]);
    }
    if(chainflag==1)
        printf("success rot with new dt %d %f\n",chainflag,p->dt);
    // generate trial positions of all beads, COM chromosomes and Gpatchpositions
    
}

int checkStretch(struct SYSTEM *box, struct PARAINPUT *p)
{
    
    // subroutine to check stretching for the new positions
    // we use updated trial patch positions obtained after translation and rotation steps
    // cycle through 
    
    
	int chainflag,ch,i,n1,ii,j,jj;
	chainflag=0;
	double diffX,diffY,diffZ,rij,rijsq,r0;
    
    for(ii=0;ii<p->nChrom;ii++)
    {
        i = box->chromList[ii];
// 		printf("checking sep chrom %d par %d has %d patches\n",ii,i,box->patchList[ii][0]);
        for(jj=1;jj<=box->patchList[ii][0];jj++)
        {
			
            j = box->patchList[ii][jj]; // which bead
            // distance of bead from the respective global patch position
            diffX = box->trialX[j] - box->trialGpatchposx[ii][jj];
            diffY = box->trialY[j] - box->trialGpatchposy[ii][jj];
            diffZ = box->trialZ[j] - box->trialGpatchposz[ii][jj];
            
            rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
            if(rijsq > 0.0)
            rij=sqrt(rijsq);
            else
            rij = 0.0;
			
            r0 = 3.0*(box->radius[i] + box->radius[j]);

            if(rij >= r0)
            {
                printf("%d %d -- %d %d (points %d) sep %2.12f %2.12f -- restoring?? %f\n",ii,i,jj,j,box->patchPointer[ii][jj],rij,r0,box->fx[j]*diffX + box->fy[j]*diffY + box->fz[j]*diffZ);
                printf("%2.12f %2.12f %2.12f -- %2.12f %2.12f %2.12f \n",box->fRx[j],box->fRy[j],box->fRz[j],box->activity[j]*box->actx[j],box->activity[j]*box->acty[j],box->activity[j]*box->actz[j]);

                chainflag=1; // extension beyond limits
                
            }
            
        }
    }
	// check wall -- NOTE: exact condition unclear, whether overlap beyond wall or centre beyond wall.

	
	return chainflag;
}

void integrateLangevinEM(struct SYSTEM *box,struct PARAINPUT *p)
{
    
    int i,j,ii,lcount;
    double displacesq;
    double rgx,rgy,rgz;
    double rgori;
	
	double diffX,diffY,diffZ,rijsq,rij,r0;
	int chainflag,n1,ch,oldchain;
    box->maxdispsq=0.0;

	p->dt = p->origdt; // reset timestep
	for(i=0;i<p->nAto;i++)
	{
	
		// noise direction for this step
		box->printflag[i] = box->type[i];
		box->rgx[i] = gauss_var_AT(0,1);
		box->rgy[i] = gauss_var_AT(0,1);
		box->rgz[i] = gauss_var_AT(0,1);

        // rotational noise for active forces

        box->argx[i] = gauss_var_AT(0,1);
		box->argy[i] = gauss_var_AT(0,1);
		box->argz[i] = gauss_var_AT(0,1);
		
        // NOTE: later set rotational noise for torque equation
		
	}
	oldchain=0;
	updateTrial(&(*box),&(*p),0); //generate first trial position
	chainflag=checkStretch(&(*box),&(*p)); // check if all is gut
// 	chainflag=0;
	while(chainflag==1) // keep redoing the timestep till the chainflag thing comes zero
	{
        printf("stretched bond!! -- %2.12f\n",p->dt);
		p->dt = p->dt*0.01;
		updateTrial(&(*box),&(*p),1); // get trial positions for current timestep value
        oldchain = chainflag;
		chainflag=checkStretch(&(*box),&(*p));
		printf("resolved?? %d\n",chainflag);
        if(p->dt < 0.0000000000000001 && chainflag ==1) // 0.0000000000000001
        {
            printf("exit\n");
		exit(0);
        }
	}
	// if no bond stretched beyond limit	
// 	printf("chain %d\n",chainflag);
// 	chainflag=0;
	if(chainflag==0)
	{
        // update real positions from trial positions -- all chromosomes and patches
		for(i=0;i<p->nAto;i++)
		{
			box->X[i] = box->trialX[i];
			box->Y[i] = box->trialY[i];
			box->Z[i] = box->trialZ[i];
			
			box->vx[i] = (box->X[i] - box->oldX[i])/p->dt;
			box->vy[i] = (box->Y[i] - box->oldY[i])/p->dt;
			box->vz[i] = (box->Z[i] - box->oldZ[i])/p->dt;
			
			box->oldX[i] = box->X[i];
			box->oldY[i] = box->Y[i];
			box->oldZ[i] = box->Z[i];

            box->oldactx[i] = box->actx[i];
            box->oldacty[i] = box->acty[i];
            box->oldactz[i] = box->actz[i];
		}

		// update rotated patch positions on each chromosome 
		for(ii=0;ii<p->nChrom;ii++)
		{
			i=box->chromList[ii];
// 		    printf("updating quats old %f %f %f %f -- trial %f %f %f %f\n",box->oldqx[i],box->oldqy[i],box->oldqz[i],box->oldqw[i],box->trialqx[i],box->trialqy[i],box->trialqz[i],box->trialqw[i]);
			box->qx[i] = box->trialqx[i]; // NOTE: redundant to put here but ensures we do not updated quaternions for patches
			box->qy[i] = box->trialqy[i];
			box->qz[i] = box->trialqz[i];
			box->qw[i] = box->trialqw[i];
			
			box->oldqx[i] = box->qx[i];
			box->oldqy[i] = box->qy[i];
			box->oldqz[i] = box->qz[i];
			box->oldqw[i] = box->qw[i];
// 			printf("updated quats now %f %f %f %f -- new old %f %f %f %f\n",box->qx[i],box->qy[i],box->qz[i],box->qw[i],box->oldqx[i],box->oldqy[i],box->oldqz[i],box->oldqw[i]);
			for(j=1;j<=box->patchList[ii][0];j++)
			{
				box->Gpatchposx[ii][j] = box->trialGpatchposx[ii][j];
				box->Gpatchposy[ii][j] = box->trialGpatchposy[ii][j];
				box->Gpatchposz[ii][j] = box->trialGpatchposz[ii][j];
// 				printf("final pos patch %d chrom %d %f %f %f\n",j,ii,box->Gpatchposx[ii][j],box->Gpatchposy[ii][j],box->Gpatchposz[ii][j]);
// 				printf("corr atom %d and pos %f %f %f\n",box->patchList[ii][j],box->X[box->patchList[ii][j]],box->Y[box->patchList[ii][j]],box->Z[box->patchList[ii][j]]);
			}
			
		}
		if(oldchain != chainflag)
            printf("successful update of all things\n");
    }
        
}


// void updateSurfTrial(struct SYSTEM *box,struct PARAINPUT *p,int chainflag)
// {
//
// }
//
// int checkSurfStretch(struct SYSTEM *box,struct PARAINPUT *p)
// {
//
//
//
// }


void integrateLangevinEMSurf(struct SYSTEM *box,struct PARAINPUT *p)
{
    // Subroutine to update the positions of the surface particles
    // we do Euler-Maruyama update of positions based on their net forces
    // generate the Brownian kicks and update position here

    double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5);
    double halfsep = 0.5*avesep;
//     double surfviscoterm = 6.0*M_PI*p->surftauinv*p->gamma*halfsep/p->surfDiffu; // 6 pi eta r x surf_visco_scaler > 1 (surf visco scaler slows down surface dynamics compared to chrom and patch dynamics
    double surfviscoterm = 6.0*M_PI*p->surftauinv*p->gamma*halfsep;

    int i,j,chainflag;
    double rgx,rgy,rgz;


//     // set surf rgx,rgy,rgz
//     for(i=0;i<p->nSurfPoints;i++)
//     {
//         box->surfrgx = gauss_var_AT(0,1);
//         box->surfrgy = gauss_var_AT(0,1);
//         box->surfrgz = gauss_var_AT(0,1);
//     }
//
//     //generate trial and check stretch
//     updateSurfTrial(&(*box),&(*p),0);
//     chainflag=checkSurfStretch(&(*box),&(*p));

    p->dt = p->origdt;
    for(i=0;i<p->nSurfPoints;i++)
    {
        // harmonic force with fundamental point
        // Brownian forces -- using surfDiffu
        rgx = gauss_var_AT(0,1);
        rgy = gauss_var_AT(0,1);
        rgz = gauss_var_AT(0,1);
//         printf("surf forces %d - %f %f %f -- random %f %f %f\n",i,box->surffx[i],box->surffy[i],box->surffz[i],sqrt(2.0*p->surfDiffu/p->dt)*rgx,sqrt(2.0*p->surfDiffu/p->dt)*rgy,sqrt(2.0*p->surfDiffu/p->dt)*rgz);
        // update is old + spring force + brownian noise
//         box->surfPointx[i] = box->surfPointx[i] + p->dt*(box->surffx[i]*p->surfDiffu + sqrt(2.0*p->surfDiffu/p->dt)*rgx);
//         box->surfPointy[i] = box->surfPointy[i] + p->dt*(box->surffy[i]*p->surfDiffu + sqrt(2.0*p->surfDiffu/p->dt)*rgy);
//         box->surfPointz[i] = box->surfPointz[i] + p->dt*(box->surffz[i]*p->surfDiffu + sqrt(2.0*p->surfDiffu/p->dt)*rgz);

        // NOTE: Re-scale point motion to be lower than geometry motion

//         box->surfPointx[i] = box->surfPointx[i] + p->dt*(box->surffx[i]*p->surfDiffu/p->surftauinv + sqrt(2.0*p->surfDiffu/(p->dt*p->surftauinv))*rgx);
//         box->surfPointy[i] = box->surfPointy[i] + p->dt*(box->surffy[i]*p->surfDiffu/p->surftauinv + sqrt(2.0*p->surfDiffu/(p->dt*p->surftauinv))*rgy);
//         box->surfPointz[i] = box->surfPointz[i] + p->dt*(box->surffz[i]*p->surfDiffu/p->surftauinv + sqrt(2.0*p->surfDiffu/(p->dt*p->surftauinv))*rgz);

        // NOTE: dX/dt = fx*D/k_BT + sqrt(2D) G(0,1) = fx*(1/surfviscoterm) + sqrt(2 k_B T / surfviscoterm) G(0,1)

        box->surfPointx[i] = box->surfPointx[i] + p->dt*(box->surffx[i]/surfviscoterm + sqrt(2.0*p->tem/(surfviscoterm*p->dt))*rgx);
        box->surfPointy[i] = box->surfPointy[i] + p->dt*(box->surffy[i]/surfviscoterm + sqrt(2.0*p->tem/(surfviscoterm*p->dt))*rgy);
        box->surfPointz[i] = box->surfPointz[i] + p->dt*(box->surffz[i]/surfviscoterm + sqrt(2.0*p->tem/(surfviscoterm*p->dt))*rgz);
	}

}

void updateGeom(struct SYSTEM *box,struct PARAINPUT *p,double currsig0x, double currsig0z)
{
    // NOTE: subroutine to update the geometry based on external stress, delta stress, membrane elasticity (surfSpringSelf), timescale (surfDiff*nSurfPoints inverse)

    double velx,vely,velz;
    double oldpx,oldpy,oldpz;

    oldpx = p->ella;
    oldpy = p->ellb;
    oldpz = p->ellc;


    // calculate velocity
//     velx = p->sig0x/p->surfEta + (p->sigdel/p->surfEta)*rgx - (p->surfSpringSelf/p->surfEta)*p->ella;
//     vely = p->sig0y/p->surfEta + (p->sigdel/p->surfEta)*rgy - (p->surfSpringSelf/p->surfEta)*p->ellb;

    velx = currsig0x*p->surftauinv/p->surfEta - (p->surfSpringSelf/p->surfEta)*p->ella;
//     vely = currsig0y/p->surfEta - (p->surfSpringSelf/p->surfEta)*p->ellb;
//     velz = p->sig0z/p->surfEta + p->sigdel/p->surfEta*rgz - (p->surfSpringSelf/p->surfEta)*p->ellc;
    velz = currsig0z*p->surftauinv/p->surfEta - (p->surfSpringSelf/p->surfEta)*p->ellc;

//     vely = p->sig0y/p->surfEta + p->sigdel/p->surfEta*rgy - (p->surfSpringSelf/p->surfEta)*p->ellb;


    // update extensions based on velocity using finite difference update

    p->ella = p->ella + p->dt*velx;
//     p->ellb = p->ellb + p->dt*vely;
    p->ellc = p->ellc + p->dt*velz;
    p->ellb = oldpx*oldpy*oldpz/(p->ella*p->ellc);

//     printf("current geom %f %f %f -- vel %f %f %f -- comp %f %f %f\n",p->ella,p->ellb,p->ellc,velx,vely,velz,p->sig0x,p->sigdel,p->surfEta);
//     genSurf(p,box);

    //MCRelaxSurf(p,box,2000);

    // rescale the surfFund positions for the new ellipsoid geometry
    scaleSurf(p,box);



}
