void placeChromMCInit(struct PARAINPUT *p,struct SYSTEM *box)
{
// NOTE: Subroutine to place chromosomes iteratively in a MC protocol
	// first we place all chromosome centroids around (0,0,0)
	// then we enter a loop for p->nMCIter iterations -- each iteration calculate total energy and accept/reject based on
	int i,j,k,ii,jj;
	double curr[3];
	double rX,rY,rZ,r;
	double norm1,norm2,spam;
	int atomdone, atomnotdone;
	double random,distCent,randx,randy,randz;
	double en=0.0,prob=0.0;
	char buffer[600],buffer2[600];
	double mspam = pow(p->massscalefact/8.0,0.3333);

	printf("entering MC iteration -- avoiding big overlaps\n");
	for(ii=0;ii<p->nDChrom;ii++)
	{


		i=box->chromList[ii]; // the global index
// 		printf("%d %d\n",i,ii);


		atomdone=0;
		while(atomdone==0)
		{
			random=drand48();
			randx=drand48();
			randy=drand48();
			randz=drand48();

			if(random<0.5)
				curr[0] = -2.2*p->constraint_freq*box->chromDist[ii]*p->ella*randx*(0.5+box->chromType[ii]/(2.0*p->nChrom));
			else
				curr[0] = 2.2*p->constraint_freq*box->chromDist[ii]*p->ella*randx*(0.5+box->chromType[ii]/(2.0*p->nChrom));

			random=drand48();

			if(random<0.5)
				curr[1] = -1.7*p->constraint_freq*box->chromDist[ii]*p->ellb*randy*(0.5+box->chromType[ii]/(2.0*p->nChrom));
			else
				curr[1] = 1.7*p->constraint_freq*box->chromDist[ii]*p->ellb*randy*(0.5+box->chromType[ii]/(2.0*p->nChrom));


			random=drand48();

			if(random<0.5)
				curr[2] = -1.2*p->constraint_freq*box->chromDist[ii]*p->ellc*randz*(0.5+box->chromType[ii]/(2.0*p->nChrom));
			else
				curr[2] = 1.2*p->constraint_freq*box->chromDist[ii]*p->ellc*randz*(0.5+box->chromType[ii]/(2.0*p->nChrom));

// 			printf("%d %d - %f %f %f - %f %f %f %f\n",ii,i,curr[0],curr[1],curr[2],box->radius[i],p->ella,p->ellb,p->ellc);
			atomnotdone=0;
			for(jj=0;jj<ii;jj++)
			{

				j = box->chromList[jj];
				// IPD energy -- find interparticle distance and desired IPD
				rX = box->X[j] - curr[0];
				rY = box->Y[j] - curr[1];
				rZ = box->Z[j] - curr[2];

				r = sqrt(rX*rX + rY*rY + rZ*rZ);
// 				printf("%d %d %f %f %f %f %f\n",ii,jj,r,curr[0],curr[1],curr[2],box->radius[i]+box->radius[j]);
				if(box->chromType[ii] == box->chromType[jj]) // more stringent for homologues
				{
					if(r <= 0.2*(mspam*box->radius[i]+mspam*box->radius[j])) // overlapping
					{
						atomnotdone=1;

					}
				}
				else
				{

					if(r <= 0.25*(mspam*box->radius[i]+mspam*box->radius[j])) // overlapping
					{
						atomnotdone=1;
					}
				}

			}
			if(fabs(curr[0])>p->ella - box->radius[i] || fabs(curr[1])>p->ellb - box->radius[i] || fabs(curr[2]) > p->ellc - box->radius[i])
			{
				atomnotdone=1;
				printf("%d - x- %f %f  -y- %f %f -z- %f %f -r %f\n",ii,curr[0],p->ella,curr[1],p->ellb,curr[2],p->ellc,box->radius[i]);
			}

			atomdone = 1 - atomnotdone;
			box->X[i] = curr[0];
			box->Y[i] = curr[1];
			box->Z[i] = curr[2];
			box->Xpre[i] = curr[0];
			box->Ypre[i] = curr[1];
			box->Zpre[i] = curr[2];
// 			printf("%d %f %f %f\n",ii,curr[0],curr[1],curr[2]);
		}


	}

}

void ChromMCStepConstraint(struct PARAINPUT *p,struct SYSTEM *box)
{
    int i,j,k,ii,jj;
	double curr[3];
	double rX,rY,rZ,r;
	double norm1,norm2,spam;
	int atomdone, atomnotdone;
	double random,distCent;
	double en=0.0,prob=0.0,deltaen=0.0,olden=0.0;
	char buffer[600],buffer2[600];
// 	printf("entering MC iteration\n");
	// now that initial locations are fixed, calculate energy once
// 	printf("calling first energy calc\n");

	int trial;
	double oldX,oldY,oldZ;


    // can later incorporate volume change etc in the MC step process

// 	printf("mc stepping\n");

    // do p->nDChrom trial displacements

    for(ii=0;ii<p->nDChrom;ii++)
    {
        // pick a chrom at random to displace
        random=drand48();
        trial = (int) (p->nDChrom*drand48());
        if(trial >= p->nDChrom)
            trial = p->nDChrom-1;
        if(trial < 0)
            trial = 0;

        i = box->chromList[ii];
        oldX = box->X[i];
        oldY = box->Y[i];
        oldZ = box->Z[i];

		// en pre move
// 		printf("chrom %d storing\n",i);
// 		olden = chromEnergy(p,box,i);

		olden = chromIPDBias(p,box,ii); // old IPD bias energy for this chromosome type
		olden += surfEnergy(p,box,i);

		// displace a chromosome
        random = 2.0*drand48()-1.0;

        box->X[i] = box->X[i] + p->mcdisp*random;

        random = 2.0*drand48()-1.0;

        box->Y[i] = box->Y[i] + p->mcdisp*random;

        random = 2.0*drand48()-1.0;

        box->Z[i] = box->Z[i] + p->mcdisp*random;

        // calculate energy with these new positions
// 			printf("chrom %d, id %d, type %d displaced to %f %f %f\n",ii,i,box->chromType[ii],box->X[i],box->Y[i],box->Z[i]);
//         en = chromEnergy(p,box,i);

		// recalc all biasvalues for this displacement
		en = chromIPDBias(p,box,ii);
		en += surfEnergy(p,box,i);

		deltaen = en - olden;

        prob = exp(-deltaen/p->mctemp);
// 		printf("deltae prob %f %f %f\n",deltaen,prob,box->curren);
// 		printf("old %f %f %f new %f %f %f\n",oldX,oldY,oldZ,box->X[i],box->Y[i],box->Z[i]);

        random=drand48();
        if(random >= prob) // reject
        {
            box->X[i] = oldX;
            box->Y[i] = oldY;
            box->Z[i] = oldZ;
            box->currW = box->currW;
        }
        else // do nothing
        {
            box->currW = box->currW + deltaen;
        }
//         printf("old en %f -- new en %f -- den %f -- fact %f -- current %f\n",olden,en,deltaen,prob,box->curren);
//         printf("trial %d chrom %d par %d -- max %d\n",ii,trial,i,p->nAto-1);

    }

// 	printf("exiting mc step\n");


}

void ChromMCStepLinfConstraint(struct PARAINPUT *p,struct SYSTEM *box)
{
	int i,j,k,ii,jj;
	double curr[3];
	double rX,rY,rZ,r;
	double norm1,norm2,spam;
	int atomdone, atomnotdone;
	double random,distCent;
	double en=0.0,prob=0.0,deltaen=0.0,olden=0.0,spamen=0.0;
	char buffer[600],buffer2[600];
// 	printf("entering MC iteration\n");
	// now that initial locations are fixed, calculate energy once
// 	printf("calling first energy calc\n");

	int trial;
	double oldX,oldY,oldZ;

	for(ii=0;ii<p->nDChrom;ii++)
    {
        // pick a chrom at random to displace
        random=drand48();
        trial = (int) (p->nDChrom*drand48());
        if(trial >= p->nDChrom)
            trial = p->nDChrom-1;
        if(trial < 0)
            trial = 0;

        i = box->chromList[trial];

// 		printf("trial chrom %d global %d\n",trial,i);
        oldX = box->X[i];
        oldY = box->Y[i];
        oldZ = box->Z[i];

		// calculate deviations and identify the L_infty
		spamen = chromIPDBias(p,box,trial); // old IPD bias energy for this chromosome type
		spamen = surfEnergy(p,box,i);

		calcDeviations(p,box,0);

		olden = p->k_ipd*box->l_infty + spamen;
// 		printf("chrom %d (g %d) max dev %f -- surfen %f -- olden %f\n",ii,i,box->l_infty,spamen,olden);

		// displace a chromosome
        random = 2.0*drand48()-1.0;
        box->X[i] = box->X[i] + p->mcdisp*random;

        random = 2.0*drand48()-1.0;
        box->Y[i] = box->Y[i] + p->mcdisp*random;

        random = 2.0*drand48()-1.0;
        box->Z[i] = box->Z[i] + p->mcdisp*random;


		// recalculate the constraint deviations and identify the

		spamen = chromIPDBias(p,box,trial);
		spamen = surfEnergy(p,box,i);

		calcDeviations(p,box,0);

		en = p->k_ipd*box->l_infty + spamen;
// 		printf("chrom %d (g %d) max dev %f -- surfen %f -- newen %f\n",ii,i,box->l_infty,spamen,en);

		deltaen = en - olden;

        prob = exp(-deltaen/p->mctemp);
// 		printf("deltae prob %f %f %f\n",deltaen,prob,box->curren);
// 		printf("old %f %f %f new %f %f %f\n",oldX,oldY,oldZ,box->X[i],box->Y[i],box->Z[i]);

        random=drand48();
        if(random >= prob) // reject
        {
            box->X[i] = oldX;
            box->Y[i] = oldY;
            box->Z[i] = oldZ;
				spamen = chromIPDBias(p,box,trial);
				spamen = surfEnergy(p,box,i);

			calcDeviations(p,box,0);

            box->currW = box->currW;
        }
        else // do nothing
        {
            box->currW = box->currW + deltaen;
        }
	}
// 	printf("did an MC sweep\n");

}


void ChromMCStepElast(struct PARAINPUT *p,struct SYSTEM *box)
{
    int i,j,k,ii,jj;
	double curr[3];
	double rX,rY,rZ,r;
	double norm1,norm2,spam;
	int atomdone, atomnotdone;
	double random,distCent;
	double en=0.0,prob=0.0,deltaen=0.0,olden=0.0;
	char buffer[600],buffer2[600];
// 	printf("entering MC iteration\n");
	// now that initial locations are fixed, calculate energy once
// 	printf("calling first energy calc\n");

	int trial;
	double oldX,oldY,oldZ;


    // can later incorporate volume change etc in the MC step process

// 	printf("mc stepping\n");

    // do p->nDChrom trial displacements

    for(ii=0;ii<p->nDChrom;ii++)
    {
        // pick a chrom at random to displace
        random=drand48();
        trial = (int) (p->nDChrom*drand48());
        if(trial >= p->nDChrom)
            trial = p->nDChrom-1;
        if(trial < 0)
            trial = 0;

        i = box->chromList[trial];
        oldX = box->X[i];
        oldY = box->Y[i];
        oldZ = box->Z[i];

		// en pre move
// 		printf("chrom %d storing\n",i);
		olden = chromEnergy(p,box,i);

// 		olden = chromIPDBias(p,box,ii); // old IPD bias energy for this chromosome type

		// displace a chromosome
        random = 2.0*drand48()-1.0;

        box->X[i] = box->X[i] + p->mcdisp*random;

        random = 2.0*drand48()-1.0;

        box->Y[i] = box->Y[i] + p->mcdisp*random;

        random = 2.0*drand48()-1.0;

        box->Z[i] = box->Z[i] + p->mcdisp*random;

        // calculate energy with these new positions
// 			printf("chrom %d, id %d, type %d displaced to %f %f %f\n",ii,i,box->chromType[ii],box->X[i],box->Y[i],box->Z[i]);
        en = chromEnergy(p,box,i);

		// recalc all biasvalues for this displacement
// 		en = chromIPDBias(p,box,ii);

		deltaen = en - olden;

        prob = exp(-deltaen/p->mctemp);
// 		printf("deltae prob %f %f %f\n",deltaen,prob,box->curren);
// 		printf("old %f %f %f new %f %f %f\n",oldX,oldY,oldZ,box->X[i],box->Y[i],box->Z[i]);

        random=drand48();
        if(random >= prob) // reject
        {
            box->X[i] = oldX;
            box->Y[i] = oldY;
            box->Z[i] = oldZ;
            box->curren = box->curren;
        }
        else // do nothing
        {
            box->curren = box->curren + deltaen;
        }
//         printf("old en %f -- new en %f -- den %f -- fact %f -- current %f\n",olden,en,deltaen,prob,box->curren);
//         printf("trial %d chrom %d par %d -- max %d\n",ii,trial,i,p->nAto-1);

    }

// 	printf("exiting mc step\n");


}

void ChromGradStep(struct PARAINPUT *p,struct SYSTEM *box,int max)
{
	/// update all positions for indices corresponding to chromosomes, so loop only over chromosome indices

	// note that since the IPD/constraint entry is an unweighted average over 4 distances, we need the gradient transferred to the chromosomes themselves
	int ii,i;
	double curr[3];
	double rX,rY,rZ,r;
	double norm1,norm2,spam;
	double en=0.0,prob=0.0,deltaen=0.0,olden=0.0,spamen=0.0;

	int trial;
	double oldX,oldY,oldZ;
	double normf;

	/// calculate gradients for this step
// 	printf("calculating gradients -- %d\n",max);
	calcAllGradients(p,box,max);

	// we now have all the gradients -- if max==1 then only the forces for maxi,maxj related chrom are non-zero and only they will move
	for(ii=0;ii<p->nDChrom;ii++)
	{
		// move each chromosome by gradient x mcdisp
		i = box->chromList[ii];
		trial=ii;
// 		printf("%d %d %d\n",ii,i,trial);
        oldX = box->X[i];
        oldY = box->Y[i];
        oldZ = box->Z[i];

		spamen = chromIPDBias(p,box,trial); // old IPD bias energy for this chromosome type
		spamen = surfEnergy(p,box,i);

		calcDeviations(p,box,0);

		olden = p->k_ipd*box->l_infty + spamen;
		normf = box->fx[i]*box->fx[i] + box->fy[i]*box->fy[i] + box->fz[i]*box->fz[i];

		//if(normf > 0.0)
		//{
		//	normf = sqrt(normf);
		//	normf = MIN(1.0,normf);
		//	box->X[i] = box->X[i] + p->mcdisp*box->fx[i]/normf;
		//	box->Y[i] = box->Y[i] + p->mcdisp*box->fy[i]/normf;
        	//	box->Z[i] = box->Z[i] + p->mcdisp*box->fz[i]/normf;
		//}
		box->X[i] = box->X[i] + p->mcdisp*box->fx[i];
                box->Y[i] = box->Y[i] + p->mcdisp*box->fy[i];
		box->Z[i] = box->Z[i] + p->mcdisp*box->fz[i];


		// recalculate the constraint deviations and identify the

		spamen = chromIPDBias(p,box,trial);
		spamen = surfEnergy(p,box,i);

		calcDeviations(p,box,0);

		en = p->k_ipd*box->l_infty + spamen;
// 		printf("chrom %d (g %d) max dev %f -- surfen %f -- newen %f\n",ii,i,box->l_infty,spamen,en);

		deltaen = en - olden;

	}
// 	printf("updated positions -- %d -- %f\n",max,box->l_infty);


}

void ChromSwapStep(struct PARAINPUT *p,struct SYSTEM *box)
{
	int ii,i,jj,j;
	double curr[3];
	double rX,rY,rZ,r,random;
	double norm1,norm2,spam;
	double en=0.0,prob=0.0,deltaen=0.0,olden=0.0,newen=0.0,spamen=0.0;

	int trial1,trial2;
	double oldX,oldY,oldZ;


	random=drand48();
	trial1 = (int) (p->nDChrom*drand48());

	if(trial1 >= p->nDChrom)
		trial1 = p->nDChrom-1;
	if(trial1 < 0)
		trial1 = 0;

	i = box->chromList[trial1];
	oldX = box->X[i];
	oldY = box->Y[i];
	oldZ = box->Z[i];

	// en pre move
// 		printf("chrom %d storing\n",i);
// 	olden = chromEnergy(p,box,i);


	/// next candidate
	random=drand48();
	trial2 = (int) (p->nDChrom*drand48());
	while(box->chromType[trial2] == box->chromType[trial1])
	{
		random=drand48();
		trial2 = (int) (p->nDChrom*drand48());
	}

	if(trial2 >= p->nDChrom)
		trial2 = p->nDChrom-1;
	if(trial2 < 0)
		trial2 = 0;

	j = box->chromList[trial2];

// 	olden += chromEnergy(p,box,j);
// 	spamen = chromIPDBias(p,box,trial1); // old IPD bias energy for this chromosome type
//
// 	spamen = chromIPDBias(p,box,trial2); // old IPD bias energy for this chromosome type
//
//
// 	spamen = surfEnergy(p,box,i);
// 	spamen += surfEnergy(p,box,j);
//
// 	calcDeviations(p,box,0);
	calcBias(p,box);
	printf("%d %d old %f --  maxdev %d, ij %d %d,devdir %f, ex %f\n",box->chromType[trial1],box->chromType[trial2],box->l_infty,
		  box->maxdev,box->maxi,box->maxj,box->devdir,box->expectedMat[box->maxi][box->maxj]);
	olden = p->k_ipd*box->l_infty;// + spamen;

	// execute swap
	box->X[i] = box->X[j];
	box->Y[i] = box->Y[j];
	box->Z[i] = box->Z[j];

	box->X[j] = oldX;
	box->Y[j] = oldY;
	box->Z[j] = oldZ;

	// check new deviations

// 	spamen = chromIPDBias(p,box,trial1); // old IPD bias energy for this chromosome type
//
// 	spamen = chromIPDBias(p,box,trial2); // old IPD bias energy for this chromosome type
//
//
// 	spamen = surfEnergy(p,box,i);
// 	spamen += surfEnergy(p,box,j);
//
//
// 	calcDeviations(p,box,0);
	calcBias(p,box);
	printf("%d %d new %f --  maxdev %d, ij %d %d,devdir %f, ex %f\n",box->chromType[trial1],box->chromType[trial2],box->l_infty,
		  box->maxdev,box->maxi,box->maxj,box->devdir,box->expectedMat[box->maxi][box->maxj]);
	newen = p->k_ipd*box->l_infty;// + spamen;

	if(newen <= olden) // accept
	{

// 		spamen = chromIPDBias(p,box,trial1); // old IPD bias energy for this chromosome type
//
// 		spamen = chromIPDBias(p,box,trial2); // old IPD bias energy for this chromosome type
//
//
// 		spamen = surfEnergy(p,box,i);
// 		spamen += surfEnergy(p,box,j);
//
// 		calcDeviations(p,box,0);
		calcBias(p,box);
		printf("swapping %d %d -- %f %f L_inf %f, L1 %f, maxdev %d, ij %d %d,devdir %f, ex %f n constraint %d, proxLinf %f\n",box->chromType[trial1],box->chromType[trial2],olden,newen,box->l_infty,box->l1,
		  box->maxdev,box->maxi,box->maxj,box->devdir,box->expectedMat[box->maxi][box->maxj],(int)((p->nChrom+1)*(p->nChrom+2)/2),box->proxLinf);
	}

	if(newen > olden) // reject -- i.e., revert
	{

		box->X[j] = box->X[i];
		box->Y[j] = box->X[i];
		box->Z[j] = box->X[i];

		box->X[i] = oldX;
		box->Y[i] = oldY;
		box->Z[i] = oldZ;

// 		spamen = chromIPDBias(p,box,trial1); // old IPD bias energy for this chromosome type
//
// 		spamen = chromIPDBias(p,box,trial2); // old IPD bias energy for this chromosome type
//
//
// 		spamen = surfEnergy(p,box,i);
// 		spamen += surfEnergy(p,box,j);
//
// 		calcDeviations(p,box,0);
		calcBias(p,box);


	}
}

