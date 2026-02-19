//overlap energy and boundary energy contribute to MC energy

double surfEnergy(struct PARAINPUT *p,struct SYSTEM *box,int i)
{
    int j;
    double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	double surfcut=0.0;

    double rX,rY,rZ,r,en=0.0;
//     printf("calc surf en for %d ch %d type %d at %f %f %f\n",i,box->whichChrom[i],box->chromType[box->whichChrom[i]],box->X[i],box->Y[i],box->Z[i]);
    double mspam = 1.0; //pow(p->massscalefact/8.0,0.3333);
    for(j=0;j<p->nSurfPoints;j++)
    {
        rX = box->X[i] - box->surfPointx[j];
        rY = box->Y[i] - box->surfPointy[j];
        rZ = box->Z[i] - box->surfPointz[j];

        r = sqr(rX) + sqr(rY) + sqr(rZ);

        // check if distance within WCA cut-off
//         surfcut = 0.5*avesep + mspam*box->radius[i];
        surfcut = 0.5*avesep + 0.9*box->radius[i];

//         if(r <= sqr(pow(2.0,1.0/6.0)*surfcut))
//         {
//             r = sqrt(r);
//             en += 4.0*( pow(surfcut/r,12.0) - pow(surfcut/r,6.0)) + 1.0;
//         }
        if(r <= surfcut*surfcut) // overlapping
        {
            r = sqrt(r);
            en += p->eps_hertz*pow((1-r/surfcut),2.5);
        }

    }
    if (en < 0.0)
        printf("DANGER DANGER neg surf en %f\n",en);
    return en;

}


double chromEnergy(struct PARAINPUT *p,struct SYSTEM *box,int i)
{
    // subroutine to calculate chromosome energy from purely overlap and boundary interactions (not OP interactions)
	int j,k,jj;
	double curr[3];
	double rX,rY,rZ,r,rij,patchcut;

	double en = 0.0;
	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	double surfcut=0.0;

//     i = box->chromList[ii];
//     printf("chrom %d calc\n",i);
    for(jj=0;jj<p->nDChrom;jj++)
    {
//         ii = box->whichChrom[i];
//         if(jj<=ii)
//             continue;
        j = box->chromList[jj];
        // IPD energy -- find interparticle distance and desired IPD
        rX = box->X[j] - box->X[i];
        rY = box->Y[j] - box->Y[i];
        rZ = box->Z[j] - box->Z[i];

        rij = sqrt(rX*rX + rY*rY + rZ*rZ);

        // find lamin patch indices and get their radii to add here
        patchcut = 1.4*(box->radius[i] + box->radius[j]);
        // apply Hertz force to i and j
        if(rij <= patchcut) // overlapping
        {
            en += p->eps_hertz*pow((1-rij/patchcut),2.5);
        }

    }
    // now calculate the boundary energy
    en += surfEnergy(p,box,i);

	return en;
}


double fullEnergy(struct PARAINPUT *p,struct SYSTEM *box)
{
    // subroutine to calculate chromosome energy from purely overlap and boundary interactions (not OP interactions)
	int i,j,k,jj,ii;
	double curr[3];
	double rX,rY,rZ,r,rij,patchcut;

	double en = 0.0;

//     double mspam = pow(p->massscalefact/8.0,0.3333);
    for(ii=0;ii<p->nDChrom;ii++)
    {

        i = box->chromList[ii];
        for(jj=0;jj<p->nDChrom;jj++)
        {
            if(jj<=ii)
                continue;
            j = box->chromList[jj];
            // IPD energy -- find interparticle distance and desired IPD
            rX = box->X[j] - box->X[i];
            rY = box->Y[j] - box->Y[i];
            rZ = box->Z[j] - box->Z[i];

            rij = sqrt(rX*rX + rY*rY + rZ*rZ);

            // find lamin patch indices and get their radii to add here

            patchcut = 1.4*(box->radius[i] + box->radius[j]);
            // apply Hertz force to i and j
            if(rij <= patchcut) // overlapping
            {
                en += p->eps_hertz*pow((1-rij/patchcut),2.5);
            }

        }
        // now calculate the boundary energy
        en +=surfEnergy(p,box,i);

    }

	return en;
}


double enDiffEps(struct PARAINPUT *p,struct SYSTEM *box,double epsrecv,double *Xrecv,double *Yrecv,double *Zrecv,double *surfXrecv,double *surfYrecv,double *surfZrecv)
{
    // subroutine to calculate chromosome energy from purely overlap and boundary interactions (not OP interactions)
	int i,j,k,jj,ii;
	double curr[3];
	double rX,rY,rZ,r,rij,patchcut;

	double en = 0.0;
    double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	double surfcut=0.0;
    double mspam = 1.0; //pow(p->massscalefact/8.0,0.3333);
    for(ii=0;ii<p->nDChrom;ii++)
    {

        i = box->chromList[ii];
        for(jj=0;jj<p->nDChrom;jj++)
        {
            if(jj<=ii)
                continue;
            j = box->chromList[jj];
            // IPD energy -- find interparticle distance and desired IPD
            rX = Xrecv[j] - Xrecv[i];
            rY = Yrecv[j] - Yrecv[i];
            rZ = Zrecv[j] - Zrecv[i];

            rij = sqrt(rX*rX + rY*rY + rZ*rZ);

            // find lamin patch indices and get their radii to add here

            patchcut = box->radius[i] + box->radius[j];
            // apply Hertz force to i and j
            if(rij <= patchcut) // overlapping
            {
                en += epsrecv*pow((1-rij/patchcut),2.5);
            }

        }
        // now calculate the boundary energy
        for(j=0;j<p->nSurfPoints;j++)
        {
            rX = Xrecv[i] - surfXrecv[j];
            rY = Yrecv[i] - surfYrecv[j];
            rZ = Zrecv[i] - surfZrecv[j];

            r = sqr(rX) + sqr(rY) + sqr(rZ);

            // check if distance within WCA cut-off
            surfcut = 0.5*avesep + mspam*box->radius[i];

            if(r <= sqr(pow(2.0,1.0/6.0)*surfcut))
            {
                r = sqrt(r);
                en += 4.0*( pow(surfcut/r,12.0) - pow(surfcut/r,6.0)) + 1.0;
            }

        }

    }

	return en;
}


