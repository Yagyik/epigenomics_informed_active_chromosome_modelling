double distCalc(double r1, double r2)
{
    double spam;

    spam = MIN(r1 - r2,0);

    return spam;
}


double monitorIPD(struct PARAINPUT *p,struct SYSTEM *box, int chi, int chj)
{
    // subroutine to calculate the deviation of the chi,chj pair from the provided IPD matrix entry
    // given the chromosome type number chi, chj -- we want 4 pairwise distances
    int ii,i,jj,j;
    double rX,rY,rZ,r;
    double avgdist = 0.0;
    int countmatches = 0;
    int typei = box->chromType[chi];
    int typej = box->chromType[chj];
    double spammax;
    for(ii = 0; ii < p->nDChrom;ii++)
    {
        i = box->chromList[ii];
        if(box->chromType[ii] == typei) // true twice
        {
            for(jj = 0; jj < p->nDChrom;jj++)
            {
                j = box->chromList[jj];

                if(box->chromType[jj] == typej) // true twice
                {
                    // we now know that i and j are the chrom positions
                    // we also know that chi and chj are not homologues
                    rX = box->X[j] - box->X[i];
                    rY = box->Y[j] - box->Y[i];
                    rZ = box->Z[j] - box->Z[i];

                    r = sqrt(rX*rX + rY*rY + rZ*rZ);
//                     avgdist += 0.5*r/p->ella;
//                     avgdist += sqrt(sqr(0.5*rX/p->ella) + sqr(0.5*rY/p->ellb) + sqr(0.5*rZ/p->ellc)); // raw distance between the chrom centres in units of the nucleus dimensions


//                     spammax = MAX(fabs(0.5*rX/p->ella),fabs(0.5*rY/p->ellb));
//                     spammax = MAX(spammax,fabs(0.5*rZ/p->ellc));
//                     spammax = 0.33*(fabs(0.5*rX/p->ella) + fabs(0.5*rY/p->ellb) + fabs(0.5*rZ/p->ellc));
                    spammax = sqrt(sqr(0.5*rX/p->ella) + sqr(0.5*rY/p->ellb) + sqr(0.5*rZ/p->ellc));
//                     spammax = sqrt((p->ella/p->ella)*sqr(0.5*rX/p->ella) + (p->ellb/p->ella)*sqr(0.5*rY/p->ellb) + (p->ellc/p->ella)*sqr(0.5*rZ/p->ellc));
                    avgdist += spammax;
//                     avgdist += sqrt(sqr(0.5*rX/p->ella) + sqr(0.5*rY/p->ellb)); // + sqr(0.5*rZ/p->ellc));
                    countmatches +=1;
                }
            }
        }
    }
    avgdist /= 4.0; // the average of the 4 distances
//     printf("ch %d (t %d) and ch %d (t %d) matched %d times with avg dis %f\n",chi,typei,chj,typej,countmatches,avgdist);
    if(box->ipddev[chi][chj] > 0.0)
    {
//         double track = (avgdist - box->ipd[chi][chj])/box->ipddev[chi][chj]; //; *sqrt(p->k_ipd*p->mctemp)/box->ipddev[chi][chj];
        double track = (avgdist - box->ipd[chi][chj])*sqrt(p->k_ipd*p->mctemp)/box->ipddev[chi][chj];
        if(fabs(track) >= box->maxavgIPDdist)
            box->maxavgIPDdist = fabs(track);
    }
    return avgdist;

}

double monitorHomologues(struct PARAINPUT *p,struct SYSTEM *box, int chi)
{
    int ii,i,jj,j;
    double rX,rY,rZ,r;
    double avgdist = 0.0;
    int typei = box->chromType[chi];
    double spammax;
    for(ii = 0; ii < p->nChrom;ii++)
    {
        i = box->chromList[ii];
        if(box->chromType[ii] == typei) // true once
        {
            for(jj = p->nChrom; jj < p->nDChrom;jj++)
            {
                j = box->chromList[jj];
                if(box->chromType[jj] == typei) // true once
                {
                    // we now know that i and j are the chrom positions
                    // we also know that i and j are both homologues (for chi)
                    rX = box->X[j] - box->X[i];
                    rY = box->Y[j] - box->Y[i];
                    rZ = box->Z[j] - box->Z[i];

                    r = sqrt(rX*rX + rY*rY + rZ*rZ);
//                     avgdist += 0.5*r/p->ella;
//                     avgdist += sqrt(sqr(0.5*rX/p->ella) + sqr(0.5*rY/p->ellb) + sqr(0.5*rZ/p->ellc)); // raw distance between the chrom centres in units of the nucleus dimensions


//                     spammax = MAX(fabs(0.5*rX/p->ella),fabs(0.5*rY/p->ellb));
//                     spammax = MAX(spammax,fabs(0.5*rZ/p->ellc));
//                     spammax = 0.33*(fabs(0.5*rX/p->ella) + fabs(0.5*rY/p->ellb) + fabs(0.5*rZ/p->ellc));
                    spammax = sqrt(sqr(0.5*rX/p->ella) + sqr(0.5*rY/p->ellb) + sqr(0.5*rZ/p->ellc));
                    avgdist += spammax;

//                     avgdist += sqrt(sqr(0.5*rX/p->ella) + sqr(0.5*rY/p->ellb)); // + sqr(0.5*rZ/p->ellc));
                }
            }
        }
    }
//     avgdist /= 1.0; // the average of the 4 distances
//     double track = (avgdist - box->ipdhomo[ii][1]); //*sqrt(p->k_homologue*p->mctemp);
    double track = (avgdist - box->ipdhomo[ii][1])*sqrt(p->k_homologue*p->mctemp);
    if(fabs(track) >= box->maxavghomologuedist)
        box->maxavghomologuedist = fabs(track);
    return avgdist;
}

double monitorCentroid(struct PARAINPUT *p,struct SYSTEM *box, int chi)
{
    // subroutine to calculate the deviation of the chromosome chi from its desired centroid position
    int ii,i;
    double r = 0.0;
    int typei = box->chromType[chi];
    double spammax;
    for(ii = 0; ii<p->nDChrom;ii++)
    {

        i = box->chromList[ii];
        if(box->chromType[ii] == typei) // true twice
        {
//             r += 0.5*sqrt(box->X[i]*box->X[i] + box->Y[i]*box->Y[i] + box->Z[i]*box->Z[i])/p->ella;
//             r += sqrt(sqr(0.5*box->X[i]/p->ella) + sqr(0.5*box->Y[i]/p->ellb) + sqr(0.5*box->Z[i]/p->ellc));


//             spammax = MAX(fabs(0.5*box->X[i]/p->ella),fabs(0.5*box->Y[i]/p->ellb));
//             spammax = MAX(spammax,fabs(0.5*box->Z[i]/p->ellc));
//             spammax = 0.33*(fabs(0.5*box->X[i]/p->ella) + fabs(0.5*box->Y[i]/p->ellb) + fabs(0.5*box->Z[i]/p->ellc));
            spammax = sqrt(sqr(0.5*box->X[i]/p->ella) + sqr(0.5*box->Y[i]/p->ellb) + sqr(0.5*box->Z[i]/p->ellc));
            r += spammax;
//             r += sqrt(sqr(0.5*box->X[i]/p->ella) + sqr(0.5*box->Y[i]/p->ellb)); // + sqr(0.5*box->Z[i]/p->ellc));
        }
    }
    r /= 2.0; // average over homologues
//     double track = (r - box->chromDist[chi])/box->chromDistDev[chi]; //*sqrt(p->k_centroid*p->mctemp)/box->chromDistDev[chi];
    double track = (r - box->chromDist[chi])*sqrt(p->k_centroid*p->mctemp)/box->chromDistDev[chi];
    if(fabs(track) >= box->maxavgcentroiddist)
        box->maxavgcentroiddist = fabs(track);
    return r;
}



double IPDBias(struct PARAINPUT *p,struct SYSTEM *box,int flag)
{
    // subroutine to calculate the energy cost of all the IPD deviations
    int i,j,k,ii,jj;
    double r = 0.0;
    double en = 0.0;
    double netspring = 0.0;

    for(ii=0;ii<p->nChrom;ii++)
	{

		for(jj=0;jj<p->nChrom;jj++)
		{
			if(jj<=ii)
				continue;

            // calculate the deviation of the average distance of this pair of homologues from each other

            if( box->ipdhomo[ii][0] != jj && box->ipdhomo[jj][0] != ii) // not homologue
			{
                netspring = 0.5*p->k_ipd/(box->ipddev[ii][jj]*box->ipddev[ii][jj]);
                r = monitorIPD(p,box,ii,jj); // distance is scaled to chromosome sizes
// 				printf("chrom %d and %d -- desired %f current %f -- spring const %f -- cost %f\n",ii,jj,box->ipd[ii][jj],r,netspring,netspring*(r - box->ipd[ii][jj])*(r - box->ipd[ii][jj]));
                if(flag==1)
                {
                    box->constraintMat[ii][jj] = r;
                    box->constraintMat[jj][ii] = r;

                }

				en += netspring*(r - box->ipd[ii][jj])*(r - box->ipd[ii][jj]);
                box->cumuldev += fabs(r - box->ipd[ii][jj]);
                box->cumuldenom += 1.0;

			}
        }

    }
    box->countConstraints +=1.0;

    return en;


}


double HomologueBias(struct PARAINPUT *p,struct SYSTEM *box,int flag)
{
    // subroutine to calculate the energy cost of all the IPD deviations
    int i,j,k,ii,jj;
    double r = 0.0;
    double en = 0.0;
    for(ii=0;ii<p->nChrom;ii++)
	{
        r = monitorHomologues(p,box,ii); // distance is scaled to chromosome sizes
        en += (0.5*p->k_homologue)*sqr(r - box->ipdhomo[ii][1]);
//         printf("homologues %d %d type %d -- desired %f current %f -- spring constant %f -- cost %f\n",ii,ii+p->nChrom,box->chromType[ii],box->ipdhomo[ii][1],r,0.5*p->k_homologue,0.5*p->k_homologue*sqr(r - box->ipdhomo[ii][1]));
        box->cumuldev += fabs(r - box->ipdhomo[ii][1]);
        box->cumuldenom += 1.0;

        if(flag==1)
        box->constraintMat[ii][ii] = r;
    }
    return en;

}


double CentroidBias(struct PARAINPUT *p,struct SYSTEM *box,int flag)
{
    // subroutine to calculate deviation of each chrom from the centroid
    int ii,i;
    double r,en=0.0;
    double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
    for(int ii=0;ii<p->nChrom;ii++)
    {

        r = monitorCentroid(p,box,ii); // distance is scale to nuclear dimensions -- averaged over the two homolouges of this chromosome
//         printf("chrom %d (r %f) s %f  is supposed to be %f away from centre but is %f -- spring constant %f -- cost %f\n",ii,box->radius[box->chromList[ii]],0.5*avesep,box->chromDist[ii],r,0.5*p->k_centroid/(box->chromDistDev[ii]*box->chromDistDev[ii]),0.5*p->k_centroid/(box->chromDistDev[ii]*box->chromDistDev[ii])*sqr(r - box->chromDist[ii]));
        en += 0.5*p->k_centroid/(box->chromDistDev[ii]*box->chromDistDev[ii])*sqr(r - box->chromDist[ii]);
        box->cumuldev += fabs(r - box->chromDist[ii]);
        box->cumuldenom += 1.0;
        if(flag==1)
        {
            box->constraintMat[ii][p->nChrom] = r;
            box->constraintMat[p->nChrom][ii] = r;
        }

    }
    return en;
}


double chromIPDBias(struct PARAINPUT *p,struct SYSTEM *box,int chi)
{
    int jj,typej;
    double netspring,r,en=0.0;
    // IPD contribution to energy

    // chi is the chromosome index -- diploid range, it's type is monoploid range typei
    // ipddev, ipd, ipd homo etc are all diploid range with redundant entries so using typei vs chi should not make a difference
    // monitor ipd takes diploid range entries and searches through all diploid ids to match types
    // constraint mat is (monoploid+1)^2 size matrix so needs to have monoploid+1 range indices
    int typei=box->chromType[chi];
//     printf("check redundancy %d %d --homoi %f %f -- ipdij %f %f -- ipddevij %f %f\n",chi,typei,box->ipdhomo[chi][0],box->ipdhomo[typei][0],box->ipd[chi][chi+1],box->ipd[typei][typei+1],box->ipddev[chi][typei+1],box->ipddev[typei][typei+1]);
    for(jj=0;jj<p->nChrom;jj++)
    {
        // jj and type of jj should be same
//         printf("%d %d %d\n",chi,jj,box->chromType[jj]);
//         if(jj==typei)
//             continue;

        // calculate the deviation of the average distance of this pair of homologues from each other

        if( box->ipdhomo[chi][0] != jj && box->ipdhomo[jj][0] != chi) // not homologue
        {
            netspring = 0.5*p->k_ipd/(box->ipddev[chi][jj]*box->ipddev[chi][jj]);

            // monitorIPD takes diploid chom id
            r = monitorIPD(p,box,chi,jj); // distance is scaled to chromosome sizes

            // update the constraint matrix
//             printf("cmat %d %d %f %f\n",typei,jj,box->constraintMat[typei][jj],box->constraintMat[jj][typei]);
            box->constraintMat[typei][jj] = r;
            box->constraintMat[jj][typei] = r;

            en += netspring*(r - box->ipd[chi][jj])*(r - box->ipd[chi][jj]);
        }
    }
    return en;

}

double chromHomologueBias(struct PARAINPUT *p,struct SYSTEM *box,int chi)
{
    // homologue contribution to energy
    // chi is diploid range index
    // typei is monoploid range type specifier
    // monitorHomologues takes diploid range entry as argument
    // constraint mat takes monoploid+1 range indices
    int typei = box->chromType[chi];
    double r,en=0.0;
    r = monitorHomologues(p,box,chi); // distance is scaled to chromosome sizes
//     printf("check homologue entries %d %d - %f %f -- %f %f\n",chi,typei,box->ipdhomo[chi][0],box->ipdhomo[typei][0],box->ipdhomo[chi][1],box->ipdhomo[typei][1]);
    en += (0.5*p->k_homologue)*sqr(r - box->ipdhomo[chi][1]);
    box->constraintMat[typei][typei] = r;

    return en;
}

double chromCentroidBias(struct PARAINPUT *p,struct SYSTEM *box,int chi)
{
    // centroid contribution to energy
    // chi is diploid range entry
    // typei is monoploid range type specifier
    // monitor centroid takes diploid range entry and searches through diploid indices to return the separation from centroid
    // constraint mat takes monoploid+1 range indices
    int typei = box->chromType[chi];
    double r,en=0.0;
    r = monitorCentroid(p,box,chi);
//     printf("check symm centroid ch type %d %d -- dist %f %f dist dev %f %f\n",chi,typei,box->chromDist[chi],box->chromDist[typei],box->chromDistDev[chi],box->chromDistDev[typei]);
    en += 0.5*p->k_centroid/(box->chromDistDev[chi]*box->chromDistDev[chi])*sqr(r - box->chromDist[chi]);
    box->constraintMat[typei][p->nChrom] = r;
    box->constraintMat[p->nChrom][typei] = r;
    return en;
}


void calcDeviations(struct PARAINPUT *p, struct SYSTEM *box,int step)
{
    // calculate the cosine similarity and the L_infty
    double linf=0.0,cosine_sim=0.0;
    double compdist;
    int ii,jj;
    double maxcompdist=0.0;
    double normconstraint,normexpected=0.0;
    int maxdeviant=0,maxi=0,maxj=0;
    int constrcount=0;
    int compdistbin=0;
    double deltacd = 0.0;
    double totDev=0.0;
    double dotab=0.0;
    box->devdir = 0.0;

    for(ii=0;ii<p->nChrom+1;ii++)
	{
		for(jj=0;jj<p->nChrom+1;jj++)
		{
            if(jj>=ii && box->toCalc[ii][jj] == 1) // && ii != p->nChrom-1 && jj != p->nChrom-1) // condition to ensure that only constraints to be calculated are considered
            {
                if(ii==p->nChrom-1 || jj==p->nChrom-1)
                {
                    printf("%d %d %d %f\n",ii,jj,box->toCalc[ii][jj],box->constraintMat[ii][jj]);
                }

                constrcount = (p->nChrom+1)*ii - (int)(ii*(ii-1)/2) + jj - ii;
//                 printf("testing constrcount %d %d %d\n",ii,jj,constrcount);
//                 compdist = fabs(box->constraintMat[ii][jj] - box->expectedMat[ii][jj]);
                compdist = fabs(distCalc(box->constraintMat[ii][jj],box->expectedMat[ii][jj]));
                totDev += fabs(box->constraintMat[ii][jj] - box->expectedMat[ii][jj]);
                // now compute cosine similarity -- dot product(a,b)/norm(a)norm(b)

                dotab += box->constraintMat[ii][jj]*box->expectedMat[ii][jj];

                if ( compdist > maxcompdist)
                {
                    maxcompdist = compdist;
                    box->devdir = box->constraintMat[ii][jj] - box->expectedMat[ii][jj];
                    maxdeviant = constrcount; // ii*(p->nChrom+1) + jj; // index of largest deviant
                    maxi = ii;
                    maxj = jj;
                }
                normconstraint+= box->constraintMat[ii][jj]*box->constraintMat[ii][jj];
                normexpected+= box->expectedMat[ii][jj]*box->expectedMat[ii][jj];



            }

        }

    }

    box->l_infty=maxcompdist;
//     printf("constr count %d should be %d\n",constrcount+1,(int)((p->nChrom+1)*(p->nChrom+2)/2));
    box->l1=totDev; ///(constrcount+1);
    box->maxdev=maxdeviant;
    box->maxi = maxi;
    box->maxj = maxj;
    normconstraint = sqrt(normconstraint);
    normexpected = sqrt(normexpected);
    cosine_sim = dotab/(normconstraint*normexpected);
    box->cosine_sim = cosine_sim;

//     // reset L1 histo
//
//     for(ii=0;ii<p->nbinl1histo;ii++)
//         box->l1histo[ii]=0.0;

    // and update the L1 histo and L_inf proximity tracker if beyond initial equilibration phase
    double proxLinf=0.0;
    if(step>(int)(0.2*p->MCIter))
    {
        for(ii=0;ii<p->nChrom+1;ii++)
        {
            for(jj=0;jj<p->nChrom+1;jj++)
            {
                if(jj>=ii && box->toCalc[ii][jj] == 1)
                {


//                     compdist = fabs(box->constraintMat[ii][jj] - box->expectedMat[ii][jj]);
                    compdist = fabs(distCalc(box->constraintMat[ii][jj],box->expectedMat[ii][jj]));
                    // bin the given constraint

                    deltacd = 1.0/(1.0*p->nbinl1histo);
                    compdistbin = (int) ((box->l_infty - compdist)/deltacd);
                    // update the histo

                    if(compdist < box->l_infty)
                        box->l1histo[compdistbin] += 1.0;

                    if(box->l_infty - compdist < 0.2*box->l_infty)
                        proxLinf +=1.0;



                }
            }
        }


        box->proxLinf = proxLinf;


    //  spot normalising factor of l1histo
        box->norml1histo=0.0;
        for(ii=0;ii<p->nbinl1histo;ii++)
        {
            box->norml1histo += box->l1histo[ii];
        }
//         printf("norm %f\n",box->norml1histo);
    }

}


void findParticipants(struct PARAINPUT *p, struct SYSTEM *box,int chi, int chj)
{
    // subroutine to update box->i1,i2,j1,j2 based on input chi,chj

    int typei,typej,ii,jj,i,j;
//     printf("finding participants for constr %d %d -- 1D %d\n",chi,chj,(p->nChrom+1)*chi - (int)(chi*(chi-1)/2) + chj - chi);
    if(chi<p->nChrom && chj<p->nChrom && (chi==chj || box->chromType[chi]==box->chromType[chj])) // homologues from diagonal -- neither is the centroid
    {
//         printf("Homologues!!! %d %d\n",chi,chj);
        typei = box->chromType[chi];
        // find the two homologues
        // NOTE: In this loop, chi1 (=chi2) is assigned the first homologue and chj1 (=chj2) is assigned the second
        // NOTE: This is done because the 4 pairwise distances we calculate are i1j1,i1j2,i2j1,i2j2 (never i1i2,jij2)
        box->chi1 = p->nDChrom;
        box->chj1 = -1;
        for(ii = 0; ii < p->nDChrom;ii++)
        {
            i = box->chromList[ii];
            if(box->chromType[ii] == typei) // true twice
            {
                box->chi1 = MIN(box->chi1,ii); // box->chi1 is first homologue -- will always be updated when first encountering typei, will not be updated at second encounter
                // update chj1 with changed chi1 -- first update chj1=chi1 at first encounter, second update chi2=ii at second encounter
                box->chj1 = MAX(box->chi1,ii);
            }
        }
        box->chi2 = box->chi1;
        box->chj2 = box->chj1;
    }
    else
    {
        if(chi<p->nChrom)
        {
            typei = box->chromType[chi];
            // find the two homologues
            box->chi1 = p->nDChrom;
            box->chi2 = -1;
            for(ii = 0; ii < p->nDChrom;ii++)
            {
                i = box->chromList[ii];
                if(box->chromType[ii] == typei) // true twice
                {
                    box->chi1 = MIN(box->chi1,ii); // box->chi1 is first homologue -- will always be updated when first encountering typei, will not be updated at second encounter
                    // update chi2 with changed chi1 -- first update chi2=chi1 at first encounter, second update chi2=ii at second encounter
                    box->chi2 = MAX(box->chi1,ii);
                }
            }
        }
        if(chi==p->nChrom) // chi marks centroid so is nChrom
        {
            box->chi1 = p->nDChrom; // not indexed in chromList
            box->chi2 = p->nDChrom;
        }
        if(chj<p->nChrom)
        {
            typej = box->chromType[chj];
            // find the two homologues
            box->chj1 = p->nDChrom;
            box->chj2 = -1;
            for(jj = 0; jj < p->nDChrom;jj++)
            {
                j = box->chromList[jj];
                if(box->chromType[jj] == typej) // true twice
                {
                    box->chj1 = MIN(box->chj1,jj); // box->chj11 is first homologue -- will always be updated when first encountering typej, will not be updated at second encounter
                    // update chi2 with changed chi1 -- first update chj2=chj1 at first encounter, second update chj2=jj at second encounter
                    box->chj2 = MAX(box->chj1,jj);
                }
            }
        }
        if(chj==p->nChrom)
        {
            box->chj1 = p->nDChrom; // not indexed in chromList
            box->chj2 = p->nDChrom;

        }
    }
//     printf("participants %d %d are %d %d %d %d\n",chi,chj,box->chi1,box->chi2,box->chj1,box->chj2);



}


void updateGrad(struct PARAINPUT *p, struct SYSTEM *box, int ii1, int jj1, int ii, int jj, int max)
{
    // test code to calc dx, dy,dz and update forces based on i1,j1 values -- i1,j1 are calculated based on chi1,chi2
    int i1,j1,j;
    double pxi1,pyi1,pzi1,pxj1,pyj1,pzj1;
    double dx,dy,dz,rij,invrij;
    double ff;
    double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	double surfcut=0.0;

    if(ii1 < p->nDChrom)
    {
        i1 = box->chromList[ii1];
        pxi1 = box->X[i1];
        pyi1 = box->Y[i1];
        pzi1 = box->Z[i1];
    }
    else
    {
        i1 = -1;
        pxi1 = 0.0;
        pyi1 = 0.0;
        pzi1 = 0.0;
    }

    if(jj1 < p->nDChrom)
    {
        j1 = box->chromList[jj1];
        pxj1 = box->X[j1];
        pyj1 = box->Y[j1];
        pzj1 = box->Z[j1];
    }
    else
    {
        j1 = -1;
        pxj1 = 0.0;
        pyj1 = 0.0;
        pzj1 = 0.0;
    }

    dx = 0.5*(pxi1 - pxj1)/p->ella;
    dy = 0.5*(pyi1 - pyj1)/p->ellb;
    dz = 0.5*(pzi1 - pzj1)/p->ellc;

    rij = sqrt(dx*dx + dy*dy + dz*dz);
    if(rij > 0.0)
        invrij = 1.0/rij;
    else
        invrij = 0.0;

    /// use the toCalc array to decide if force calc wrt this pair of constraints needs to be made
    // this is sufficient when gradients for all chromosomes are calculated, with a simultaneous optimisation update
    // for the max gradient optimisation, we need to ensure that toCalc flag has been used while selecting the maxi and maxj constraints to use here
    // making this check below redundant when max==1 (but we leave it in anyway for a sanity check
    ff = 0.0;
//     if(box->toCalc[ii][jj] == 0 && max ==1)
//     {
//         printf("%d %d -- mat entry %d %d\n",ii,jj,box->toCalc[ii][jj],box->toCalc[jj][ii]);
//         printf("this should not ever ever happen\n");
//         exit(0);
//     }
    if (box->toCalc[ii][jj] ==1)
    {
//     ff = 0.5*p->k_ipd*(rij - box->expectedMat[ii][jj]); // 2.0/4.0 because averaged over 4
    ff = 0.5*p->k_ipd*distCalc(rij,box->expectedMat[ii][jj]);
    }
        // third, update the forces to the particle indices i12,j12 based on the individual gradients and the individual distances

    if(i1 >=0)
    {
        box->fx[i1] += -ff*dx*invrij;
        box->fy[i1] += -ff*dy*invrij;
        box->fz[i1] += -ff*dz*invrij;
    }
    if(j1 >=0)
    {
        box->fx[j1] += ff*dx*invrij;
        box->fy[j1] += ff*dy*invrij;
        box->fz[j1] += ff*dz*invrij;
    }

    // calculate surface gradient contribution for i1 and j1 -- this is always calculated
//     double r=0.0;
//     double patchcut=0.0;
    for(j=0;j<p->nSurfPoints;j++)
    {

        // first for i1
        dx = box->surfPointx[j] - pxi1;
        dy = box->surfPointy[j] - pyi1;
        dz = box->surfPointz[j] - pzi1;

//         rij = sqrt(dx*dx + dy*dy + dz*dz);
        rij = dx*dx + dy*dy + dz*dz;

        // check if distance within WCA cut-off
        surfcut = 0.5*avesep + box->radius[i1];
//         printf("chrom %d -- dist %f - cut %f pos %f %f %f -- surf %f %f %f - r %f - pdps %d %f\n",i1,sqrt(rij),surfcut,pxi1,pyi1,pzi1,
//               box->surfPointx[j],box->surfPointy[j],box->surfPointz[j],box->radius[i1],j,avesep);
        ff=0.0;
//         if(rij <= sqr(pow(2.0,1.0/6.0)*surfcut))
//         {
//             rij = sqrt(rij);
//             invrij = 1.0/rij;
// //             en += 4.0*( pow(surfcut/r,12.0) - pow(surfcut/r,6.0)) + 1.0;
//             ff = -4.0*(-12.0*pow(surfcut/rij,12.0)*invrij + 6.0*pow(surfcut/rij,6.0)*invrij);
//         }
        if(rij <= surfcut*surfcut)
        {
            rij = sqrt(rij);
            invrij = 1.0/rij;
            ff = p->eps_hertz*2.5*pow((1 - rij/surfcut),1.5)/surfcut;
        }
        if(i1>=0)
        {
        box->fx[i1] += -ff*dx*invrij;
        box->fy[i1] += -ff*dy*invrij;
        box->fz[i1] += -ff*dz*invrij;
        }

        // now for j1

        dx = box->surfPointx[j] - pxj1;
        dy = box->surfPointy[j] - pyj1;
        dz = box->surfPointz[j] - pzj1;

//         rij = sqrt(dx*dx + dy*dy + dz*dz);
        rij = dx*dx + dy*dy + dz*dz;

        // check if distance within WCA cut-off
        surfcut = 0.5*avesep + box->radius[j1];
//         printf("chrom %d -- dist %f - cut %f pos %f %f %f -- surf %f %f %f - r %f - pdps %d %f\n",j1,sqrt(rij),surfcut,pxj1,pyj1,pzj1,
//               box->surfPointx[j],box->surfPointy[j],box->surfPointz[j],box->radius[j1],j,avesep);
        ff=0.0;
//         if(rij <= sqr(pow(2.0,1.0/6.0)*surfcut))
//         {
//             rij = sqrt(rij);
//             invrij = 1.0/rij;
// //             en += 4.0*( pow(surfcut/r,12.0) - pow(surfcut/r,6.0)) + 1.0;
//             ff = -4.0*(-12.0*pow(surfcut/rij,12.0)*invrij + 6.0*pow(surfcut/rij,6.0)*invrij);
//         }
        if(rij <= surfcut*surfcut)
        {
            rij = sqrt(rij);
            invrij = 1.0/rij;
            ff = p->eps_hertz*2.5*pow((1 - rij/surfcut),1.5)/surfcut;
        }
        if(j1>=0)
        {
        box->fx[j1] += -ff*dx*invrij;
        box->fy[j1] += -ff*dy*invrij;
        box->fz[j1] += -ff*dz*invrij;
        }
    }
}

void calcAllGradients(struct PARAINPUT *p, struct SYSTEM *box, int max)
{
    // subroutine to calculate the gradient of each assuming harmonic deviation
    int ii,jj,i1,j1,i2,j2,i;
    int constrcount=0;
    int maxdeviant=0,maxi=0,maxj=0;
    double dx,dy,dz,rij,invrij;
    double compdist=0.0;
    double maxcompdist=0.0;
    double random;
    box->devdir = 0.0;

    // reset forces
    for(ii=0;ii<p->nAto;ii++)
    {
        box->fx[ii]=0.0;
        box->fy[ii]=0.0;
        box->fz[ii]=0.0;
    }
//     printf("reset forces in gradient calc - %d\n",max);
    for(ii=0;ii<p->nChrom;ii++)
	{
		for(jj=0;jj<p->nChrom+1;jj++)
		{
            if(jj>=ii)
            {

                if(box->toCalc[ii][jj]==0) // || ii == p->nChrom-1 || jj == p->nChrom-1) // skip the X chrom for gradient descent
                    continue;
                constrcount = (p->nChrom+1)*ii - (int)(ii*(ii-1)/2) + jj - ii;
//                 printf("testing constrcount %d %d %d\n",ii,jj,constrcount);


                // first, find the i1,i2,j1,j2 that correspond to ii,jj -- generalise for either of ii,jj to be centroid, homologue etc
                findParticipants(p,box,ii,jj);
                // findParticipants assigns box->chij12 (4 vars) to either [0,NDChrom) or NDChrom if centroid
                // i1=i2 and/or j1=j2 if homologue distance
                //

                // NOTE: Code block below will have to be replicated -- moved to assignPositions

                // second, the distance between the pair i1,j1, i2,j2, gives the gradient contribution -- the gradient of algebraic sum is equivalent to algebraic sum of gradients
//                 assignPositions(p,box);

                /// replicate code block below for each of the contributions -- this contri i1,j1

                // NOTE: Subroutine updateGrad takes a given combination and updates the gradient contribution on the corresponding chromosomes (0 update for centroid)
                // we therefore call updateGrad 4 times for each of the 4 contributions
                // note that in the case of homologues and centroid distances we will have either 1 or 2 redundant pairs but the weight of 1/4 on the algrbraic sum of gradients accounts

                if(max==0 && box->toCalc[ii][jj] == 1) // update all gradients
                {
                    updateGrad(p,box,box->chi1,box->chj1,ii,jj,max);
                    updateGrad(p,box,box->chi1,box->chj2,ii,jj,max);
                    updateGrad(p,box,box->chi2,box->chj1,ii,jj,max);
                    updateGrad(p,box,box->chi2,box->chj2,ii,jj,max);


                }

                /// now calc max gradient
                compdist = fabs(box->constraintMat[ii][jj] - box->expectedMat[ii][jj]);
                if( compdist > maxcompdist && box->toCalc[ii][jj] == 1) // need to make sure maxi and maxj are never constraints that are flagged out of consideration
                {
                    maxcompdist = compdist;
                    box->devdir = box->constraintMat[ii][jj] - box->expectedMat[ii][jj];
                    maxdeviant = constrcount; // ii*(p->nChrom+1) + jj; // index of largest deviant
                    maxi=ii;
                    maxj=jj;
                }
            }
        }
    }
    box->maxi = maxi;
    box->maxj = maxj;
//     printf("%d %d max constraints\n",maxi,maxj);
    // max==1 means only the gradients of chromosomes participating in the largest deviant are updated -- forces are set to zero for everything else
    // note that imposing the toCalc flag on the constraints while identifying maxi and maxj ensure that unwanted constraints are never identified as the max
    // this is done despite the fact that updateGrad has a similar condition in calculating force gradients because setting force gradient to zero there is insufficient
    // to ensure that only valid, important constraints are considered while identifying the max gradient
    if(max==1)
    {


        if(box->toCalc[maxi][maxj] == 1)
        {
        findParticipants(p,box,maxi,maxj);

        random=drand48();

//         if(random>0 && random < 0.25)
        updateGrad(p,box,box->chi1,box->chj1,maxi,maxj,max);

//         if(random>=0.25 && random < 0.5)
        updateGrad(p,box,box->chi1,box->chj2,maxi,maxj,max);

//         if(random>=0.5 && random < 0.75)
        updateGrad(p,box,box->chi2,box->chj1,maxi,maxj,max);

//         if(random>=0.75 && random < 1)
        updateGrad(p,box,box->chi2,box->chj2,maxi,maxj,max);
        }
    }

//     for(ii=0;ii<p->nDChrom;ii++)
//     {
//         i = box->chromList[ii];
//         if(ii==box->chi1 || ii==box->chj1 || ii==box->chi2 || ii==box->chj2)
// //         printf("%d %d -ch %d %d %d %d -- f-- %f %f %f -- p --%f %f %f\n",ii,i,box->chi1,box->chi2,box->chj1,box->chj2,box->fx[i],box->fy[i],box->fz[i],box->X[i],box->Y[i],box->Z[i]);
//     }

}


void calcBias(struct PARAINPUT *p, struct SYSTEM *box)
{
    int ii;

    box->currW = IPDBias(p,box,1);

    // calcDeviations decides if a particular constraint can contribute to l_infty based on toCalc
	calcDeviations(p,box,0);

    // only l_inf deviation en is retained for currW

	box->currW = p->k_ipd*box->l_infty;

	for(ii=0;ii<p->nDChrom;ii++)
	{
		box->currW += surfEnergy(p,box,ii);
	}

}
