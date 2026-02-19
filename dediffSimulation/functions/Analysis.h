void getProxyHiC(struct PARAINPUT *p,struct SYSTEM *box)
{
    // this subroutine updates a proxy HiC matrix based on sampled proximity between components of each chromosome

    int i,ii,j,jj,ci,cj,ii1,i1,jj1,j1;
    double diffX,diffY,diffZ;
    double rij,rijsq,invrij,r0,ff,patchcut,separation;

    for(ii=0;ii<p->nChrom-1;ii++)
    {
        i = box->chromList[ii];
        for(jj=ii+1;jj<p->nChrom;jj++)
        {
            j = box->chromList[jj];

            // check they are not same type -- box->type[i] == chromtype for chrom (0 for patch -- use as sanity check)

            if(box->type[i] == box->type[j])
                continue;

            ci = box->type[i]-1;
            cj = box->type[j]-1;
//             printf("sanity check chrom %d - %d (raw %d - %d) and type %d - %d\n",ii,jj,i,j,ci,cj);
            // check separation chroms i and j

            diffX=box->X[j]-box->X[i];

			diffY=box->Y[j]-box->Y[i];

			diffZ=box->Z[j]-box->Z[i];

			rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;

            r0 = box->radius[i] + box->radius[j];

            separation = sqrt(rijsq/(r0*r0));
            // NOTE: Here we can introduce an activation function for separation
            // NOTE: For now we are taking a tanh function


            // assign matrix entry to (ci,cj) entry of HiC matrix
//             printf("sizes %d\n",(int)(0.5*p->nChrom));
//             box->HiCMat[ci][cj] += 0.5 - 0.5*tanh((separation-1)/(0.5*r0));
//             box->HiCMat[cj][ci] += 0.5 - 0.5*tanh((separation-1)/(0.5*r0));


            // now loop over all patches of chrom i and chrom j

            for(ii1=1;ii1<=box->patchList[ii][0];ii1++)
            {
                i1 = box->patchList[ii][ii1];

                if(box->lamin[i1]>0)
                    continue;

                for(jj1=1;jj1<=box->patchList[jj][0];jj1++)
                {
                    j1 = box->patchList[jj][jj1];
                    if(box->lamin[j1]>0)
                    continue;
//                     printf("second sanity check -- %d %d %d with %d %d %d type %d %d\n",ii,ii1,i1,jj,jj1,j1,ci,cj);

                    // find separation

                    diffX=box->X[j1]-box->X[i1];

                    diffY=box->Y[j1]-box->Y[i1];

                    diffZ=box->Z[j1]-box->Z[i1];

                    rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;

                    r0 = box->radius[i1] + box->radius[j1];

                    separation = sqrt(rijsq) ; // /(r0*r0));

                    box->HiCMat[ci][cj] += 0.5 - 0.5*tanh((separation-p->hic_cut)/(p->hic_decay*p->hic_cut));
                    box->HiCMat[cj][ci] += 0.5 - 0.5*tanh((separation-p->hic_cut)/(p->hic_decay*p->hic_cut));
                }
            }
        }
    }
//     printf("done computing HiC\n");

}

void getIPD(struct PARAINPUT *p,struct SYSTEM *box)
{
    // subroutine to calculate the IPD matrix output based on the centroids of each chrom
    int ii,jj,nhap,constrcount,i,j,ci,cj;
    int ci1,ci2,cj1,cj2,i1,i2,j1,j2;

    double diffX,diffY,diffZ,avgdist;

    nhap = (int)(0.5*p->nChrom);

    for(ii=0;ii<nhap+1;ii++)
    {

        for(jj=ii;jj<nhap+1;jj++)
        {
            constrcount = (nhap+1)*ii - (int) (ii*(ii-1)/2) + jj - ii;
            avgdist = 0.0;

            // check if ii is centroid
            if(ii==nhap)
            {
                if(jj==ii)
                    box->constraintMat[constrcount] = 0.0;
                continue;
            }

            // check for ii chrom, jj == ii (diploid)
            if(ii!=nhap && ii==jj)
            {
                ci = ii;
                cj = nhap+ii;
                i = box->chromList[ci];
                j = box->chromList[cj];

                diffX=box->X[j]-box->X[i];
                diffY=box->Y[j]-box->Y[i];
                diffZ=box->Z[j]-box->Z[i];

//                 rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;

                // calculate the scaled distance ( current geom)
                avgdist = sqrt(sqr(0.5*diffX/p->ella) + sqr(0.5*diffY/p->ellb) + sqr(0.5*diffZ/p->ellc));

                box->constraintMat[constrcount] = avgdist;
            }

            // check centroid dist
            if(ii!=nhap && jj==nhap)
            {
                ci1 = ii;
                ci2 = nhap + ii;
                i = box->chromList[ci1];
                j = box->chromList[ci2];

                diffX = box->X[i];
                diffY = box->Y[i];
                diffZ = box->Z[i];

                avgdist = sqrt(sqr(0.5*diffX/p->ella) + sqr(0.5*diffY/p->ellb) + sqr(0.5*diffZ/p->ellc));

                diffX = box->X[j];
                diffY = box->Y[j];
                diffZ = box->Z[j];

                avgdist += sqrt(sqr(0.5*diffX/p->ella) + sqr(0.5*diffY/p->ellb) + sqr(0.5*diffZ/p->ellc));

                box->constraintMat[constrcount] = 0.5*avgdist;
            }

            // check IPD dist
            if(ii!=nhap && jj!=nhap && ii!=jj)
            {
                ci1 = ii;
                ci2 = nhap + ii;

                cj1 = jj;
                cj2 = nhap + jj;

                i1 = box->chromList[ci1];
                i2 = box->chromList[ci2];

                j1 = box->chromList[cj1];
                j2 = box->chromList[cj2];

                // now compute 4 distance combinations

                diffX=box->X[j1]-box->X[i1];
                diffY=box->Y[j1]-box->Y[i1];
                diffZ=box->Z[j1]-box->Z[i1];
                avgdist = sqrt(sqr(0.5*diffX/p->ella) + sqr(0.5*diffY/p->ellb) + sqr(0.5*diffZ/p->ellc));

                diffX=box->X[j2]-box->X[i1];
                diffY=box->Y[j2]-box->Y[i1];
                diffZ=box->Z[j2]-box->Z[i1];
                avgdist += sqrt(sqr(0.5*diffX/p->ella) + sqr(0.5*diffY/p->ellb) + sqr(0.5*diffZ/p->ellc));

                diffX=box->X[j1]-box->X[i2];
                diffY=box->Y[j1]-box->Y[i2];
                diffZ=box->Z[j1]-box->Z[i2];
                avgdist += sqrt(sqr(0.5*diffX/p->ella) + sqr(0.5*diffY/p->ellb) + sqr(0.5*diffZ/p->ellc));

                diffX=box->X[j2]-box->X[i2];
                diffY=box->Y[j2]-box->Y[i2];
                diffZ=box->Z[j2]-box->Z[i2];
                avgdist += sqrt(sqr(0.5*diffX/p->ella) + sqr(0.5*diffY/p->ellb) + sqr(0.5*diffZ/p->ellc));


                box->constraintMat[constrcount] = 0.25*avgdist;
            }
        }
    }

    // relegate analysis to analysis scripts
}
