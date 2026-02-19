double gaussVar(double mu,double sig)
{
    double gauss_r,v1,v2,gauss_var;
    gauss_r=2.0;
    while(gauss_r>=1.0)
    {
//         v1=2.0*genrand_real1(&(*a))-1.0;
//         v2=2.0*genrand_real1(&(*a))-1.0;
		v1 = 2.0*drand48()-1.0;
		v2 = 2.0*drand48()-1.0;
//         printf("v1 v2 %f %f\n",v1,v2);
        gauss_r = v1*v1 + v2*v2;
    }
//     printf("gr log gr %f %f\n",gauss_r,log(gauss_r));
    gauss_var=v1*sqrt(-2.0*log(gauss_r)/gauss_r);
    gauss_var=mu+sig*gauss_var;
    
    return gauss_var;
}

void readChromInfo(FILE *fp,struct PARAINPUT *p,struct SYSTEM *box)
{
	// Expected Chrom_props columns:
	// npatch mass chrom_type ec_fraction hc_fraction lamin_fraction
	int i,j,npatch,ctype;
	double spamm,spamhc,spamec,spamlamin;
	int globCount=0;
	int nlines=0;
	char line[1024];
	int npatch_rows[512];
	int ctype_rows[512];
	double mass_rows[512];
	double ec_rows[512];
	double hc_rows[512];
	double lamin_rows[512];
	printf("reading chrom properties file\n");

	while(fgets(line,sizeof(line),fp)!=NULL)
	{
		int vals;
		if(line[0]=='#' || line[0]=='\n')
			continue;
		vals = sscanf(line,"%d %lf %d %lf %lf %lf",&npatch_rows[nlines],&mass_rows[nlines],&ctype_rows[nlines],&ec_rows[nlines],&hc_rows[nlines],&lamin_rows[nlines]);
		if(vals!=6)
		{
			printf("Chrom_props parse error at line %d: expected 6 columns\n",nlines+1);
			exit(0);
		}
		nlines++;
	}

	if(nlines != p->nChrom && nlines != p->nDChrom)
	{
		printf("expected %d or %d lines in chrom properties file but got %d\n",p->nChrom,p->nDChrom,nlines);
		exit(0);
	}

	int patchcountcheck=0;
	int patchpointer=0;
	int exit_flag=0;
	for(i=0;i<p->nDChrom;i++)
	{
		int src_idx = (nlines==p->nDChrom) ? i : (i % p->nChrom);
		npatch = npatch_rows[src_idx];
		spamm = mass_rows[src_idx];
		ctype = ctype_rows[src_idx];
		spamec = ec_rows[src_idx];
		spamhc = hc_rows[src_idx];
		spamlamin = lamin_rows[src_idx];

		box->chromList[i] = globCount;

		box->chrompatchmass[i] = spamm*spamec/(1.0*npatch);
		box->chromcoremass[i] = spamm*spamhc;
		box->chromlaminmass[i] = spamm*spamlamin;

		box->mass[globCount]=box->chromcoremass[i];
		box->radius[globCount]=pow(box->chromcoremass[i]/p->massscalefact,0.3333);

		box->whichChrom[globCount]=i;
		box->whichPatch[globCount]=-1;
		box->chromType[i] = ctype - 1;
		// Centroid constraints are inactive; keep placeholders for file-shape compatibility.
		box->chromDist[i] = 0.0;
		box->chromDistDev[i] = 1.0;
		globCount++;

		if(patchcountcheck > 2*p->nPatchTot)
		{
			printf("danger danger max %d (non-lamin) but given %d -- so far %d chrom\n",2*p->nPatchTot,patchcountcheck,i+1);
			exit_flag=1;
		}
		patchcountcheck += npatch;
		// Patch arrays use index 0 to store count (+1 for lamin patch).
		box->patchposx[i][0]=npatch+1;
		box->patchposy[i][0]=npatch+1;
		box->patchposz[i][0]=npatch+1;
		box->patchList[i][0]=npatch+1;
		box->patchPointer[i][0]=npatch+1;

		if(npatch != p->nChrom-1)
		{
			printf("Mismatch between patch count and number of chromosome types for chrom %d -- expected %d patches based on %d chrom types but got (%d,%d) respectively\n",i,p->nChrom-1,p->nChrom,npatch,p->nChrom);
			exit_flag=1;
		}

		// Populate each chromosome's patch and lamin entries consumed by SimStart/IntMat output.
		for(j=1;j<=box->patchList[i][0];j++)
		{
			box->patchList[i][j]=globCount;
			box->whichChrom[globCount]=i;
			box->whichPatch[globCount]=j;

			if(j<box->patchList[i][0])
			{
				box->lamin[i][j]=0;
				box->mass[globCount]=box->chrompatchmass[i];
				box->radius[globCount]=pow(box->chrompatchmass[i]/p->massscalefact,0.3333);
			}
			else
			{
				box->lamin[i][j]=1;
				box->mass[globCount]=box->chromlaminmass[i];
				box->radius[globCount]=pow(box->chromlaminmass[i]/p->massscalefact,0.3333);
			}
			globCount++;
		}

	}
	for(i=0;i<p->nDChrom;i++)
	{
		patchpointer=1;
		for(j=0;j<p->nChrom;j++) // only for each chrom type
		{
			// loop over chromosomes to assign the patch identity
			if(box->chromType[j] != box->chromType[i]) // skip self
			{
// 				printf("chrom %d is type %d -- patch %d (%d) points to %d\n",i,box->chromType[i],j,patchpointer,box->chromType[j]);
				box->patchPointer[i][patchpointer] = box->chromType[j];
				patchpointer++;
			}

		}
		box->patchPointer[i][box->patchList[i][0]]=p->nDChrom+1;
	}
	if(exit_flag==1)
	{
		printf("Encountered pathology, exiting....\n");
		exit(0);
	}
}

void readIPDInfo(FILE *fp,struct PARAINPUT *p,struct SYSTEM *box)
{
	// subroutine to read pairwise distances between chromosomes, as well as the stdev
	int spamci1,spamci2;
	double spamdist,spamdistdev;
	printf("reading IPD matrix\n");
	spamdistdev=p->spamipddev;
	int last_read=0;
	while(!feof(fp))
	{
		fscanf(fp,"%d %d %lf\n",&spamci1,&spamci2,&spamdist);

		box->ipd[spamci1-1][spamci2-1] = p->IPD_scale*spamdist;
		box->ipddev[spamci1-1][spamci2-1] = p->IPD_scale*spamdistdev; //*spamdist;

		box->ipd[spamci2-1][spamci1-1] = p->IPD_scale*spamdist;
		box->ipddev[spamci2-1][spamci1-1] = p->IPD_scale*spamdistdev; //*spamdist;

		if(spamci1==23)
			box->ipd[spamci1-1][spamci2-1] = 0.0;
		if(spamci1==23)
			box->ipd[spamci2-1][spamci1-1] = 0.0;

		// specify entries for the second copy as well

		last_read=spamci2-1;

	}


	// last_read+1 is the X/Y chromosome -- it's entry is random with a large shell



	int homo_i;
	int homo_j;
	int i,j;
	// one unread chrom
	if(last_read < p->nChrom - 1)
	{
		// check here for consistency that 2*(last_read+2) == size of chromList -- nChrom etc.
		if(2*(last_read+2) != p->nDChrom)
		{
			printf("error with IPD matrix size and (diploid) chromosome number mismatch, expecting 2x%d to be equal to %d but getting %d\n",last_read+2,p->nDChrom,2*(last_read+2));
			exit(0);
		}
		for(i=0;i<=last_read;i++)
		{
// 			printf("ipd XY at %d for %d\n",last_read+1,i);
			box->ipd[last_read+1][i] = 0.35; // box->ipd[i][6]; //0.25;
			box->ipd[i][last_read+1] = 0.35; //box->ipd[i][6]; //0.25;

			box->ipddev[last_read+1][i] = 0.2;
			box->ipddev[i][last_read+1] = 0.2;
		}
		last_read+=1;
	}
	else
	{
		if(2*(last_read+1) != p->nDChrom)
		{
			printf("error with IPD matrix size and (diploid) chromosome number mismatch, expecting 2x%d to be equal to %d but getting %d\n",last_read+1,p->nDChrom,2*last_read+2);
			exit(0);
		}
	}

	// last_read+2-- 2*(last_read+2) -1 -- indices of other chroms, loop over and use IPD values
	for(i=last_read+1;i<p->nDChrom;i++)
	{


		homo_i = i - last_read - 1;
// 		printf("fix ipd of %d, homo %d with ref split at %d\n",i,homo_i,last_read);

		for(j=last_read+1;j<p->nDChrom;j++)
		{
			homo_j = j - last_read - 1;
// 			printf("with %d, homo %d\n",j,homo_j);
			if(j>i)
			{

				box->ipd[i][homo_j] = box->ipd[homo_i][homo_j];
				box->ipd[i][j] = box->ipd[homo_i][homo_j];
				box->ipd[homo_i][j] = box->ipd[homo_i][homo_j];

				box->ipd[homo_j][i] = box->ipd[homo_j][homo_i];
				box->ipd[j][i] = box->ipd[homo_j][homo_i];
				box->ipd[j][homo_i] = box->ipd[homo_j][homo_i];

				box->ipddev[i][homo_j] = box->ipddev[homo_i][homo_j];
				box->ipddev[i][j] = box->ipddev[homo_i][homo_j];
				box->ipddev[homo_i][j] = box->ipddev[homo_i][homo_j];

				box->ipddev[homo_j][i] = box->ipddev[homo_j][homo_i];
				box->ipddev[j][i] = box->ipddev[homo_j][homo_i];
				box->ipddev[j][homo_i] = box->ipddev[homo_j][homo_i];

			}


		}

	}



}

void readIPDHomo(FILE *fp,struct PARAINPUT *p,struct SYSTEM *box)
{
	int spamci1;
	double spam;
	int last_read=0;
	while(!feof(fp))
	{
		fscanf(fp,"%d %lf\n",&spamci1,&spam);
// 		printf("read entries %d %f\n",spamci1,spam);
		box->ipdhomo[spamci1-1][0]=0;
		box->ipdhomo[spamci1-1][1]=p->IPD_scale*spam;

		last_read = spamci1-1;
	}

	// execute check nChrom == 2*(last_read+2)
	if(2*(last_read+1) != p->nDChrom)
	{
		printf("error with IPD Homo size and (diploid) chromosome number mismatch, expecting 2x%d to be equal to %d but getting %d\n",last_read+2,p->nChrom,2*last_read+2);
		exit(0);
	}
// 	box->ipdhomo[last_read+1][0] = p->nChrom-1;
// 	box->ipdhomo[last_read+1][1] = 0.5;
	last_read += 1;
	// loop over and assign ipdhomo first entry and second entry for the rest of the array
	int i,homo_i;
	for(i=0;i<p->nDChrom;i++)
	{
		if(i<last_read) // 23
		{
			homo_i = i + last_read;
// 			printf("small i %d %d -- %d\n",i,homo_i,last_read);
			box->ipdhomo[i][0] = homo_i;
			// ipd_homo[i][1] unchanged

			box->ipdhomo[homo_i][0] = i;
			box->ipdhomo[homo_i][1] = box->ipdhomo[i][1];
		}
		// now check
		if(i>= last_read) //23
		{
			homo_i = i - last_read;
// 			printf("large i %d %d\n",i,homo_i);
			if(box->ipdhomo[i][0] != homo_i || box->ipdhomo[homo_i][0] != i)
			{
// 				printf("something wrong with homologous chrom pairing -- expected %d for %d but got %d\n",homo_i,i,(int)box->ipdhomo[i][0]);
				exit(0);
			}

		}
	}


}

int readIntMatrix(FILE *fp,struct PARAINPUT *p,struct SYSTEM *box)
{
	// subroutine to read out the NxM matrix of interactions
	// first read chromi1 patchj1 and chromi2 patchj2 then read the epsilon value
	// interaction matrix is 0 for all others
	int spamci1,spamci2,spampj1,spampj2,j1,j2,i,j;
	double spameps1,spameps2,maxeps,mineps;
	int npairs=0;
	printf("we ok pre-rewind\n");
// 	rewind(fp);
	printf("reading interaction matrix\n");
	maxeps = 0.0;
	mineps = 10.0;
	while(!feof(fp))
	{
// 		fscanf(fp,"%d %d %lf %lf\n",&spamci1,&spamci2,&spameps1,&spameps2);
		fscanf(fp,"%d %d %lf\n",&spamci1,&spamci2,&spameps1);
// 		printf("read int mat for %d %d as %f\n",spamci1,spamci2,spameps1);
		if(spameps1<0.0)
			spameps1=0.0;
		if(spamci1 <=0 || spamci2 <=0)
		{
			printf("chrom indices must be 1 or, will be mod in code to index correctly, higher check interactions file\n");
			exit(0);
		}
		// spamci1 and spamci2 specify the chromosome-chromosome interaction

		// we need to map these numbers to patch numbers

		// revert to indices from read numbers
		spamci1 -=1;
		spamci2 -=1;

		// NOTE: Need to work out what to do with X/Y chromosome patch energies

// 		spampj1 = box->patchList[spamci1][spamci2+1];
// 		spampj2 = box->patchList[spamci2][spamci1+1];

// 		spampj1 = box->patchPointer[spamci1][spamci2+1];
// 		spampj2 = box->patchPointer[spamci2][spamci1+1];

		spampj1 = spamci2;
		spampj2 = spamci1;

// 		printf("main %d %d %d %d\n",spamci1,spamci2,spampj1,spampj2);
		box->eps_ij[spampj1][spampj2] = spameps1;
		box->eps_ij[spampj2][spampj1] = spameps1;


// 		printf("%d,%d -- 1x transformed %d,%d -- raw %d %d\n",spamci1,spamci2,spamci1,spamci2,spampj1,spampj2);

// 		spampj1 = box->patchList[spamci1][spamci2+24];
// 		spampj2 = box->patchList[spamci2+23][spamci1+1];

// 		printf("half diag %d %d %d %d\n",spamci1,spamci2,spampj1,spampj2);
// 		box->eps_ij[spampj1][spampj2] = spameps1;
// 		box->eps_ij[spampj2][spampj1] = spameps1;
// // 		printf("%d,%d -- 2x transformed %d,%d -- raw %d %d\n",spamci1,spamci2,spamci1,22+spamci2,spampj1,spampj2);
//
// 		spampj1 = box->patchList[spamci1+23][spamci2+1];
// 		spampj2 = box->patchList[spamci2][spamci1+24];
// 		printf("other half diag %d %d %d %d\n",spamci1,spamci2,spampj1,spampj2);
// 		box->eps_ij[spampj1][spampj2] = spameps1;
// 		box->eps_ij[spampj2][spampj1] = spameps1;
// // 		printf("%d,%d -- 3x transformed %d,%d -- raw %d %d\n",spamci1,spamci2,22+spamci1,spamci2,spampj1,spampj2);
//
// 		spampj1 = box->patchList[spamci1+23][spamci2+24];
// 		spampj2 = box->patchList[spamci2+23][spamci1+24];
// 		printf("full diag %d %d %d %d\n",spamci1,spamci2,spampj1,spampj2);
// 		box->eps_ij[spampj1][spampj2] = spameps1;
// 		box->eps_ij[spampj2][spampj1] = spameps1;
// 		printf("%d,%d -- 4x transformed %d,%d -- raw %d %d\n",spamci1,spamci2,22+spamci1,22+spamci2,spampj1,spampj2);
		npairs++;


		if(spameps1 > maxeps)
			maxeps=spameps1;
		if(spameps1 < mineps && spameps1 > 0)
			mineps=spameps1;

		spameps1=0.0;
	}

	// rescale interactions
	for(i=0;i<p->nChrom-1;i++)
	{
		for(j=i+1;j<p->nChrom;j++)
		{
			if(box->eps_ij[i][j]>0.0)
			{
// 				printf("old eps ij for %d %d -- %f to %f\n",i,j,box->eps_ij[i][j],10.0*(box->eps_ij[i][j] - mineps)/(maxeps-mineps));
				box->eps_ij[i][j] = p->epsij_rescale*(box->eps_ij[i][j] - mineps)/(maxeps-mineps);
				box->eps_ij[j][i] = box->eps_ij[i][j];
			}
		}
	}


	return npairs;
}

void placePatchesPropPot(struct PARAINPUT *p,struct SYSTEM *box)
{
	// subroutine to place patches such that pairwise intermingling is consistent.
	printf("skipping this!! Stay tuned for updates\n");
}

int placePatchesVectorised(struct PARAINPUT *p,struct SYSTEM *box)
{
	// subroutine to place patches on the surface of the chromosome, along the vector connecting the two centroids

	int i,j,ii,jj,homo_i,homo_j,k,imin,jmin,l,pi,pj;
	double x1,y1,z1;
	double rX,rY,rZ,r,dmin;
	double curr[3];
	double norm1,norm2;
	int d[4];
	int globCount=0;
	// for each pair, check not homologue, read positions, find vector, use half norm to determine positions wrt centre
	int donepatches=0;

	for(i=0;i<p->nChrom;i++)
	{
		homo_i = i+p->nChrom; // homologue id
		for(j=0;j<p->nChrom;j++)
		{
			homo_j = j+p->nChrom;

			if(box->chromType[i] == box->chromType[j] || j <= i) // should not place patches if same type
				continue;

			// if heterologous chromosomes -- find all pairwise distances, identify minimum
			dmin = 2.0*p->ella; // max possible

			ii = box->chromList[i];
			jj = box->chromList[j];
			rX = box->X[jj] - box->X[ii];
			rY = box->Y[jj] - box->Y[ii];
			rZ = box->Z[jj] - box->Z[ii];

			r = sqr(rX) + sqr(rY) + sqr(rZ);
			if(fabs(sqrt(r) - (box->radius[ii] + box->radius[jj]))  <= dmin) // find the pair whose separation is closest to sigma_a + sigma_b
			{
				dmin = fabs(sqrt(r) - (box->radius[ii] + box->radius[jj]));
				imin = i;
				jmin = j;
				for(k=0;k<4;k++) // reset flags
					d[k]=0;
				d[0]=1; // flag the correct pair
			}

			ii = box->chromList[homo_i];
			jj = box->chromList[j];
			rX = box->X[jj] - box->X[ii];
			rY = box->Y[jj] - box->Y[ii];
			rZ = box->Z[jj] - box->Z[ii];

			r = sqr(rX) + sqr(rY) + sqr(rZ);
			if(fabs(sqrt(r) - (box->radius[ii] + box->radius[jj])) <= dmin) // find the pair whose separation is closest to sigma_a + sigma_b
			{
				dmin = fabs(sqrt(r) - (box->radius[ii] + box->radius[jj]));
				imin = homo_i;
				jmin = j;
				for(k=0;k<4;k++) // reset flags only if current distance less than previous min
					d[k]=0;
				d[1]=1; // flag the correct pair
			}

			ii = box->chromList[homo_i];
			jj = box->chromList[homo_j];
			rX = box->X[jj] - box->X[ii];
			rY = box->Y[jj] - box->Y[ii];
			rZ = box->Z[jj] - box->Z[ii];

			r = sqr(rX) + sqr(rY) + sqr(rZ);
			if(fabs(sqrt(r) - (box->radius[ii] + box->radius[jj])) <= dmin) // find the pair whose separation is closest to sigma_a + sigma_b
			{
				dmin = fabs(sqrt(r) - (box->radius[ii] + box->radius[jj]));
				imin = homo_i;
				jmin = homo_j;
				for(k=0;k<4;k++) // reset flags only if current distance less than previous min
					d[k]=0;
				d[2]=1; // flag the correct pair
			}

			ii = box->chromList[i];
			jj = box->chromList[homo_j];
			rX = box->X[jj] - box->X[ii];
			rY = box->Y[jj] - box->Y[ii];
			rZ = box->Z[jj] - box->Z[ii];

			r = sqr(rX) + sqr(rY) + sqr(rZ);
			if(fabs(sqrt(r) - (box->radius[ii] + box->radius[jj])) <= dmin) // find the pair whose separation is closest to sigma_a + sigma_b
			{
				dmin = fabs(sqrt(r) - (box->radius[ii] + box->radius[jj]));
				imin = i;
				jmin = homo_j;
				for(k=0;k<4;k++) // reset flags
					d[k]=0;
				d[3]=1; // flag the correct pair
			}

			// find vector connecting our COMs for the flagged pair and place -- other pairs place randomly on surface

			for(k=0;k<4;k++)
			{
				if(d[k] == 1)
				{
					ii = box->chromList[imin];
					jj = box->chromList[jmin];

					rX = box->X[jj] - box->X[ii];
					rY = box->Y[jj] - box->Y[ii];
					rZ = box->Z[jj] - box->Z[ii];

					norm1 = sqrt(sqr(rX) + sqr(rY)  + sqr(rZ));

					rX /= norm1;
					rY /= norm1;
					rZ /= norm1;

					// place patch at X[ii] + rX, Y[ii] + rY, Z[ii] + rZ -- id which patch
					for(l=1;l<=box->patchPointer[imin][0];l++)
					{

						if(box->patchPointer[imin][l] == box->chromType[jmin] && box->lamin[imin][l] != 1) // patch l on imin points to jmin and is not lamin
						{
							pj = box->patchList[imin][l]; // global id of patch
							box->X[pj] = box->X[ii] + rX*box->radius[ii];
							box->Y[pj] = box->Y[ii] + rY*box->radius[ii];
							box->Z[pj] = box->Z[ii] + rZ*box->radius[ii];

							box->patchposx[imin][l] = rX;
							box->patchposy[imin][l] = rY;
							box->patchposz[imin][l] = rZ;
							box->patchDone[imin][l] = 1;
							donepatches++;

							// put active forces
							box->actox[imin][l] = gaussVar(0,1);
							box->actoy[imin][l] = gaussVar(0,1);
							box->actoz[imin][l] = gaussVar(0,1);
							norm2 = sqrt(box->actox[imin][l]*box->actox[imin][l] + box->actoy[imin][l]*box->actoy[imin][l] + box->actoz[imin][l]*box->actoz[imin][l]);
							box->actox[imin][l] /= norm2;
							box->actoy[imin][l] /= norm2;
							box->actoz[imin][l] /= norm2;

						}

					}

					// place patch at X[jj] - rX, Y[jj] - rY, Z[jj] - rZ -- id which patch
					for(l=1;l<=box->patchPointer[jmin][0];l++)
					{
						if(box->patchPointer[jmin][l] == box->chromType[imin] && box->lamin[jmin][l] != 1) // patch l on jmin points to imin and is not lamin
						{
							pi = box->patchList[jmin][l]; // global id of patch
							box->X[pi] = box->X[jj] - rX*box->radius[jj];
							box->Y[pi] = box->Y[jj] - rY*box->radius[jj];
							box->Z[pi] = box->Z[jj] - rZ*box->radius[jj];

							box->patchposx[jmin][l] = -rX;
							box->patchposy[jmin][l] = -rY;
							box->patchposz[jmin][l] = -rZ;
							box->patchDone[jmin][l] = 1;
							donepatches++;

							// put active forces
							box->actox[jmin][l] = gaussVar(0,1);
							box->actoy[jmin][l] = gaussVar(0,1);
							box->actoz[jmin][l] = gaussVar(0,1);
							norm2 = sqrt(box->actox[jmin][l]*box->actox[jmin][l] + box->actoy[jmin][l]*box->actoy[jmin][l] + box->actoz[jmin][l]*box->actoz[jmin][l]);
							box->actox[jmin][l] /= norm2;
							box->actoy[jmin][l] /= norm2;
							box->actoz[jmin][l] /= norm2;
						}
					}


				}


			}

		}
	}
// 	printf("done placing vectorised %d -- should be %d\n",donepatches,p->nChrom*(p->nChrom-1));
	return donepatches; // should be 2xp->nChrom
}


int placePatchesRandom(struct PARAINPUT *p,struct SYSTEM *box)
{
	// subroutine to place patches randomly with min separation radius patch
	int i,ii,j,jj,k;
	double curr[3];
	int atomnotdone,atomdone;
	double rX,rY,rZ,r,spam;
	double norm1,norm2;
	int globCount=0;
	int donepatches=0;
	for(ii=0;ii<p->nDChrom;ii++)
	{
		i=box->chromList[ii];
// 		printf("placing random %d %d \n",ii,box->patchList[ii][0]);
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
// 			printf("chrom %d (g %d) patch %d (g %d) out of %d -- rad %f -- dist %f -- chrom rad %f\n",ii,i,jj,box->patchList[ii][jj],box->patchList[ii][0],box->radius[box->patchList[ii][jj]],box->radius[box->patchList[ii][jj]]/box->radius[i],box->radius[i]);
			atomdone=0;
// 			j = box->patchList[ii][jj];
			if(box->patchDone[ii][jj]==1 || box->lamin[ii][jj] == 1) // skip if placed vectorised or lamin patch
				continue;
			while(atomdone==0)
			{
				curr[0] = gaussVar(0,1);
				curr[1] = gaussVar(0,1);
				curr[2] = gaussVar(0,1);
				norm1 = sqrt(curr[0]*curr[0] + curr[1]*curr[1] + curr[2]*curr[2]);
				curr[0] /= norm1;
				curr[1] /= norm1;
				curr[2] /= norm1;


				atomnotdone=0;
				for(k=1;k<jj;k++)
				{
					rX=(curr[0]-box->patchposx[ii][k]);
					rY=(curr[1]-box->patchposy[ii][k]);
					rZ=(curr[2]-box->patchposz[ii][k]);

					r=rX*rX+rY*rY+rZ*rZ;
					// here patch placement has to have no consistency condition wrt interaction so just distance
					if(sqrt(r) < box->radius[box->patchList[ii][jj]]/box->radius[i])
						atomnotdone=1;
				}


				atomdone=1-atomnotdone;



			}
			box->patchposx[ii][jj]=curr[0];
			box->patchposy[ii][jj]=curr[1];
			box->patchposz[ii][jj]=curr[2];

			globCount = box->patchList[ii][jj];
			box->X[globCount]=box->X[i] + box->radius[i]*curr[0];
			box->Y[globCount]=box->Y[i] + box->radius[i]*curr[1];
			box->Z[globCount]=box->Z[i] + box->radius[i]*curr[2];
			donepatches++;
// 			printf("chrom %d raw %d patch %d raw %d\n",ii,i,jj,globCount);
// 			printf("chrom %d patch %d -- raw %d  corrob  patchList %d whichPatch %d whichChrom %d\n",ii,jj,globCount,box->patchList[ii][jj],box->whichPatch[globCount],box->whichChrom[globCount]);
			globCount++;

			// put active forces
			box->actox[ii][jj] = gaussVar(0,1);
			box->actoy[ii][jj] = gaussVar(0,1);
			box->actoz[ii][jj] = gaussVar(0,1);
			norm2 = sqrt(box->actox[ii][jj]*box->actox[ii][jj] + box->actoy[ii][jj]*box->actoy[ii][jj] + box->actoz[ii][jj]*box->actoz[ii][jj]);
			box->actox[ii][jj] /= norm2;
			box->actoy[ii][jj] /= norm2;
			box->actoz[ii][jj] /= norm2;
		}
	}
	printf("done placing random %d -- should be %d\n",donepatches,p->nChrom*(p->nChrom-1));
	return donepatches;
}

int placeLaminPatches(struct PARAINPUT *p,struct SYSTEM *box)
{
	int i,ii,j,jj,k,jmin;
	double curr[3];
	int atomnotdone,atomdone;
	double rX,rY,rZ,r,spam;
	double norm1,norm2;
	int globCount=0;
	int donepatches=0;
	double dmin=0.0;
	for(ii=0;ii<p->nDChrom;ii++)
	{
		i = box->chromList[ii];
		jj = box->patchList[ii][0]; // the last patch is lamin
		globCount = box->patchList[ii][jj]; // id of lamin patch
		if(box->lamin[ii][jj] != 1)
		{
			printf("something wrong with lamin patch assignment, please check for chrom %d\n",ii);
			exit(0);
		}

		if(box->patchDone[ii][jj] ==1)
		{
			printf("lamin patch already placed, something wrong, %d (raw %d) - %f %f %f\n",jj,globCount,box->X[globCount],box->Y[globCount],box->Z[globCount]);
			exit(0);
		}
		// check closest surface point to chrom COM
		dmin = 2.0*p->ella;
		for(j=0;j<p->nSurfPoints;j++)
		{
			rX = box->X[i] - box->surfPointx[j];
			rY = box->Y[i] - box->surfPointy[j];
			rZ = box->Z[i] - box->surfPointz[j];

			r = sqrt(sqr(rX) + sqr(rY) + sqr(rZ));

			if(r < dmin)
			{
				dmin=r;
				jmin=j;
			}
		}

		// for the min sep, find the vector

		rX = box->surfPointx[jmin] - box->X[i];
		rY = box->surfPointy[jmin] - box->Y[i];
		rZ = box->surfPointz[jmin] - box->Z[i];

		norm1 = sqrt(sqr(rX) + sqr(rY)  + sqr(rZ));

		rX /= norm1;
		rY /= norm1;
		rZ /= norm1;

		box->X[globCount] = box->X[i] + rX*box->radius[i];
		box->Y[globCount] = box->Y[i] + rY*box->radius[i];
		box->Z[globCount] = box->Z[i] + rZ*box->radius[i];

		box->patchposx[ii][jj] = rX;
		box->patchposy[ii][jj] = rY;
		box->patchposz[ii][jj] = rZ;
		box->patchDone[ii][jj] = 1;

// 		printf("chrom %d (dist bound %f) has lamin patch %d (raw %d) placed at %f %f %f wrt %f %f %f (glob %f %f %f)\n",ii,dmin,jj,globCount,rX,rY,rZ,box->X[i],box->Y[i],box->Z[i],box->X[globCount],box->Y[globCount],box->Z[globCount]);

		donepatches++;

		// put active forces
		box->actox[ii][jj] = gaussVar(0,1);
		box->actoy[ii][jj] = gaussVar(0,1);
		box->actoz[ii][jj] = gaussVar(0,1);
		norm2 = sqrt(box->actox[ii][jj]*box->actox[ii][jj] + box->actoy[ii][jj]*box->actoy[ii][jj] + box->actoz[ii][jj]*box->actoz[ii][jj]);
		box->actox[ii][jj] /= norm2;
		box->actoy[ii][jj] /= norm2;
		box->actoz[ii][jj] /= norm2;


		// specify orientation vector of the chromosome as that of the body fixed lamin patch positions
		box->ox[i] = rX;
		box->oy[i] = rY;
		box->oz[i] = rZ;

	}
	printf("placed lamin patches %d -- should be %d\n",donepatches,p->nDChrom);
	return donepatches;
}
