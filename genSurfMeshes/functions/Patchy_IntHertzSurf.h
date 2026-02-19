double calcPatchAnchor(struct SYSTEM *box,struct PARAINPUT *p)
{
	int ii,n1,j,i,neipar,ch;
	double spamx,spamy,spamz,leverx,levery,leverz,diffX,diffY,diffZ;
	
	double rij,rijsq,invrij,r0,ff,patchcut;
	double en=0.0;
	// get patch force interactions, calculate com force on chromosomes, torque on chromosomes
	for(ii=0;ii<p->nChrom;ii++)
	{
		n1 = box->chromList[ii]; // raw index of chromosome
		for(j=1;j<=box->patchList[ii][0];j++) // loop through patches
		{
			i = box->patchList[ii][j]; // which patch
			spamx = box->Gpatchposx[ii][j]; // connected to (in global coordinates)?!
			spamy = box->Gpatchposy[ii][j];
			spamz = box->Gpatchposz[ii][j];
			
			// 
// 			printf("cal sep for %d %f %f %f -- %f %f %f -- diff %f %f %f\n",i,box->X[i],box->Y[i],box->Z[i],spamx,spamy,spamz,box->X[i]-spamx,box->Y[i]-spamy,box->Z[i]-spamz);
			diffX=box->X[i]-spamx;
// 			diffX-=box->lx*lround(diffX/box->lx); // this line applies PBC
			
			diffY=box->Y[i]-spamy;
// 			diffY-=box->ly*lround(diffY/box->ly);
			
			diffZ=box->Z[i]-spamz;
// 			diffZ-=box->lz*lround(diffZ/box->lz);
		
			rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
			if(rijsq > 0.0)
			{
				rij = sqrt(rijsq);
				invrij = 1.0/rij;
			}
			else
			{
				rij = 0.0;
				invrij = 0.0;
			}
// 			rij=sqrt(rijsq);
// 			invrij = 1.0/rij;
			/*ch = box->type[i] + box->type[n1];
			r0 = 2.0*p->sigma[ch]; // check diam/ radius scenes here	*/	
			
			r0=3.0*(box->radius[n1] + box->radius[i]);
			
			// spring force COM on both
// 			en = en + -0.5*box->isBond[i]*p->k_spring*r0*r0*log(1.0 - pow(rij/r0,2.0));
// 			box->en[i]+=-0.5*0.5*p->k_spring*r0*r0*log(1.0 - pow(rij/r0,2.0));
// 			box->en[n1]+=-0.5*0.5*p->k_spring*r0*r0*log(1.0 - pow(rij/r0,2.0));
// 			printf("force ingredients ch %d par %d, patch %d -- sp %f bond %f r0 %f rij %f ff %f\n",ii,n1,i,box->k_spring[i],box->isBond[i],r0,rij,-box->k_spring[i]*box->isBond[i]*rij*r0*r0/(r0*r0 - rijsq));
//
			
// 			en += box->k_spring[i]*box->isBond[i]*rijsq;
//
// 			ff = -box->k_spring[i]*box->isBond[i]*rij;
			if(rijsq > r0*r0)
			{
				printf("WHOA WHOA c %d %d -- p %d %d %f %f\n",ii,n1,j,i,rijsq,r0*r0);
				exit(0);
			}

			ff = -box->k_spring[i]*box->isBond[i]*rij*r0*r0/(r0*r0 - rijsq); // NOTE: be careful here, only patch-types have non-zero k_spring - chromosomes have 0
			en += -0.5*box->k_spring[i]*box->isBond[i]*r0*r0*log(1.0 - rijsq/(r0*r0));

			
			box->fx[i] += ff*diffX*invrij;
			box->fy[i] += ff*diffY*invrij;
			box->fz[i] += ff*diffZ*invrij;
			box->fx[n1] += -ff*diffX*invrij;
			box->fy[n1] += -ff*diffY*invrij;
			box->fz[n1] += -ff*diffZ*invrij;
// 			printf("%d ch %d p -- spring force %f %f %f towards %f %f %f\n",n1,i,-ff*diffX*invrij,-ff*diffY*invrij,-ff*diffZ*invrij,diffX,diffY,diffZ);
			// torque on chrom -- use force and patchpos in global coords
			
			leverx = spamx - box->X[n1];
			levery = spamy - box->Y[n1];
			leverz = spamz - box->Z[n1];
// 			printf("par %d sep %f -- ff %f and force %f %f %f\n",i,rij,ff,box->fx[i],box->fy[i],box->fz[i]);
			box->torquex[n1] += -1.0*(levery*ff*diffZ*invrij - leverz*ff*diffY*invrij); // NOTE: minus sign important here because force on COM is neg of ff dependent
			box->torquey[n1] += -1.0*(leverz*ff*diffX*invrij - leverx*ff*diffZ*invrij);
			box->torquez[n1] += -1.0*(leverx*ff*diffY*invrij - levery*ff*diffX*invrij);
// 			printf("resultant torque %f %f %f \n",box->torquex[n1],box->torquey[n1],box->torquez[n1]);
// 			(v1).x = (v2).y * (v3).z - (v2).z * (v3).y,
// 			(v1).y = (v2).z * (v3).x - (v2).x * (v3).z,
// 			(v1).z = (v2).x * (v3).y - (v2).y * (v3).x

			
		}
	}
// 	en=0;
	if(isnan(en) || isnan(-en))
		printf(" patch anchor nan energy\n");
	return en;
}


double calcPatchPatch(struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,j,ii,n1,jj,neipar,ch;
	int ci,cj,pi,pj;
	double en=0.0;
    
	double spamx,spamy,spamz,leverx,levery,leverz,diffX,diffY,diffZ;
	
	double rij,rijsq,invrij,r0,ff,rijanch;
	
	double patchcut,alphij,invwidthwell,tanhfact;

	// get patch - patch forces -- repulsive blob + attractive Gaussian thing
	
	for(i=0;i<p->nAto;i++)
	{
// 		if(box->type[i] ==0)
// 		{
// 			constant force, divide by gamma, unit vector
// 			printf("patch active forces not implemented\n");
// 		}
// 		if(box->type[i]==2)
// 		{
// 			if wall patch check wall separation
// 			printf("patch wall interactions not implemented\n");
// 		}
		for(j=0;j<p->nAto;j++)
		{
			// check if patch-patch
			
			// introduce subsequent checks -- same chromosome or not
			// introduce subsequent features -- interaction matrix
			neipar = j;
			ff=0.0;
			if(box->type[i]==0 && box->type[neipar]==0 && neipar > i)
			{
			// if patch, compute separation from patchpos (check coordinate reference)
			// calculate spring along straight vector connecting them
			
			// patch-patch interactions
// 				ch=box->type[i]+box->type[neipar];
				diffX=box->X[neipar]-box->X[i];
// 				diffX-=box->lx*lround(diffX/box->lx); // this line applies PBC
				
				diffY=box->Y[neipar]-box->Y[i];
// 				diffY-=box->ly*lround(diffY/box->ly);
				
				diffZ=box->Z[neipar]-box->Z[i];
// 				diffZ-=box->lz*lround(diffZ/box->lz);
		
				rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
				
				patchcut = (box->radius[i] + box->radius[neipar]);



				
// 				if(rijsq < 9.0*patchcut*patchcut && box->eps_ij[i][neipar] !=0.0)
				if(rijsq < 4.0*patchcut*patchcut)
				{
// 					rij = sqrt(rijsq);
// 					invrij = 1.0/rij;
					if(rijsq > 0.0)
					{
						rij = sqrt(rijsq);
						invrij = 1.0/rij;
					}
					else
					{
						rij = 0.0;
						invrij = 0.0;
					}

					
					// blob repulsion
					alphij = 1.0/(1.04*patchcut);
					ff = 2.0*p->eps_blob*alphij*rij*exp(-alphij*rijsq);
					

					
					// attraction --- tanh function derivative
					invwidthwell = 1.0/(0.25*patchcut); // inverse of 0.1*patchcut
					tanhfact = tanh((rij - 1.2*patchcut)*invwidthwell);
					ff += -box->eps_ij[i][neipar]*(1.0 - tanhfact*tanhfact)*invwidthwell;
					
					en += p->eps_blob*exp(-alphij*rijsq) + 0.5*box->eps_ij[i][neipar]*tanhfact - 0.5*box->eps_ij[i][neipar];

// 					box->IMEne += p->eps_blob*exp(-alphij*rijsq) + 0.5*box->eps_ij[i][neipar]*tanh((rij - 1.05*patchcut)*invwidthwell) - 0.5*box->eps_ij[i][neipar];
					box->IMEne += 0.5*box->eps_ij[i][neipar]*tanhfact - 0.5*box->eps_ij[i][neipar];


					if(box->eps_ij[i][neipar] > 0.0)
					{
						ci = box->whichChrom[i];
						cj = box->whichChrom[neipar];
						pi = box->whichPatch[i];
						pj = box->whichPatch[neipar];

						spamx = box->Gpatchposx[ci][pi] - box->Gpatchposx[cj][pj];
						spamy = box->Gpatchposy[ci][pi] - box->Gpatchposy[cj][pj];
						spamz = box->Gpatchposz[ci][pi] - box->Gpatchposz[cj][pj];

						rijanch = spamx*spamx + spamy*spamy + spamz*spamz;
						if(rijanch > 0.0)
							rijanch = sqrt(rijanch);

// 						box->isIntermingling[i] = 1.0 - MIN(1.0,rijanch/(4.0*patchcut));
// 						box->isIntermingling[neipar] = 1.0 - MIN(1.0,rijanch/(4.0*patchcut));
						if(box->checkIntermingling[i] ==0.0)
						{
						box->isIntermingling[i] = 0.5 - 0.5*tanh((rijanch-2.0*patchcut)/(0.5*patchcut));
						box->checkIntermingling[i] = 1.0;
						}
						if(box->checkIntermingling[neipar] == 0.0)
						{
						box->isIntermingling[neipar] = 0.5 - 0.5*tanh((rijanch-2.0*patchcut)/(0.5*patchcut));
						box->checkIntermingling[neipar] = 1.0;
						}

					}
					

// 					if(box->eps_ij[i][neipar] > 0.0)
// 					printf("patch raw %d on chrom %d patch %d with patch %d on chrom %d patch %d -- %f %f %f %f \n",i,box->whichChrom[i],box->whichPatch[i],neipar,box->whichChrom[neipar],box->whichPatch[neipar],sqrt(rijsq),box->eps_ij[i][neipar],p->eps_blob*exp(-alphij*rijsq),0.5*box->eps_ij[i][neipar]*tanh((rij - 1.05*patchcut)*invwidthwell) - 0.5*box->eps_ij[i][neipar]);
// 					ff += -box->eps_ij[i][neipar]*( -12.0*pow(0.15*patchcut/rij,12.0)/rij + 6.0*pow(0.15*patchcut/rij,6.0)/rij);
					
					
// 					ff = 2.0*p->eps_excl[ch]*p->alpha_excl[ch]*rij*exp(-p->alpha_excl[ch]*rijsq);
// 					ff += +2.0*p->eps_hc_in*rij*(1.0 + rij*p->alpha_hc*(rmin -rij))*exp(-p->alpha_hc*pow(rmin - rij,2.0));
// 					ff = -box->eps_ij[i][neipar]*rij;
					
// 					ff=0.0;
					
					box->fx[i] += -ff*diffX*invrij;
					box->fy[i] += -ff*diffY*invrij;
					box->fz[i] += -ff*diffZ*invrij;
					
					box->fx[neipar] += ff*diffX*invrij;
					box->fy[neipar] += ff*diffY*invrij;
					box->fz[neipar] += ff*diffZ*invrij;
// 					printf("forces on par %d %d at sep %f w eps %f is %f -- disp %f %f %f\n",i,neipar,rij,box->eps_ij[i][neipar],ff,diffX,diffY,diffZ);
// 					printf("forces reso par %d -- %f %f %f and %d %f %f %f\n",i,-ff*diffX*invrij,-ff*diffY*invrij,-ff*diffZ*invrij,neipar,ff*diffX*invrij,ff*diffY*invrij,ff*diffZ*invrij);
				}
			}
		
			
		}
	}
// 	en = 0.0;
	if(isnan(en) || isnan(-en))
		printf(" patch patch nan energy\n");
	return en;
}

double calcChromInt(struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,j,ii,n1,jj,neipar,ch;
	double en=0.0;
    
	double spamx,spamy,spamz,leverx,levery,leverz,diffX,diffY,diffZ;
	
	double rij,rijsq,invrij,r0,ff,patchcut,alphij;
	
	// NOTE: Calculate chrom-chrom Hertz repulsion
	for(ii=0;ii<p->nChrom-1;ii++)
	{
		i=box->chromList[ii];
		for(jj=ii+1;jj<p->nChrom;jj++)
		{
			j = box->chromList[jj];
			diffX=box->X[j]-box->X[i];
			
			diffY=box->Y[j]-box->Y[i];
			
			diffZ=box->Z[j]-box->Z[i];
	
			rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;
			
			
			
			patchcut = 1.2*(box->radius[i] + box->radius[j]);
			// find overlap


			// apply Hertz force to i and j
			if(rijsq <= patchcut*patchcut) // overlapping
			{
				if(rijsq > 0.0)
				{
					rij = sqrt(rijsq);
					invrij = 1.0/rij;
				}
				else
				{
					rij = 0.0;
					invrij = 0.0;
				}

				ff = p->eps_hertz*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
				
				en += p->eps_hertz*pow((1-rij/patchcut),2.5);
// 				printf("chroms %d raw %d, rad %f and %d raw %d, rad %f separated by %f -- repulsion %f\n",ii,i,box->radius[i],jj,j,box->radius[j],rij,p->eps_hertz*pow((1-rij/patchcut),2.5));

					
				box->fx[i] += -ff*diffX*invrij;
				box->fy[i] += -ff*diffY*invrij;
				box->fz[i] += -ff*diffZ*invrij;
				
				box->fx[j] += ff*diffX*invrij;
				box->fy[j] += ff*diffY*invrij;
				box->fz[j] += ff*diffZ*invrij;
			}
			
		}
		
	}
	
	// NOTE: Calculate Chrom Patch Gaussian core repulsion
	for(ii=0;ii<p->nChrom;ii++)
	{
		i = box->chromList[ii];
		for(j=0;j<p->nAto;j++)
		{
			// make sure j is not self
			if(j==i)
				continue;

			ff=0.0;
			// make sure j is not another chrom
			if(box->type[j]==0)
			{
				// calculate a soft core interaction here
// 				continue;
				neipar=j;
				diffX=box->X[neipar]-box->X[i];
// 				diffX-=box->lx*lround(diffX/box->lx); // this line applies PBC

				diffY=box->Y[neipar]-box->Y[i];
// 				diffY-=box->ly*lround(diffY/box->ly);

				diffZ=box->Z[neipar]-box->Z[i];
// 				diffZ-=box->lz*lround(diffZ/box->lz);

				rijsq=diffX*diffX+diffY*diffY+diffZ*diffZ;

				patchcut = (box->radius[i] + box->radius[neipar]);

// 				if(rijsq < 9.0*patchcut*patchcut && box->eps_ij[i][neipar] !=0.0)
				if(rijsq < 4.0*patchcut*patchcut)
				{
					if(rijsq > 0.0)
					{
						rij = sqrt(rijsq);
						invrij = 1.0/rij;
					}
					else
					{
						rij = 0.0;
						invrij = 0.0;
					}


					// blob repulsion
					alphij = 1.0/(1.04*patchcut);
					ff = 2.0*p->eps_blob*alphij*rij*exp(-alphij*rijsq);

					en += p->eps_blob*exp(-alphij*rijsq);

// 					printf("chrom %d raw %d rad %f and patch %d on chrom %d as %d rad %f -- sep %f repuls %f\n",ii,i,box->radius[i],j,box->whichChrom[j],box->whichPatch[j],box->radius[j],rij,p->eps_blob*exp(-alphij*rijsq));



					// update forces on i and j



					box->fx[i] += -ff*diffX*invrij;
					box->fy[i] += -ff*diffY*invrij;
					box->fz[i] += -ff*diffZ*invrij;

					box->fx[neipar] += ff*diffX*invrij;
					box->fy[neipar] += ff*diffY*invrij;
					box->fz[neipar] += ff*diffZ*invrij;
				}

			}

		}
	}
	
// 	en=0.0;
	if(isnan(en) || isnan(-en))
		printf(" chrom-chrom patch nan energy\n");
	return en;
	
}

double calcWallInt(struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,j,ii,n1,jj,neipar,ch,kk,k;
	double en=0.0;
    
	double diffX,diffY,diffZ;
	double spamx,spamy,spamz,fundx,fundy,fundz;
	double rij,rijsq,invrij,r0,ff,rijanch;
	double patchcut,alphij,invwidthwell,tanhfact,lamintspam;

	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5);
	for(i=0;i<p->nSurfPoints;i++)
	{
		box->surffx[i]=0.0;
		box->surffy[i]=0.0;
		box->surffz[i]=0.0;
	}
		
	
	for(i=0;i<p->nSurfPoints;i++)
	{
		// find self interaction
		spamx = box->surfPointx[i];
        spamy = box->surfPointy[i];
        spamz = box->surfPointz[i];

        fundx = box->surfFundx[i];
        fundy = box->surfFundy[i];
        fundz = box->surfFundz[i];
// 		printf("surf %f %f %f -- fund %f %f %f\n",spamx,spamy,spamz,fundx,fundy,fundz);

        rijsq = (spamx-fundx)*(spamx-fundx) + (spamy-fundy)*(spamy-fundy) + (spamz-fundz)*(spamz-fundz);
		rij=0.0;
		invrij = 0.0;
        if(rijsq > 0)
		{
			rij=sqrt(rijsq);
			invrij = 1.0/rij;
		}
		
		ff = -p->surfSpringSelf*(p->surftauinv)*rij;
		en += 0.5*p->surfSpringSelf*(p->surftauinv)*rijsq;

// 		ff = -p->surfSpringSelf*(p->surftauinv)*rij;
// 		en += 0.5*p->surfSpringSelf*(p->surftauinv)*rijsq;

// 		ff = 0.0;
// 		r0 = 3.0*avesep;
// 		if(rijsq >= r0*r0)
// 		{
// 			printf("GRAVE PATHOLOGY, exiting...\n");
// 			exit(0);
// 		}
// 		ff = -37.75*p->surfSpringSelf*rij*r0*r0/(r0*r0 - rijsq); // entropic spring for surface particles
// 		en += -0.5*37.75*p->surfSpringSelf*r0*r0*log(1.0 - rijsq/(r0*r0));

// 		printf("sf %d -- sep %f %f %f -- rij %f -- ff %f -- components %f %f %f\n",i,spamx-fundx,spamy-fundy,spamz-fundz,rij,ff,ff*(spamx-fundx)*invrij,ff*(spamy-fundy)*invrij,ff*(spamz-fundz)*invrij);
		box->surffx[i] += ff*(spamx-fundx)*invrij;
		box->surffy[i] += ff*(spamy-fundy)*invrij;
		box->surffz[i] += ff*(spamz-fundz)*invrij;
		// find neighbours of i for cross interaction
		for(jj=1;jj<=box->surfNeigh[i][0];jj++)
		{
			// avoid double counting
			j = box->surfNeigh[i][jj];

			patchcut = 1.2*avesep;


			if(j>i)
			{
				// find cross interaction
				
				diffX = box->surfPointx[j] - box->surfPointx[i];
				diffY = box->surfPointy[j] - box->surfPointy[i];
				diffZ = box->surfPointz[j] - box->surfPointz[i];
				
				rijsq = diffX*diffX + diffY*diffY + diffZ*diffZ;
				if(rijsq > 0)
				{
					rij=sqrt(rijsq);
					invrij = 1.0/rij;
				}
				else
				{
					rij = 0.0;
					invrij = 0.0;
				}

// 				printf("%d is %d neighbour (out of %d) of surf par %d -- d %f sep %f \n",j,jj,box->surfNeigh[i][0],i,rij,patchcut);
				




				ff=0.0;
				if(rij < patchcut)
				{
// 					ff = p->surfSpringCross*(p->surftauinv/p->surfDiffu)*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
// 					en += p->surfSpringCross*(p->surftauinv/p->surfDiffu)*pow((1-rij/patchcut),2.5);

					ff = p->surfSpringCross*(p->surftauinv)*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
					en += p->surfSpringCross*(p->surftauinv)*pow((1-rij/patchcut),2.5);



				}

// 				en += 37.75*0.5*p->surfSpringCross*rijsq;
// 				ff += -37.75*p->surfSpringCross*rij;
// 				en += 0.5*p->surfSpringCross*(p->surftauinv/p->surfDiffu)*rijsq;
// 				ff += -p->surfSpringCross*(p->surftauinv/p->surfDiffu)*rij;

				en += 0.5*p->surfSpringCross*(p->surftauinv)*rijsq;
				ff += -p->surfSpringCross*(p->surftauinv)*rij;


// 				r0 = 1.5*avesep;
// 				ff += -37.75*p->surfSpringCross*rij*r0*r0/(r0*r0 - rijsq); // entropic spring for surface particles
// 				en += -0.5*37.75*p->surfSpringCross*r0*r0*log(1.0 - rijsq/(r0*r0));

// 				spamx = box->surfFundx[j] - box->surfFundx[i];
// 				spamy = box->surfFundy[j] - box->surfFundy[i];
// 				spamz = box->surfFundz[j] - box->surfFundz[i];
//
// 				rijsq = (diffX - spamx)*(diffX - spamx) + (diffY - spamy)*(diffY - spamy) + (diffZ - spamz)*(diffZ - spamz);
// 				rij = sqrt(rijsq);
// 				invrij = 1.0/rij;
//
// 				en += 37.75*0.5*p->surfSpringCross*rijsq;
//
// 				ff = -37.75*p->surfSpringCross*rij;
//
// 				box->surffx[i] -= ff*(diffX - spamx)*invrij;
// 				box->surffy[i] -= ff*(diffY - spamy)*invrij;
// 				box->surffz[i] -= ff*(diffZ - spamz)*invrij;
//
// 				box->surffx[j] += ff*(diffX - spamx)*invrij;
// 				box->surffy[j] += ff*(diffY - spamy)*invrij;
// 				box->surffz[j] += ff*(diffZ - spamz)*invrij;
//

// 				r0 = 3.0*avesep;
// 				ff += -37.75*p->surfSpringCross*rij*r0*r0/(r0*r0 - rijsq); // entropic spring for surface particles
// 				en += -0.5*37.75*p->surfSpringCross*r0*r0*log(1.0 - rijsq/(r0*r0));


				
				box->surffx[i] -= ff*(diffX)*invrij;
				box->surffy[i] -= ff*(diffY)*invrij;
				box->surffz[i] -= ff*(diffZ)*invrij;

				box->surffx[j] += ff*(diffX)*invrij;
				box->surffy[j] += ff*(diffY)*invrij;
				box->surffz[j] += ff*(diffZ)*invrij;
			
			}
		}
		box->isMembIntermingling[i] = 0.0; // reset intermingling for surface points
		

	}
// 	printf("done with wall-wall interactions -- %f \n",en);
	double minrijlam=4.0*p->ella;
	for(i=0;i<p->nSurfPoints;i++)
	{
		minrijlam=4.0*p->ella;
		// find chrom and patch interaction
		for(jj=0;jj<p->nChrom;jj++)
		{
			// note the radii of the particles
			j = box->chromList[jj];
			diffX = box->X[j] - box->surfPointx[i];
			diffY = box->Y[j] - box->surfPointy[i];
			diffZ = box->Z[j] - box->surfPointz[i];
			rijsq = diffX*diffX + diffY*diffY + diffZ*diffZ;



			patchcut = 0.5*avesep + box->radius[j];

// 			patchcut = pow(2.0,1.0/6.0)*(0.5*avesep + box->radius[j]);
// 			printf("distance of chrom %d from wall %d is %f %f \n",jj,i,rijsq,patchcut);
			if(rijsq <= sqr(pow(2.0,1.0/6.0)*patchcut)) // overlapping
			{
				rij = sqrt(rijsq);
				invrij = 1.0/rij;

// 				ff = 0.1*p->eps_hertz*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
// 				en += 0.1*p->eps_hertz*pow((1-rij/patchcut),2.5);

				ff = -4.0*(-12.0*pow(patchcut/rij,12.0)*invrij + 6.0*pow(patchcut/rij,6.0)*invrij);

				en += 4.0*(pow(patchcut/rij,12.0) - pow(patchcut/rij,6.0)) + 1.0;

// 				printf("chrom %d raw %d is %f away from boundary at %d, cut-off %f\n",jj,j,rij,i,patchcut);

				box->surffx[i] += -ff*diffX*invrij;
				box->surffy[i] += -ff*diffY*invrij;
				box->surffz[i] += -ff*diffZ*invrij;

				box->fx[j] += ff*diffX*invrij;
				box->fy[j] += ff*diffY*invrij;
				box->fz[j] += ff*diffZ*invrij;

// 				printf("updated forces for wall %d, chrom %d\n",i,j);
			}


// 			// if type patch then blob -- lengthscales from sum of radii etc
			for(kk=1;kk<=box->patchList[jj][0];kk++) // loop through patches
			{
				k = box->patchList[jj][kk];
				diffX = box->X[k] - box->surfPointx[i];
				diffY = box->Y[k] - box->surfPointy[i];
				diffZ = box->Z[k] - box->surfPointz[i];
				rijsq = diffX*diffX + diffY*diffY + diffZ*diffZ;

				patchcut = 0.5*avesep + box->radius[k];
				ff=0.0;

				if(rijsq <= sqr(2.0*patchcut))
				{
					// if within cutoff then WCA repulsion
					if(rijsq <= sqr(pow(2.0,1.0/6.0)*patchcut)) // overlapping
					{
						rij = sqrt(rijsq);
						invrij = 1.0/rij;

	// 					alphij = 1.0/(0.5*patchcut);
	// 					ff = 0.1*2.0*p->eps_blob*alphij*rij*exp(-alphij*rijsq);
	// 					en += 0.1*p->eps_blob*exp(-alphij*rijsq);
	// 					printf("patch %d is on chrom %d (%d) as %d (%d) sep %f -- contri %f\n",k,box->whichChrom[k],jj,box->whichPatch[k],kk,rij, p->eps_blob*exp(-alphij*rijsq));
						if(rijsq < 0.5*patchcut)
							printf("VERY VERY large surface overlap\n");

						ff = -4*(-12.0*pow(patchcut/rij,12.0)*invrij + 6.0*pow(patchcut/rij,6.0)*invrij);
						en += 4.0*(pow(patchcut/rij,12.0) - pow(patchcut/rij,6.0)) + 1.0;
					}


					if(box->lamin[k]>0.0) // lamin patch
					{

						invwidthwell = 1.0/(0.01*box->lamin[k]*patchcut); // inverse of 0.1*patchcut
						tanhfact = tanh((rij - (1+0.01*box->lamin[k])*patchcut)*invwidthwell);
						ff += -box->lamin[k]*(1.0 - tanhfact*tanhfact)*invwidthwell;

						en += 0.5*box->lamin[k]*tanhfact - 0.5*box->lamin[k];
// 						box->IMEne += 0.5*box->lamin[k]*tanhfact - 0.5*box->lamin[k];

						// figure out lamin modifier and add tanh attraction
						spamx = box->Gpatchposx[jj][kk] - box->surfFundx[i];
						spamy = box->Gpatchposy[jj][kk] - box->surfFundy[i];
						spamz = box->Gpatchposz[jj][kk] - box->surfFundz[i];

						rijanch = spamx*spamx + spamy*spamy + spamz*spamz;
						if(rijanch > 0.0)
							rijanch = sqrt(rijanch);


						lamintspam = 0.5 - 0.5*tanh((rijanch-1.4*patchcut)/(0.25*1.4*patchcut));

						box->isIntermingling[k] = MAX(box->isIntermingling[k],lamintspam);


// 							box->isIntermingling[k] = 1.0 - MIN(1.0,rij/(2.0*patchcut));
// 							box->isMembIntermingling[i] = 1.0 - MIN(1.0,rij/(2.0*patchcut));
// 						if(box->checkIntermingling[k] == 0.0)
// 						{
// 							box->isIntermingling[k] = 0.5 - 0.5*tanh((rijanch-1.4*patchcut)/(0.25*1.4*patchcut));
// 							box->checkIntermingling[k] = 1.0;
// 						}
						box->isMembIntermingling[i] = 0.5 - 0.5*tanh((rijanch-1.4*patchcut)/(0.25*1.4*patchcut));
						if(rij < minrijlam)
							minrijlam = rij;


// 						printf("patch %d on chrom %d as %d is lamin, sep %f -- contri %f\n",k,box->whichChrom[k],box->whichPatch[k],rij,0.5*10.0*box->lamin[k]*tanh((rij - 1.05*patchcut)*invwidthwell) - 0.5*10.0*box->lamin[k]);
					}

				}



				box->surffx[i] += -ff*diffX*invrij;
				box->surffy[i] += -ff*diffY*invrij;
				box->surffz[i] += -ff*diffZ*invrij;

				box->fx[k] += ff*diffX*invrij;
				box->fy[k] += ff*diffY*invrij;
				box->fz[k] += ff*diffZ*invrij;


			}
		}
	}
// 	en=0.0;
	if(isnan(en) || isnan(-en))
		printf(" lamin nan energy\n");
	return en;
}


double energyForce(struct SYSTEM *box,struct PARAINPUT *p)
{
	// calculate interactions and forces here
	 
	// the main ones are patch chrom FENE, patch active, patch wall
	
	// then for the chromosomes we calculate net COM force and torque 
	
	int i,j,ii,n1,jj,neipar,ch;
	double en=0.0;
    
	double spamx,spamy,spamz,leverx,levery,leverz,diffX,diffY,diffZ;
	
	double rij,rijsq,invrij,r0,ff,patchcut;
	box->IMEne=0.0;
	// reset forces NOTE: when adding active force be careful of order
	
	for(i=0;i<p->nAto;i++)
	{
		box->fx[i]=0.0;
		box->fy[i]=0.0;
		box->fz[i]=0.0;
		box->torquex[i]=0.0;
		box->torquey[i]=0.0;
		box->torquez[i]=0.0;

		box->isIntermingling[i] = 0.0;
		box->checkIntermingling[i] = 0.0;
	}	
// 	printf("reset forces\n");
	// NOTE: Calculate patch-anchor interactions and compute torque, COM forces
	en = en + calcPatchAnchor(box,p);
	
// 	printf("calculated anchor forces and torques %f\n",en);
	// NOTE: Calculate Patch-Patch attractive + repulsive interactions
	
	en = en + calcPatchPatch(box,p);
// 	printf("calculated patch patch interactions %f\n",en);

	// NOTE: Calculate Chrom-Chrom, Chrom-Patch (chrom-wall calculated in calcWallInt)
	
	en = en + calcChromInt(box,p);
// 	printf("calculated chromosome repulsion %f\n",en);

	// NOTE: Calculate Wall interactions, surfPoint-self, surfPoint-surpoint, surfPoint-lamin patch, surfPoint-Chrom
// 	printf("calculating wall interactions %f\n",en);
	en = en + calcWallInt(box,p);
	

// 	printf("done with wall interactions %f\n",en);
	
	// NOTE: at this point, we have first reset all forces and torques -- then we call the individual modules to compute these interactions
	
	// Now we have the final set of forces and torques before integration. Positions + forces can be used here to calculate virial pressure etc
	// If we need the above subroutines to return additional output, pass as arg (static array/variable) as elsewhere
	
	return en;
}

void kinEnergy(struct SYSTEM *box,struct PARAINPUT *p)
{
	// NOTE: Rotational kinetic energy needs to be implemented
	
    int i;
    box->KE = 0.0; // reset for call on this time step
    for(i=0;i<p->nAto;i++)
    {
        box->KE += 0.5*(box->vx[i]*box->vx[i] + box->vy[i]*box->vy[i] + box->vz[i]*box->vz[i]);
    }

//     box->KE /= p->nAto;
}

