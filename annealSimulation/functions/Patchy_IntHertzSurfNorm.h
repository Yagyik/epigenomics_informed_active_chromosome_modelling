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
				printf("WHOA WHOA c %d %d -- p %d %d %f %f -- rads %f %f\n",ii,n1,j,i,sqrt(rijsq),r0,box->radius[n1],box->radius[i]);
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
	double ath,wth;

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
			// check if patch-pat0.5*ch
			
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
// 					wth = 0.4*box->eps_ij[i][neipar]/p->eps_blob;
					wth = 0.4*p->epsij_scale/p->eps_blob;
					invwidthwell = 1.0/(wth*patchcut); // inverse of 0.1*patchcut

// 					ath = 1.0 + 0.4*box->eps_ij[i][neipar]/p->eps_blob;
					ath = 1.0 + 0.4*p->epsij_scale/p->eps_blob;
					tanhfact = tanh((rij - ath*patchcut)*invwidthwell);

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

// double calcWallInt(struct SYSTEM *box,struct PARAINPUT *p)
// {
//
//
// }


double calcSurfSurf(struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,j,ii,n1,jj,neipar,ch,kk,k;
	double en=0.0;

	double diffX,diffY,diffZ;
	double spamx,spamy,spamz,fundx,fundy,fundz;
	double rij,rijsq,invrij,r0,ff,rijanch;
	double patchcut,alphij,invwidthwell,tanhfact,lamintspam;

	double tmpexp,tmpax,tmpbx;

	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5);
	for(i=0;i<p->nSurfPoints;i++)
	{
		box->surffx[i]=0.0;
		box->surffy[i]=0.0;
		box->surffz[i]=0.0;
	}


	for(i=0;i<p->nSurfPoints;i++)
	{
		// find self interaction

		if(box->surfActive[i]==0)
			continue;
		spamx = box->surfPointx[i];
        spamy = box->surfPointy[i];
        spamz = box->surfPointz[i];

        fundx = box->surfFundx[i];
        fundy = box->surfFundy[i];
        fundz = box->surfFundz[i];

        rijsq = (spamx-fundx)*(spamx-fundx) + (spamy-fundy)*(spamy-fundy) + (spamz-fundz)*(spamz-fundz);
		rij=0.0;
		invrij = 0.0;
        if(rijsq > 0)
		{
			rij=sqrt(rijsq);
			invrij = 1.0/rij;
		}

// 		if(rijsq > 0.05*avesep)
// 		{
// 			rij = sqrt(rijsq);
// 			en += p->surfSpringSelf*exp(1.0/(0.5*avesep - rij))/(2*avesep - rij); // constrained to be within 0.1 of avesep away from ellipsoid surface
// 		}


		// harmonic here so that we never have a divergence because of changes in global geometry


		ff = -p->surfSpringSelf/p->surfDiffu*rij;
		en += 0.5*p->surfSpringSelf/p->surfDiffu*rijsq;

// 		ff = -p->surfSpringSelf*rij;
// 		en += 0.5*p->surfSpringSelf*rijsq;


		// negative forcing pointing at surf so pulling from surf to fund
		box->surffx[i] += ff*(spamx-fundx)*invrij;
		box->surffy[i] += ff*(spamy-fundy)*invrij;
		box->surffz[i] += ff*(spamz-fundz)*invrij;


		// find neighbours of i for cross interaction
		for(jj=1;jj<=box->surfNeigh[i][0];jj++)
		{
			// avoid double counting
			j = box->surfNeigh[i][jj];




			if(j>i && box->surfActive[j]==1)
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

				ff=0.0;


				// attractive exponential
// 				if(rij > 1.2*avesep && rij < 1.4*avesep)
// 				{
// 					en += p->surfSpringSelf*exp(1.0/(1.2*avesep - rij))/(1.4*avesep - rij);
// 					tmpexp = p->surfSpringSelf*exp(1.0/(1.2*avesep - rij))/(1.4*avesep - rij);
// 					tmpax = (1.2*avesep - rij);
// 					tmpbx = (1.4*avesep - rij);
//
//
// 					ff += - tmpexp*( 1.0/(tmpax*tmpax) + 1.0/tmpbx);
// 				}
//
// 				if(rij < 0.8*avesep && rij > 0.6*avesep)
// 				{
// 					en += p->surfSpringSelf*exp(1.0/(rij - 0.8*avesep))/(rij - 0.6*avesep);
//
// 					tmpexp = p->surfSpringSelf*exp(1.0/(rij - 0.8*avesep))/(rij - 0.6*avesep);
// 					tmpax = (rij - 0.8*avesep);
// 					tmpbx = (rij - 0.6*avesep);
//
// 					ff += - tmpexp*( -1.0/(tmpax*tmpax) - 1.0/tmpbx);
// 				}

// 				if(rij < 0.6*avesep || rij > 1.4*avesep)
// 				{
// 					printf("outside bounds for %d %d -- fixing with other -- should be (%f,%f) -- is %f \n",i,j,0.6*avesep,1.4*avesep,rij);
// 				}
//
// 				if(rij < 0.6*avesep)
// 				{
// 					patchcut = 0.8*avesep;
// 					ff += p->surfSpringCross*(p->surftauinv/p->surfDiffu)*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
// 					en += p->surfSpringCross*(p->surftauinv/p->surfDiffu)*pow((1-rij/patchcut),2.5);
//
// 				}
//
// 				if(rij > 1.4*avesep)
// 				{
// 					en += 0.5*p->surfSpringCross*rijsq;
// 					ff += -p->surfSpringCross*rij;
// 				}

				patchcut = 1.0*avesep;
				if(rij < patchcut)
				{
// 					ff = p->surfSpringCross*(p->surftauinv/p->surfDiffu)*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
// 					en += p->surfSpringCross*(p->surftauinv/p->surfDiffu)*pow((1-rij/patchcut),2.5);

					ff = p->surfSpringCross*(p->surftauinv)*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
					en += p->surfSpringCross*(p->surftauinv)*pow((1-rij/patchcut),2.5);



				}

// 				en += 0.5*p->surfSpringCross*(p->surftauinv/p->surfDiffu)*rijsq;
// 				ff += -p->surfSpringCross*(p->surftauinv/p->surfDiffu)*rij;

				en += 0.5*p->surfSpringCross*(p->surftauinv)*rijsq;
				ff += -p->surfSpringCross*(p->surftauinv)*rij;



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
	return en;
}

double calcSurfPatch(struct SYSTEM *box,struct PARAINPUT *p)
{
	// 	printf("done with wall-wall interactions -- %f \n",en);
	int i,j,ii,n1,jj,neipar,ch,kk,k;
	double en=0.0;

	double diffX,diffY,diffZ;
	double spamx,spamy,spamz,fundx,fundy,fundz;
	double rij,rijsq,invrij,r0,ff,rijanch;
	double patchcut,alphij,invwidthwell,tanhfact,lamintspam;

	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5);
	double minrijlam=4.0*p->ella;

	double surfplane[3];
	double surfplanenorm,surfplaneconst,planepconst,planepx,planepy,planepz;


	int closesurf;
	double minrijclose = 4.0*p->ella;
	double ath, wth;
	for(i=0;i<p->nAto;i++)
	{
		minrijlam=4.0*p->ella; // reset ato by ato
		minrijclose = 4.0*p->ella;
		box->closesurf[i] = -1;
		for(j=0;j<p->nSurfPoints;j++)
		{

			if(box->surfActive[j]==0)
				continue;
			// find the normal of the ellipse at surfFund
// 			surfplane[0] = 2.0*box->surfFundx[i]/(p->ella*p->ella);
// 			surfplane[1] = 2.0*box->surfFundy[i]/(p->ellb*p->ellb);
// 			surfplane[2] = 2.0*box->surfFundz[i]/(p->ellc*p->ellc);
// 			surfplanenorm = sqrt(surfplane[0]*surfplane[0]+surfplane[1]*surfplane[1]+surfplane[2]*surfplane[2]);
//
// 			surfplane[0] /= surfplanenorm;
// 			surfplane[1] /= surfplanenorm;
// 			surfplane[2] /= surfplanenorm;
//
// 			// find the constant specifying the plane parallel to original but passing through surfpoint
// 			surfplaneconst = (surfplane[0]*box->surfPointx[i] + surfplane[1]*box->surfPointy[i] + surfplane[2]*box->surfPointz[i]);
//
// 			// find distance between chrom centre and plane
// 			j = box->chromList[jj];
//
// 			rij = (surfplane[0]*box->X[j] + surfplane[1]*box->Y[j] + surfplane[2]*box->Z[j] - surfplaneconst)/sqrt(surfplane[0]*surfplane[0] + surfplane[1]*surfplane[1] + surfplane[2]*surfplane[2]);
//
//
// 			// find closest point on plane to chrom centre
// 			planepconst = ((surfplane[0]*box->X[j] + surfplane[1]*box->Y[j] + surfplane[2]*box->Z[j])- surfplaneconst);
// 			planepconst /= (surfplane[0]*surfplane[0] + surfplane[1]*surfplane[1] + surfplane[2]*surfplane[2]);
//
// 			planepx = box->X[j] + planepconst*surfplane[0];
// 			planepy = box->Y[j] + planepconst*surfplane[1];
// 			planepz = box->Z[j] + planepconst*surfplane[2];
//
// 			// get dx,dy,dz
// 			diffX = box->X[j] - planepx;
// 			diffY = box->Y[j] - planepy;
// 			diffZ = box->Z[j] - planepz;




			diffX = box->X[i] - box->surfPointx[j];
			diffY = box->Y[i] - box->surfPointy[j];
			diffZ = box->Z[i] - box->surfPointz[j];

			rijsq = diffX*diffX + diffY*diffY + diffZ*diffZ;

			// note the radii of the particles
			if(box->lamin[i] > 0.0 || box->type[i] > 0)
				patchcut = 0.5*avesep + box->radius[i];
			else if(box->lamin[i] == 0 && box->type[i]==0)
				patchcut = 1.0*(0.5*avesep + box->radius[i]);


// 			patchcut = pow(2.0,1.0/6.0)*(0.5*avesep + box->radius[j]);
// 			printf("distance of chrom %d from wall %d is %f %f \n",jj,i,rijsq,patchcut);
			if(rijsq <= sqr(pow(2.0,1.0/6.0)*patchcut)) // overlapping
			{
				rij = sqrt(rijsq);
				invrij = 1.0/rij;

				if(box->type[i] ==0 && rij < minrijclose)
				{
					minrijclose = rij;
					box->closesurf[i] = j;
// 					printf("%d %d %f %f %d\n",i,j,rij,minrijclose,box->closesurf[i]);
				}



// 				ff = 0.1*p->eps_hertz*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
// 				en += 0.1*p->eps_hertz*pow((1-rij/patchcut),2.5);

				ff = -4.0*(-12.0*pow(patchcut/rij,12.0)*invrij + 6.0*pow(patchcut/rij,6.0)*invrij);

				en += 4.0*(pow(patchcut/rij,12.0) - pow(patchcut/rij,6.0)) + 1.0;

// 				printf("chrom %d raw %d is %f away from boundary at %d, cut-off %f\n",jj,j,rij,i,patchcut);



// 				printf("updated forces for wall %d, chrom %d\n",i,j);

				/// check for lamin
// 				k = box->whichPatch[i]; -1 for chrom,

				if(box->lamin[i] > 0.0)
				{
// 					invwidthwell = 1.0/(0.01*box->lamin[i]*patchcut); // inverse of 0.1*patchcut
// 					tanhfact = tanh((rij - (1+0.01*box->lamin[i])*patchcut)*invwidthwell);

// 					invwidthwell = (p->lamin_scale/(box->lamin[i]*0.01*patchcut)) ;
// 					tanhfact = tanh((rij - (1 + (0.01*box->lamin[i])/p->lamin_scale)*patchcut)*invwidthwell);

					ath = 1.0 + 0.07*(p->lamin_scale/box->lamin[i]); // lamin_scale = box->lamin[i] for all i always
					wth = 0.07*(p->lamin_scale/box->lamin[i]);

					invwidthwell = 1.0/(wth*patchcut);
					tanhfact = tanh((rij - ath*patchcut)*invwidthwell);

					ff += -box->lamin[i]*(1.0 - tanhfact*tanhfact)*invwidthwell;

					en += 0.5*box->lamin[i]*tanhfact - 0.5*box->lamin[i];
// 						box->IMEne += 0.5*box->lamin[k]*tanhfact - 0.5*box->lamin[k];

					// figure out lamin modifier and add tanh attraction
// 					spamx = box->Gpatchposx[jj][kk] - box->surfFundx[i];
// 					spamy = box->Gpatchposy[jj][kk] - box->surfFundy[i];
// 					spamz = box->Gpatchposz[jj][kk] - box->surfFundz[i];
//
// 					rijanch = spamx*spamx + spamy*spamy + spamz*spamz;
// 					if(rijanch > 0.0)
// 						rijanch = sqrt(rijanch);

					rijanch = rij;

					lamintspam = 0.5 - 0.5*tanh((rijanch-1.4*patchcut)/(0.25*1.4*patchcut));

					box->isIntermingling[i] = MAX(box->isIntermingling[i],lamintspam);


// 							box->isIntermingling[k] = 1.0 - MIN(1.0,rij/(2.0*patchcut));
// 							box->isMembIntermingling[i] = 1.0 - MIN(1.0,rij/(2.0*patchcut));
// 						if(box->checkIntermingling[k] == 0.0)
// 						{
// 							box->isIntermingling[k] = 0.5 - 0.5*tanh((rijanch-1.4*patchcut)/(0.25*1.4*patchcut));
// 							box->checkIntermingling[k] = 1.0;
// 						}
					box->isMembIntermingling[j] = 0.5 - 0.5*tanh((rijanch-1.4*patchcut)/(0.25*1.4*patchcut));
					if(rij < minrijlam)
						minrijlam = rij;
				}


				box->surffx[j] += -ff*diffX*invrij;
				box->surffy[j] += -ff*diffY*invrij;
				box->surffz[j] += -ff*diffZ*invrij;

				box->fx[i] += ff*diffX*invrij;
				box->fy[i] += ff*diffY*invrij;
				box->fz[i] += ff*diffZ*invrij;

			}
		}


	}
	// additional osmotic/hydrostatic pressure keeping chromosomes and their patches internal
	box->osmoflag=0; // NOTE: Hard-coded prevention of code block below -- included here in case we lose track
	if(box->osmoflag==1)
	{
	for(i=0;i<p->nAto;i++)
	{

		if(box->type[i] == 0 && box->closesurf[i] != -1)
		{
			// find the normal of the ellipse at surfFund
			closesurf = box->closesurf[i];
			surfplane[0] = 2.0*box->surfFundx[closesurf]/(p->ella*p->ella);
			surfplane[1] = 2.0*box->surfFundy[closesurf]/(p->ellb*p->ellb);
			surfplane[2] = 2.0*box->surfFundz[closesurf]/(p->ellc*p->ellc);
			surfplanenorm = sqrt(surfplane[0]*surfplane[0]+surfplane[1]*surfplane[1]+surfplane[2]*surfplane[2]);

			surfplane[0] /= surfplanenorm;
			surfplane[1] /= surfplanenorm;
			surfplane[2] /= surfplanenorm;

			// find the constant specifying the plane parallel to original but passing through surfpoint
// 			surfplaneconst = (surfplane[0]*box->surfPointx[closesurf] + surfplane[1]*box->surfPointy[closesurf] + surfplane[2]*box->surfPointz[closesurf]);
			// here using plane through fund
			surfplaneconst = (surfplane[0]*box->surfFundx[closesurf] + surfplane[1]*box->surfFundy[closesurf] + surfplane[2]*box->surfFundz[closesurf]);

// 			// find distance between chrom centre and plane
// 			j = box->chromList[jj];

			rij = (surfplane[0]*box->X[i] + surfplane[1]*box->Y[i] + surfplane[2]*box->Z[i] - surfplaneconst)/sqrt(surfplane[0]*surfplane[0] + surfplane[1]*surfplane[1] + surfplane[2]*surfplane[2]);


			// find closest point on plane to chrom centre
			planepconst = ((surfplane[0]*box->X[i] + surfplane[1]*box->Y[i] + surfplane[2]*box->Z[i])- surfplaneconst);
			planepconst /= (surfplane[0]*surfplane[0] + surfplane[1]*surfplane[1] + surfplane[2]*surfplane[2]);

			planepx = box->X[i] + planepconst*surfplane[0];
			planepy = box->Y[i] + planepconst*surfplane[1];
			planepz = box->Z[i] + planepconst*surfplane[2];

// 			printf("%d %d -- %f %f %f -- %f %f %f -- %f -- %f %f %f\n",i,closesurf,box->X[i],box->Y[i],box->Z[i],planepx,planepy,planepz,planepconst,surfplane[0],surfplane[1],surfplane[2]);

			// get dx,dy,dz pointing away from plane
			diffX = box->X[i] - planepx;
			diffY = box->Y[i] - planepy;
			diffZ = box->Z[i] - planepz;

			rijsq = diffX*diffX + diffY*diffY + diffZ*diffZ;
// 			printf("dist %f %f -- %f %f %f\n",rij,sqrt(rijsq),diffX,diffY,diffZ);

			// note the radii of the particles
// 			patchcut = (0.5*avesep + box->radius[i]);

			patchcut = 0.5*box->radius[i];

// 			patchcut = pow(2.0,1.0/6.0)*(0.5*avesep + box->radius[j]);
// 			printf("distance of chrom %d from wall %d is %f %f \n",jj,i,rijsq,patchcut);
			ff = 0.0;
// 			if(rijsq <= sqr(pow(2.0,1.0/6.0)*patchcut)) // overlapping
			if(rijsq <= sqr(patchcut)) // overlapping
			{

				rij = sqrt(rijsq);
				invrij = 1.0/rij;

				ff = p->eps_hertz*2.5*pow((1 - rij/patchcut),1.5)/patchcut;
				en += p->eps_hertz*pow((1-rij/patchcut),2.5);

// 				ff = -4.0*(-12.0*pow(patchcut/rij,12.0)*invrij + 6.0*pow(patchcut/rij,6.0)*invrij);
//
// 				en += 4.0*(pow(patchcut/rij,12.0) - pow(patchcut/rij,6.0)) + 1.0;

				if(rijsq < sqr(0.4*box->radius[i]))
				{
					printf("!!!!!! OSMO !!! chrom %d patch %d pathological distance %f -- also %f  (raw pos %f %f %f - %f %f %f (fund %f %f %f))\n",i,closesurf,rij,sqr(box->X[i]/p->ella) + sqr(box->Y[i]/p->ellb) + sqr(box->Z[i]/p->ellc) -1,box->X[i],box->Y[i],box->Z[i],box->surfPointx[closesurf],box->surfPointy[closesurf],box->surfPointz[closesurf],box->surfFundx[closesurf],box->surfFundy[closesurf],box->surfFundz[closesurf]);
				}
			}

			box->fx[i] += ff*diffX*invrij;
			box->fy[i] += ff*diffY*invrij;
			box->fz[i] += ff*diffZ*invrij;

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
// 	en = en + calcWallInt(box,p);
	box->chromPot = en;

	en = en + calcSurfSurf(box,p);

	en = en + calcSurfPatch(box,p);
	box->surfPot = en - box->chromPot;

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

