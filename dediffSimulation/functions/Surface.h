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

void vnlist(struct PARAINPUT *p, struct SYSTEM *box)
{
	int i,j;
	int neigcount_i,neigcount_j;

	double rX,rY,rZ,r;

	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	printf("recreating full surf neighlist\n");

	for(i=0;i<p->nSurfPoints;i++)
	{
		box->surfNeigh[i][0]=0;
	}
	// NOTE: populate surface neighlist
	for(i=0;i<p->nSurfPoints;i++)
	{
		neigcount_i=box->surfNeigh[i][0];

		for(j=0;j<p->nSurfPoints;j++)
		{
			if(j>i)
			{
				rX = box->surfFundx[j] - box->surfFundx[i];
				rY = box->surfFundy[j] - box->surfFundy[i];
				rZ = box->surfFundz[j] - box->surfFundz[i];

				r=rX*rX+rY*rY+rZ*rZ;

				//if(r < sqr(avesep+SKIN)) // close enough??
				if(r<sqr(1.5*avesep))
				{
					neigcount_i++;
					box->surfNeigh[i][neigcount_i]=j;
					neigcount_j = box->surfNeigh[j][0] + 1;
					box->surfNeigh[j][neigcount_j]=i;
					box->surfNeigh[j][0] = neigcount_j;
				}
			}
		}
		box->surfNeigh[i][0]=neigcount_i;
// 		print("setting neig %d %d\n",i,box->surfNeig[i][0]);


		box->surfdispX[i]=0.0;
		box->surfdispY[i]=0.0;
		box->surfdispZ[i]=0.0;
	}
	printf("finished surf neighlist\n");

	box->maxsurfdispsq=0.0;
}

void vnlistpar(int par,struct SYSTEM *box,double delx,double dely,double delz,struct PARAINPUT *p)
{
	int j,lcount=0,neigcount;
	double rijsq,Xij,Yij,Zij;
	double displacesq, neipar;
	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)

	displacesq=sqr(box->surfdispX[par]+delx)+sqr(box->surfdispY[par]+dely)+sqr(box->surfdispZ[par]+delz);

	if(displacesq>sqr(0.1*0.5*avesep))
	{

		for(j=0;j<p->nSurfPoints;j++)
		{
			if(j!=par)
			{
				Xij=box->surfFundx[j]-box->surfFundx[par];
// 				Xij-=box->lx*lround(Xij/box->lx);
				Yij=box->surfFundy[j]-box->surfFundy[par];
// 				Yij-=box->ly*lround(Yij/box->ly);
				Zij=box->surfFundz[j]-box->surfFundz[par];
// 				Zij-=box->lz*lround(Zij/box->lz);

				rijsq=Xij*Xij+Yij*Yij+Zij*Zij;
				//if(r < sqr(avesep+SKIN)) // close enough??
				if(rijsq<sqr(1.5*avesep))
				{
					lcount++;
					box->surfneigpar[lcount]=j;
				}
			}
		}
		box->surfneigpar[0]=lcount;
		box->surfdispX[par] = 0.0;
		box->surfdispY[par] = 0.0;
		box->surfdispZ[par] = 0.0;
// 		printf("surf par %d moved beyond threshold -- calculating neigh list for it -- %d neighbours\n",par,lcount);
	}
	else
	{
		neigcount=box->surfNeigh[par][0];
		for(j=1;j<=neigcount;j++)
		{
			box->surfneigpar[j]=box->surfNeigh[par][j];
		}
		box->surfneigpar[0]=neigcount;


	}
}

void update(int par,struct PARAINPUT *p,struct SYSTEM *box)
{
	double displacesq;
	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	displacesq=sqr(box->surfdispX[par])+sqr(box->surfdispY[par])+sqr(box->surfdispZ[par]);
	box->maxsurfdispsq=MAX(box->maxsurfdispsq,displacesq);

	if(box->maxsurfdispsq>sqr(0.1*0.5*avesep))
	{
		printf("updating neiglist cus %d moved %f - %f %f %f\n",par,box->maxsurfdispsq,box->surfdispX[par],box->surfdispY[par],box->surfdispZ[par]);
		box->maxsurfdispsq=0.0;

		vnlist(&(*p),&(*box));
	}
}

double FundEnSurf(struct PARAINPUT *p, struct SYSTEM *box,int par)
{

	// subroutine to calculate energy cost of current set of fundamental points
	int i,j,ii,n1,jj,neipar,ch;
	double en=0.0;

	double diffX,diffY,diffZ;
	double spamx,spamy,spamz,fundx,fundy,fundz;
	double rij,rijsq,invrij,r0,ff,patchcut;
	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)


	fundx = box->surfFundx[par];
	fundy = box->surfFundy[par];
	fundz = box->surfFundz[par];
// 		printf("surf %f %f %f -- fund %f %f %f\n",spamx,spamy,spamz,fundx,fundy,fundz);

	rijsq = fabs(sqr(fundx/p->ella) + sqr(fundy/p->ellb) + sqr(fundz/p->ellc) - 1);
	if(rijsq > avesep*avesep)
	{
		rij=sqrt(rijsq);
		invrij = 1.0/rij;
// 		printf("we should never see this, cus of new planar displacement routine - %d %f -- %f %f %f (%f %f %f)\n",par,sqrt(rijsq),fundx,fundy,fundz,p->ella,p->ellb,p->ellc);
	}
// 	en += 10*p->surfSpringSelf*rijsq;
	en += 128*p->surfSpringSelf*rijsq;
// 	printf("%d %d --self contri %f %f -- sep %f\n",par,box->surfneigpar[0],rij,p->surfSpringSelf*rijsq,avesep);
	// find neighbours of i for cross interaction
	for(jj=1;jj<=box->surfneigpar[0];jj++)
	{
		// avoid double counting
		j = box->surfneigpar[jj];
// 			printf("%d neighbour - %d of surf par %d\n",box->surfNeigh[par][0],j,par);
// 	for(j=0;j<p->nSurfPoints;j++)
// 	{
		if(j>par)
		{
			// find cross interaction

			diffX = box->surfFundx[j] - box->surfFundx[par];
			diffY = box->surfFundy[j] - box->surfFundy[par];
			diffZ = box->surfFundz[j] - box->surfFundz[par];

			rijsq = diffX*diffX + diffY*diffY + diffZ*diffZ;
			rij = sqrt(rijsq);
			invrij = 1.0/rij;
		// NOTE: WCA interact
// 			if(rij <= pow(2.0,1.0/6.0)*avesep) // WCA cut-off
// 			{
// // 				printf("wca contri %f %f %f\n",rij,avesep,p->surfSpringCross*(4.0*(pow(avesep/rij,12.0) - pow(avesep/rij,6.0)) + 1.0));
// 				en += p->surfSpringCross*(4.0*(pow(avesep/rij,12.0) - pow(avesep/rij,6.0)) + 1.0);
// 			}
			if( rij > 0.0 && rij < 1.2*avesep)
			{
// 				en += 200*p->surfSpringSelf*pow((1 - rij/avesep),2.5);
				en += 32*p->surfSpringSelf*pow((1 - rij/(1.2*avesep)),2.5);
// 				printf("neighbour contri %d %d -- %f %f\n",par,j,rij,p->surfSpringCross*pow((1 - rij/avesep),2.5));
			}
			en += 0.5*p->surfSpringSelf*tanh((rij - 1.0*avesep)/(0.25*avesep)) - 0.5*p->surfSpringSelf;
			/// add a tanh attractive component -- make sure lengthscale remains avesep

		}
	}

	return en;

}


void MCRelaxSurf(struct PARAINPUT *p, struct SYSTEM *box,int tsim)
{

	// subroutine to perform MC relaxation for surf Fund points to the new geometry

	double debugen=0.0;
	double en,olden,deltaen,rand,prob;
	double oldX,oldY,oldZ,rX,rY,rZ,r;
	double delx,dely,delz,spamx,spamy,spamdist;
	int i,ii,j,jj,random,neigcount_i,neigcount_j,novlap;

	// rotation matrix things
	double surfplane[3],rotmat[9],surfdisp[3];
	double surfplanenorm,costheta,sintheta,u1,u2;

	vnlist(p,box); // update neighbourlist

	en=olden=0.0;
	for(i=0;i<p->nSurfPoints;i++)
	{
		vnlistpar(i,box,0.0,0.0,0.0,p);
		en +=FundEnSurf(p,box,i);
//  		printf("en %d %f\n",i,en);
	}
	olden=en;
	double storeen = en;
	printf("doing MC relax, starting en %f -set old %f\n",en,olden);

	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	for(int t=0;t<tsim;t++)
	{
		for(i=0;i<p->nSurfPoints;i++)
		{
			random=(int) (p->nSurfPoints*drand48());
			if(random <0.0)
				random=0;
			if(random>=p->nSurfPoints)
			   random=p->nSurfPoints-1;

			// create vnlistpar for zero displace
			vnlistpar(random,box,0.0,0.0,0.0,p); // un-updated neiglist
			olden = FundEnSurf(p,box,random); //

			// displace this particle

			oldX = box->surfFundx[random];
			oldY = box->surfFundy[random];
			oldZ = box->surfFundz[random];


			/// naive displacement
// 			delx = 0.01*(2.0*drand48()-1.0);
// 			dely = 0.01*(2.0*drand48()-1.0);
// 			delz = 0.01*(2.0*drand48()-1.0);


			/// implement a 2D displacement in the tangent plane at the current point -- done by using a delx,dely,0 displacement and identifying the isometry with current tangent plane
			spamx = (2.0*drand48()-1.0);
			spamy = (2.0*drand48()-1.0);

			surfplane[0] = 2.0*oldX/(p->ella*p->ella);
			surfplane[1] = 2.0*oldY/(p->ellb*p->ellb);
			surfplane[2] = 2.0*oldZ/(p->ellc*p->ellc);
			surfplanenorm = sqrt(surfplane[0]*surfplane[0]+surfplane[1]*surfplane[1]+surfplane[2]*surfplane[2]);

			surfplane[0] /= surfplanenorm;
			surfplane[1] /= surfplanenorm;
			surfplane[2] /= surfplanenorm;

			// compose the rotation matrix -- costheta = plo.pln/(|plo||pln|) -- both norms are 1
			costheta = (0.0*surfplane[0] + 0.0*surfplane[1] + 1.0*surfplane[2]);
			sintheta = sqrt(1-costheta*costheta);
			// uvec = plo x pln / (|plo x pln|) - pln [0,0,1], simplifying the cross product (and the denominator which is the norm of the cross product
			u1 = surfplane[0]/(sqrt(surfplane[0]*surfplane[0] + surfplane[1]*surfplane[1]));
			u2 = surfplane[1]/(sqrt(surfplane[0]*surfplane[0] + surfplane[1]*surfplane[1]));

			rotmat[0] = costheta + u1*u1*(1-costheta);
			rotmat[1] = u1*u2*(1-costheta);
			rotmat[2] = u2*sintheta;	box->maxsurfdispsq = 0.0;

			rotmat[3] = rotmat[1];
			rotmat[4] = costheta + u2*u2*(1-costheta);
			rotmat[5] = -u1*sintheta;
			rotmat[6] = -rotmat[2];
			rotmat[7] = -rotmat[5];
			rotmat[8] = costheta;

			// R.vo gives vn
			delx = rotmat[0]*spamx + rotmat[1]*spamy;
			dely = rotmat[3]*spamx + rotmat[4]*spamy;
			delz = rotmat[6]*spamx + rotmat[7]*spamy;
			spamdist = sqrt(delx*delx + dely*dely+delz*delz);


			delx /= 100.0*spamdist;
			dely /= 100.0*spamdist;
			delz /= 100.0*spamdist;
// 			printf("%d %d -- vecnorm %f %f %f -- %f\n",t,random,delx,dely,delz,sqrt(delx*delx + dely*dely+delz*delz));

			// delx,dely,delz here are small displacements in the tangent plane of the ellipse at oldx,oldy,oldz.

// 			printf("%d %d old ellipse %f %f %f - %f\n",t,random,box->surfFundx[random],box->surfFundy[random],box->surfFundz[random],
// 				  fabs(sqr(box->surfFundx[random]/p->ella) + sqr(box->surfFundy[random]/p->ellb) + sqr(box->surfFundz[random]/p->ellc) - 1));

			delx += 0.005*(2.0*drand48()-1.0);
			dely += 0.005*(2.0*drand48()-1.0);
			delz += 0.005*(2.0*drand48()-1.0);

			box->surfFundx[random] += delx;
			box->surfFundy[random] += dely;
			box->surfFundz[random] += delz;

// 			printf("%d %d new ellipse %f %f %f - %f\n",t,random,box->surfFundx[random],box->surfFundy[random],box->surfFundz[random],
// 				  fabs(sqr(box->surfFundx[random]/p->ella) + sqr(box->surfFundy[random]/p->ellb) + sqr(box->surfFundz[random]/p->ellc) - 1));


			// update vnlistpar here
			vnlistpar(random,box,delx,dely,delz,p); // possibly updated

			en = FundEnSurf(p,box,random);

			deltaen = en - olden;
// 			printf("par %d delta E %f\n",random,deltaen);

			prob = exp(-1000*deltaen); // temp =0.1 (NOTE: Re-parameterise)

			rand = drand48();
			if(rand < prob)
			{
				// accept
				storeen +=deltaen;
// 				printf("accepting %f %f-- %f %f %d\n",olden,en,deltaen,storeen,random);
				// update dispX,dispY,dispZ here
				box->surfdispX[random] +=delx;
				box->surfdispY[random] +=dely;
				box->surfdispZ[random] +=delz;

				// check if neighlist needs to be updated
				update(random,p,box);

			}
			else
			{
				box->surfFundx[random] = oldX;
				box->surfFundy[random] = oldY;
				box->surfFundz[random] = oldZ;
			}

		}

		if(t%((int)(0.02*tsim)) == 0)
		{
			debugen=0.0;
			vnlist(p,box); // update neighbourlist
			for(i=0;i<p->nSurfPoints;i++)
			{
				vnlistpar(i,box,0.0,0.0,0.0,p);
				debugen += FundEnSurf(p,box,i);
// 				printf("en %d %f\n",i,FundEnSurf(p,box,i));
			}
			printf("%d current MC energy -- %f %f\n",t,debugen,avesep);
		}
	}
	printf("final MC energy after %d iter -- %f\n",tsim,debugen);

	vnlist(p,box); // update neighbourlist
	printf("moving on from neighlist\n");

	// at the end we (may) need to update the surf points to the fund points and reset forces to zero
}


void findGrid(struct PARAINPUT *p, struct SYSTEM *box)
{
	// subroutine to assign a rectangular grid box to each surfFund point
	int i,j;
	double boxLx=2.0*(p->ella+2);
	double boxLy=2.0*(p->ellb+2);
	double boxLz=2.0*(p->ellc+2);
	double boxcentx,boxcenty,boxcentz;

	double surfplane[3],rotmat[9],surfdisp[3];
	double surfplanenorm,costheta,sintheta,u1,u2;
	double delx,dely,delz,spamx,spamy,spamdist;
	int novlap;


	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)

	for(i=0;i<p->nSurfPoints;i++)
	{
		// find the grid box of these surfPointx

		box->igridx[i] = (int) ((box->surfFundx[i] - (-p->ella-2))/boxLx*p->nGrid);
		box->igridy[i] = (int) ((box->surfFundy[i] - (-p->ellb-2))/boxLy*p->nGrid);
		box->igridz[i] = (int) ((box->surfFundz[i] - (-p->ellc-2))/boxLz*p->nGrid);

		// find positions wrt centre of box
// 		printf("%d pos %f %f %f -- grid %d %d %d\n",i,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i],box->igridx[i],box->igridy[i],box->igridz[i]);
		// first find box centre
		boxcentx = -p->ella-2 + (box->igridx[i]+0.5)*boxLx/(1.0*p->nGrid); // xmin + (i+0.5)*Lx/nx (Lx/nx=dx -- dx is box size)
		boxcenty = -p->ellb-2 + (box->igridy[i]+0.5)*boxLy/(1.0*p->nGrid); // xmin + (i+0.5)*Lx/nx (Lx/nx=dx)
		boxcentz = -p->ellc-2 + (box->igridz[i]+0.5)*boxLz/(1.0*p->nGrid); // xmin + (i+0.5)*Lx/nx (Lx/nx=dx)

		// find distance from box centre

		box->centSurfx[i] = (box->surfFundx[i] - boxcentx)*p->nGrid/boxLx; // distance from centre as fraction of box length
		box->centSurfy[i] = (box->surfFundy[i] - boxcenty)*p->nGrid/boxLy; // distance from centre as fraction of box length
		box->centSurfz[i] = (box->surfFundz[i] - boxcentz)*p->nGrid/boxLz; // distance from centre as fraction of box length

		if(fabs(box->centSurfx[i]) >0.5 || fabs(box->centSurfy[i]) >0.5 || fabs(box->centSurfz[i]) >0.5 )
		{
			printf("error centering surf Fund %d %f %f %f\n",i,box->centSurfx[i],box->centSurfy[i],box->centSurfz[i]);
			exit(0);
		}

// 		printf("for %d this %f should be equal to (igrid x dx + min) %fx%f + %f + (cent x dx) %fx%f = %f\n",i,box->surfFundx[i],box->igridx[i]+0.5,boxLx/(1.0*p->nGrid),-p->ella,box->centSurfx[i],boxLx/(1.0*p->nGrid),-p->ella + (box->igridx[i]+0.5)*boxLx/(1.0*p->nGrid)+box->centSurfx[i]*boxLx/(1.0*p->nGrid));
// 		printf("check y %f %f\n",box->surfFundy[i],-p->ellb + (box->igridy[i]+0.5)*boxLy/(1.0*p->nGrid)+box->centSurfy[i]*boxLy/(1.0*p->nGrid));
// 		printf("check z %f %f\n",box->surfFundz[i],-p->ellc + (box->igridz[i]+0.5)*boxLz/(1.0*p->nGrid)+box->centSurfz[i]*boxLz/(1.0*p->nGrid));

	}
	for(i=0;i<p->nSurfPoints;i++)
	{


		/// naive surfpoint placement
// 		box->surfPointx[i] = box->surfFundx[i];
// 		box->surfPointy[i] = box->surfFundy[i];
// 		box->surfPointz[i] = box->surfFundz[i];


		// first we place the surfpoint outside the surf fund along the normal
		surfplane[0] = 2.0*box->surfFundx[i]/(p->ella*p->ella);
		surfplane[1] = 2.0*box->surfFundy[i]/(p->ellb*p->ellb);
		surfplane[2] = 2.0*box->surfFundz[i]/(p->ellc*p->ellc);
		surfplanenorm = sqrt(surfplane[0]*surfplane[0]+surfplane[1]*surfplane[1]+surfplane[2]*surfplane[2]);

		surfplane[0] /= surfplanenorm;
		surfplane[1] /= surfplanenorm;
		surfplane[2] /= surfplanenorm;

		box->surfPointx[i] = box->surfFundx[i] + surfplane[0]*avesep;
		box->surfPointy[i] = box->surfFundy[i] + surfplane[1]*avesep;
		box->surfPointz[i] = box->surfFundz[i] + surfplane[2]*avesep;


		// place surfpoints such that the overlaps with patches are minimised.
		// check for overlaps with internal entities and push out along normal
		novlap=1;
		while(novlap>0)
		{
			novlap=0;
			for(j=0;j<p->nAto;j++)
			{
				delx = box->surfPointx[i] - box->X[j];
				dely = box->surfPointy[i] - box->Y[j];
				delz = box->surfPointz[i] - box->Z[j];

				spamdist = delx*delx + dely*dely + delz*delz;
				if(spamdist < 1.2*(0.5*avesep + box->radius[j]))
				{
					novlap +=1;
					break;
				}
			}
			if(novlap==1)
			{
				box->surfPointx[i] += 0.5*surfplane[0]*avesep;
				box->surfPointy[i] += 0.5*surfplane[1]*avesep;
				box->surfPointz[i] += 0.5*surfplane[2]*avesep;
			}
		}


		// force reset before next integration

		box->surffx[i] = 0.0;
		box->surffy[i] = 0.0;
		box->surffz[i] = 0.0;

	}

}

void scaleSurf(struct PARAINPUT *p, struct SYSTEM *box)
{
	// subroutine to find the scale surfFund positions

	int i,j;
	double boxLx=2.0*(p->ella+2);
	double boxLy=2.0*(p->ellb+2);
	double boxLz=2.0*(p->ellc+2);

// 	printf("scaling surf points to %f %f %f %d\n",p->ella,p->ellb,p->ellc,p->nGrid);
	for(i=0;i<p->nSurfPoints;i++)
	{
// 		printf("%d old fund %f %f %f\n",i,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
		box->oldsurfFundx[i] = box->surfFundx[i];
		box->oldsurfFundy[i] = box->surfFundy[i];
		box->oldsurfFundz[i] = box->surfFundz[i];
		box->surfFundx[i] = -p->ella-2 + (box->igridx[i]+0.5)*boxLx/(1.0*p->nGrid) + box->centSurfx[i]*boxLx/(1.0*p->nGrid);
		box->surfFundy[i] = -p->ellb-2 + (box->igridy[i]+0.5)*boxLy/(1.0*p->nGrid) + box->centSurfy[i]*boxLy/(1.0*p->nGrid);
		box->surfFundz[i] = -p->ellc-2 + (box->igridz[i]+0.5)*boxLz/(1.0*p->nGrid) + box->centSurfz[i]*boxLz/(1.0*p->nGrid);

// 		printf("%d new fund %f %f %f\n",i,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}

	// need to move the surfpoints accordingly :: NOTE: Only when using a diverging potential
// 	for(i=0;i<p->nSurfPoints;i++)
// 	{
// 		box->surfPointx[i] += box->surfFundx[i] - box->oldsurfFundx[i]; // update position by moving in same direction as fund (avoiding breaking entropic spring)
// 		box->surfPointy[i] += box->surfFundy[i] - box->oldsurfFundy[i];
// 		box->surfPointz[i] += box->surfFundz[i] - box->oldsurfFundz[i];
// 	}
}

void genSurf(FILE *fp,struct PARAINPUT *p, struct SYSTEM *box)
{
	// subroutine generates uniformly distributed points on surface of ellipsoid
	// main parameters are number of points, and packing fraction on surface
	// take the ellipsoid semi axes a>=b>=c (a=b=c implies sphere)

	int i,j;
	double curr[4];
	double rX,rY,rZ,r,pkeep,random;
	int atomdone,atomnotdone;

	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	for(i=0;i<p->nSurfPoints;i++)
	{
		atomdone=0;
		while(atomdone==0)
		{

			curr[0] = gaussVar(0,1);
			curr[1] = gaussVar(0,1);
			curr[2] = gaussVar(0,1);
			curr[3] = sqrt(curr[0]*curr[0] + curr[1]*curr[1] + curr[2]*curr[2]);
			curr[0] /= curr[3];
			curr[1] /= curr[3];
			curr[2] /= curr[3];
			atomnotdone=0;
			for(j=0;j<i;j++)
			{
//                     ch=structbox->type[i]+structbox->type[j]-2;
				rX = p->ella*curr[0] - box->surfPointx[j]; // transform to ellipsoid position and check separation
				rY = p->ellb*curr[1] - box->surfPointy[j];
				rZ = p->ellc*curr[2] - box->surfPointz[j];


				r=rX*rX+rY*rY+rZ*rZ;
//                     r=rX*rX+rY*rY;
				if(r<sqr(0.5*avesep)) // too close
				{
// 						printf("%d %d %f %f\n",i,j,r,minDist2*p->sigma[ch]);
// 					printf("%d separation rejection %d %f %f\n",i,j,sqrt(r),avesep);
					atomnotdone=1;
				}
			}


			// NOTE: now do the ellipsoid scaling and discard point with probability
// 			pkeep = pow(p->ella*p->ellc*p->ellb*curr[1],2.0) + pow(p->ella*p->ellb*p->ellc*curr[2],2.0) + pow(p->ellb*p->ellc*p->ella*curr[0],0.5);
			pkeep = pow(p->ella*p->ellc*curr[1],2.0) + pow(p->ella*p->ellb*curr[2],2.0) + pow(p->ellb*p->ellc*curr[0],2.0);
			pkeep = sqrt(pkeep);
// 			printf("%d current pos %f %f %f and radius, %f - surf area %f\n",i,p->ella*curr[0],p->ellb*curr[1],p->ellc*curr[2],p->ella*curr[0]*p->ella*curr[0] + p->ellb*curr[1]*p->ellb*curr[1] + p->ellc*curr[2]*p->ellc*curr[2],pkeep);
			pkeep = pkeep/(p->ella*p->ellb); // max when a and b are the largest axes
			random=drand48();
			if(random>pkeep) // discard point
			{
				atomnotdone=1;
// 				printf("%d polar rejection %f %f\n",i,random,pkeep);
			}


			atomdone=1-atomnotdone;
// 			box->surfPointx[i] = p->ella*curr[0]; // keeps getting replaced at iterations
// 			box->surfPointy[i] = p->ellb*curr[1]; // store transformed ellipsoidal positions
// 			box->surfPointz[i] = p->ellc*curr[2];

			box->surfFundx[i] = p->ella*curr[0];
			box->surfFundy[i] = p->ellb*curr[1];
			box->surfFundz[i] = p->ellc*curr[2];


		}
// 		box->surfFundx[i] = box->surfPointx[i];
// 		box->surfFundy[i] = box->surfPointy[i];
// 		box->surfFundz[i] = box->surfPointz[i];

		box->surfPointx[i] = box->surfFundx[i];
		box->surfPointy[i] = box->surfFundy[i];
		box->surfPointz[i] = box->surfFundz[i];




		// force reset before next integration

		box->surffx[i] = 0.0;
		box->surffy[i] = 0.0;
		box->surffz[i] = 0.0;



//  		printf("%d surf %f %f %f -- fund %f %f %f\n",i,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}

	// remake neighbourlist
// 	vnlist(p,box);

	for(i=0;i<p->nSurfPoints;i++)
	{
		box->oldsurfFundx[i] = box->surfFundx[i];
		box->oldsurfFundy[i] = box->surfFundy[i];
		box->oldsurfFundz[i] = box->surfFundz[i];
// 		box->surfPointx[i] = box->surfFundx[i];
// 		box->surfPointy[i] = box->surfFundy[i];
// 		box->surfPointz[i] = box->surfFundz[i];
		fprintf(fp,"%d %d\n",i,box->surfNeigh[i][0]);
	}


}
