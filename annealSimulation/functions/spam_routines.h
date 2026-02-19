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



