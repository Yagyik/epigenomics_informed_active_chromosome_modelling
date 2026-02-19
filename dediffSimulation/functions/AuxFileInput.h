void writeChromdump(struct PARAINPUT *p, struct SYSTEM *box, struct APPOGGIO *a, int MCstep)
{
	// subroutine to write output in format that mimics LAMMPS
	// this is so that AnaCode can be used on configs printed out of this code as well
	int i,j,ii,jj;
	double spamx,spamy,spamz;
	double rotpatch[3];
	printf("writing to dump file\n");
	fprintf(box->fpLAMMPSdump,"ITEM: TIMESTEP\n");
	fprintf(box->fpLAMMPSdump,"%d\n",MCstep);		// TIMESTEP
	fprintf(box->fpLAMMPSdump,"ITEM: NUMBER OF ATOMS\n");
	fprintf(box->fpLAMMPSdump,"%d\n",2*p->nPatchTot+p->nChrom);
	fprintf(box->fpLAMMPSdump,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->lx,0.5*box->lx);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->ly,0.5*box->ly);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->lz,0.5*box->lz);
	fprintf(box->fpLAMMPSdump,"ITEM: ATOMS id type radius Transparency xu yu zu vx vy vz ox oy oz\n");
	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
// 			rotate(box->patchposx[i][j],box->patchposy[i][j],box->patchposz[i][j],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],rotpatch);
// // 			spamx = box->X[i] + box->patchposx[i][j];
// // 			spamy = box->Y[i] + box->patchposy[i][j];
// // 			spamz = box->Z[i] + box->patchposz[i][j];
//
// 			spamx = box->X[i] + rotpatch[0];
// 			spamy = box->Y[i] + rotpatch[1];
// 			spamz = box->Z[i] + rotpatch[2];

            spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
            // patches only have positions
			fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],0.5*box->radius[j],0.1,spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
            // beads have positions and velocities

            fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],box->radius[j],0.1,box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
		}
		fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->type[i]-1,box->radius[i],0.2,box->X[i],box->Y[i],box->Z[i],box->vx[i],box->vy[i],box->vz[i],box->axisx[i],box->axisy[i],box->axisz[i]);
	}
}

void writeDebug(struct PARAINPUT *p, struct SYSTEM *box,int MCstep)
{
	int i,j,ii,jj;

	// first find out how many points we print
	int nToPrint=0,countPrint;
	double spamx,spamy,spamz;
	for(i=0;i<p->nSurfPoints;i++)
	{
		if(box->surfFundy[i]>0.0 && box->surfFundx[i] > 0.0)
			nToPrint++;
	}

	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
			spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
			if(spamy>0.0 && spamx > 0.0 )
				nToPrint++;
		}
		if(box->Y[i]>0.0 && box->X[i] > 0.0)
			nToPrint++;
	}
	countPrint=0;
	fprintf(box->fpHalfPrint,"ITEM: TIMESTEP\n");
	fprintf(box->fpHalfPrint,"%d\n",MCstep);		// TIMESTEP
	fprintf(box->fpHalfPrint,"ITEM: NUMBER OF ATOMS\n");
	fprintf(box->fpHalfPrint,"%d\n",nToPrint);
	fprintf(box->fpHalfPrint,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(box->fpHalfPrint,"%f %f\n",-0.5*box->lx,0.5*box->lx);
	fprintf(box->fpHalfPrint,"%f %f\n",-0.5*box->ly,0.5*box->ly);
	fprintf(box->fpHalfPrint,"%f %f\n",-0.5*box->lz,0.5*box->lz);
	fprintf(box->fpHalfPrint,"ITEM: ATOMS id type radius xu yu zu\n");
	double print_rad=0.5;
	for(i=0;i<p->nSurfPoints;i++)
	{
		if(box->surfFundy[i]>0.0 && box->surfFundx[i]>0.0)
		{
			countPrint++;
			fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,p->nAto+4,0.5,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],0.0,0.0,0.0,0.0,0.0,0.0);

		}
	}
	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
            if(jj==ii)
				print_rad=0.01;
			else
				print_rad=0.5;
            spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
            // patches only have positions
			if(spamy>0.0 && spamx > 0.0)
			{
				countPrint++;
				fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,ii+1,print_rad,spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
				// beads have positions and velocities

				fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,ii+2,print_rad,box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
		}
		if(box->Y[i]>0.0 && box->X[i] > 0.0)
		{
			countPrint++;
			printf("chrom %d raw %d pos %f %f %f\n",ii,i,box->X[i],box->Y[i],box->Z[i]);
			fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,ii,box->radius[i],box->X[i],box->Y[i],box->Z[i],box->vx[i],box->vy[i],box->vz[i],box->axisx[i],box->axisy[i],box->axisz[i]);
		}
	}
}


void writeHalfOpen(struct PARAINPUT *p, struct SYSTEM *box,int MCstep)
{
	int i,j,ii,jj;

	// first find out how many points we print
	int nToPrint=0,countPrint;
	double spamx,spamy,spamz;
	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	for(i=0;i<p->nSurfPoints;i++)
	{
		if(box->surfFundy[i]>0.0) // & box->surfFundx[i] > 0.0)
			nToPrint++;
	}

	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
			spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
			if(spamy>0.0) // && spamx > 0.0 )
				nToPrint++;
		}
		if(box->Y[i]>0.0) // && box->X[i] > 0.0)
			nToPrint++;
	}
	countPrint=0;
	fprintf(box->fpHalfPrint,"ITEM: TIMESTEP\n");
	fprintf(box->fpHalfPrint,"%d\n",MCstep);		// TIMESTEP
	fprintf(box->fpHalfPrint,"ITEM: NUMBER OF ATOMS\n");
	fprintf(box->fpHalfPrint,"%d\n",nToPrint);
	fprintf(box->fpHalfPrint,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(box->fpHalfPrint,"%f %f\n",-0.5*box->lx,0.5*box->lx);
	fprintf(box->fpHalfPrint,"%f %f\n",-0.5*box->ly,0.5*box->ly);
	fprintf(box->fpHalfPrint,"%f %f\n",-0.5*box->lz,0.5*box->lz);
	fprintf(box->fpHalfPrint,"ITEM: ATOMS id type radius Transparency xu yu zu vx vy vz fx fy fz\n");
	double print_rad=0.5;
	for(i=0;i<p->nSurfPoints;i++)
	{
		if(box->surfFundy[i]>0.0) // && box->surfFundx[i]>0.0)
		{
			countPrint++;
			fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,p->nAto+4,0.5*avesep,0.5,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],0.0,0.0,0.0,0.0,0.0,0.0);

		}
	}
	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
			print_rad = box->radius[j];
// 			printf("printing radius for %d patch %d on chrom %d as %f\n",j,jj,ii,print_rad);
//             if(jj==ii)
// 				print_rad=0.1*box->radius[j];
// 			else
// 				print_rad=0.5*box->radius[j];
            spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
            // patches only have positions
			if(spamy>0.0) // && spamx > 0.0)
			{
				countPrint++;
				if(box->patchPointer[ii][jj]==p->nAto+2)
				{
				fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],0.5*box->radius[j],0.95,spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
				// beads have positions and velocities

				fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],box->radius[j],0.95,box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
				}
				else
				{
				fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],0.5*box->radius[j],0.8,spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
				// beads have positions and velocities

				fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],box->radius[j],0.4,box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
				}

// 				fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,box->patchPointer[ii][jj],0.5*print_rad,0.8,spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
// 				// beads have positions and velocities
//
// 				fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,box->patchPointer[ii][jj],print_rad,0.1,box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
		}
		if(box->Y[i]>0.0 ) //&& box->X[i] > 0.0)
		{
			countPrint++;
// 			printf("chrom %d raw %d pos %f %f %f\n",ii,i,box->X[i],box->Y[i],box->Z[i]);
			fprintf(box->fpHalfPrint,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",countPrint,box->type[i]-1,box->radius[i],0.4,box->X[i],box->Y[i],box->Z[i],box->vx[i],box->vy[i],box->vz[i],box->axisx[i],box->axisy[i],box->axisz[i]);
		}
	}
}

void writeTranspdump(struct PARAINPUT *p, struct SYSTEM *box, struct APPOGGIO *a, int MCstep)
{
	// subroutine to write output in format that mimics LAMMPS
	// this is so that AnaCode can be used on configs printed out of this code as well
	int i,j,ii,jj;
	double spamx,spamy,spamz;
	double rotpatch[3];

	double spam = 0.33333*(pow(p->ellai*p->ellbi,1.6) + pow(p->ellai*p->ellci,1.6) + pow(p->ellbi*p->ellci,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfBasicPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)




	printf("writing to dump file\n");
	fprintf(box->fpLAMMPSdump,"ITEM: TIMESTEP\n");
	fprintf(box->fpLAMMPSdump,"%d\n",MCstep);		// TIMESTEP
	fprintf(box->fpLAMMPSdump,"ITEM: NUMBER OF ATOMS\n");
	fprintf(box->fpLAMMPSdump,"%d\n",2*p->nPatchTot+p->nChrom + p->nSurfPoints);
	fprintf(box->fpLAMMPSdump,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->lx,0.5*box->lx);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->ly,0.5*box->ly);
	fprintf(box->fpLAMMPSdump,"%f %f\n",-0.5*box->lz,0.5*box->lz);
	fprintf(box->fpLAMMPSdump,"ITEM: ATOMS id type radius Transparency xu yu zu vx vy vz ox oy oz\n");


	for(i=0;i<p->nSurfPoints;i++)
	{
		fprintf(box->fpLAMMPSdump,"%d %d %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,p->nAto+4,0.5*avesep,0.75*(1.0 - box->isMembIntermingling[i]),box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],0.0,0.0,0.0,0.0,0.0,0.0);
// 		fprintf(box->fpSurfPoints,"%d %d %f %lf %lf %lf\n",i,200,0.5*avesep,box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}

	for(ii=0;ii<p->nChrom;ii++)
	{
        i = box->chromList[ii];


		// for this chromosome, use the axis and the angle to find each patch
		for(jj=1;jj<=box->patchList[ii][0];jj++)
		{
            j = box->patchList[ii][jj];
// 			rotate(box->patchposx[i][j],box->patchposy[i][j],box->patchposz[i][j],box->axisx[i],box->axisy[i],box->axisz[i],box->cosRotAngle[i],rotpatch);
// // 			spamx = box->X[i] + box->patchposx[i][j];
// // 			spamy = box->Y[i] + box->patchposy[i][j];
// // 			spamz = box->Z[i] + box->patchposz[i][j];
//
// 			spamx = box->X[i] + rotpatch[0];
// 			spamy = box->Y[i] + rotpatch[1];
// 			spamz = box->Z[i] + rotpatch[2];

            spamx = box->Gpatchposx[ii][jj];
			spamy = box->Gpatchposy[ii][jj];
			spamz = box->Gpatchposz[ii][jj];
            // patches only have positions
			if(box->patchPointer[ii][jj]==p->nAto+2)
			{
			fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],0.5*box->radius[j],0.95,spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
            // beads have positions and velocities

            fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],box->radius[j],0.95,box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
			else
			{
			fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],0.5*box->radius[j],0.9*(1.0-box->isIntermingling[j]),spamx,spamy,spamz,0.0,0.0,0.0,0.0,0.0,0.0);
            // beads have positions and velocities

            fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->patchPointer[ii][jj],box->radius[j],0.9*(1.0- box->isIntermingling[j]),box->X[j],box->Y[j],box->Z[j],box->vx[j],box->vy[j],box->vz[j],box->fx[j],box->fy[j],box->fz[j]);
			}
		}
		fprintf(box->fpLAMMPSdump,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i+1,box->type[i]-1,box->radius[i],0.75,box->X[i],box->Y[i],box->Z[i],box->vx[i],box->vy[i],box->vz[i],box->axisx[i],box->axisy[i],box->axisz[i]);
	}
}
