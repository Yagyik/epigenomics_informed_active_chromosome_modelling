



void readSource(FILE *fp, double **source,char **fileList,char **IADfileList,char **ChromPropsfileList,char **IPDfileList, int nSyst)
{
	char *line;
	int i,length,not_space;
	int count,vals;

	line=(char *) malloc(sizeof(char) * LINESIZE);
	if(!line)
	{
		printf("Read Source: Couldn't allocate space for line\n");
		exit(-1);
	}

	rewind(fp);
	line=fgets(line,LINESIZE,fp);
	count=0;		// Count how many lines are there

	while(line!=NULL && count<nSyst)
	{
		if(line[0]=='#')	// It is a comment
		{
			line=fgets(line,LINESIZE,fp);

		}
		else
		{
			// verify it's not a blank line
			length=strlen(line);		// string length
			printf("%s",line);
			not_space=0;
			i=0;
			while(!not_space && i<length)
			{
				if(!isspace(line[i]))	// Checks whether line[i] is a white-space character.
				{
					not_space=1;
				}
				i++;
			}

				if(not_space)
				{
					// Each row maps one MPI rank to output and input files.
					vals=sscanf(line,"%lf %lf %lf %s %s %s %s",&source[count][0],&source[count][1],&source[count][2],fileList[count],IADfileList[count],ChromPropsfileList[count],IPDfileList[count]);
					printf("%f %f %f %s %s %s %s\n",source[count][0],source[count][1],source[count][2],fileList[count],IADfileList[count],ChromPropsfileList[count],IPDfileList[count]);
					if(vals!=7)
					{
						printf("Source file: reading %d values, but expected 7 values!\n",vals);
						exit(-1);
					}
				count++;
				line=fgets(line,LINESIZE,fp);
			}
			else
			{
				line=fgets(line,LINESIZE,fp);
			}
		}
	}

	if(count!=nSyst)
	{
		printf("Source Error: expected %d lines but found %d lines\n",nSyst,count);
		exit(-1);
	}

	free(line);
}

void readChromPos(struct PARAINPUT *p,struct SYSTEM *box)
{
	int i,ii,j,jj;
	int n,spamtype;
	char datatostr[300];
	double spamrad,spamtrans,spamx,spamy,spamz;
	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	// this subroutine reads in positions of chromosomes, and surf points from a dump file
	// only sanity checks is that the number of chromosomes and the surf points have to be as we specified for this restart instance.
	// simulation subsequent to this will allow further relaxation, possibly with different parameters and placement of patches etc

	fscanf(box->fpRestart,"%s",datatostr);
	n=atoi(datatostr);
	if(n != p->nDChrom+p->nSurfPoints)
	{
		printf("mismatch in number of particles and surf points, %d %d %d\n",n,p->nDChrom,p->nSurfPoints);
	}
	fscanf(box->fpRestart,"%s",datatostr);

	// now read surf points
	for(i=0;i<p->nSurfPoints;i++)
	{
		// surfid, type, transparency, radius, x,y,z
		fscanf(box->fpRestart,"%d %d %lf %lf %lf %lf %lf\n",&n,&spamtype,&spamtrans,&spamrad,&box->surfPointx[i],&box->surfPointy[i],&box->surfPointz[i]);
		if(n!=i)
		{
			printf("surf point index mismatch %d %d\n",i,n);
			exit(0);
		}

		// cross-check raddi
		if(fabs(0.5*avesep -spamrad) >0.0001)
		{
			printf("radius issue for surfpoint %d %f should be %f -- check geom and density etc\n",n,spamrad,0.5*avesep);
			exit(0);
		}
	}
	// now read chromosome
	for(ii=0;ii<p->nDChrom;ii++)
	{
		fscanf(box->fpRestart,"%d %d %lf %lf %lf %lf %lf\n",&n,&spamtype,&spamtrans,&spamrad,&spamx,&spamy,&spamz);
		if(ii!=n || n >= p->nDChrom)
		{
			printf("mismatch chrom index %d %d\n",i,n);
			exit(0);
		}

		// check id mismatch
		i = box->chromList[ii]; // global index for this chrom
		if(box->whichChrom[i] != n)
		{
			printf("chromid mismatch %d %d glob %d\n",ii,n,i);
			exit(0);
		}

		// check type mismatch
		if(spamtype != box->chromType[n])
		{
			printf("chrom type mismatch for %d - %d %d\n",n,spamtype,box->chromType[n]);
			exit(0);
		}
		// store xyz
		box->X[i] = spamx;
		box->Y[i] = spamy;
		box->Z[i] = spamz;


	}

}

void dropSimStart(FILE *fpWrite,struct PARAINPUT *p,struct SYSTEM *box)
{

	int ii,i,j,jj;
	fprintf(fpWrite,"%lf\n",8.0*p->ella);
	double mspam = 1.0;
	for(ii=0;ii<p->nDChrom;ii++)
	{
		// first write the info for this chrom
		//fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&chromindex,&type,&spamx,&spamy,&spamz,&orix,&oriy,&oriz,&rad,&mass,&k_spring,&act,&lamin);
		i=box->chromList[ii];
// 		printf("wrote %f, writing %f %d %d\n",p->ella,box->X[i],ii,i);
		// identify rotation angle for the chrom

			// chrom type x y z a(o)x a(o)y a(o)z r m k_spring f_act e_lamin

		fprintf(fpWrite,"%d %d %f %f %f %f %f %f %f %f %f %f %f\n",ii+1,box->chromType[ii]+1,box->X[i],box->Y[i],box->Z[i],box->ox[i],box->oy[i],box->oz[i],mspam*box->radius[i],box->mass[i],0.0,0.0,0.0);
// 		printf("printing chrom info no issue %d -- %d\n",ii,box->patchList[ii][0]);
		for(j=1;j<=box->patchList[ii][0];j++)
		{
			jj=box->patchList[ii][j];
// 			printf("chrom %d raw %d patch %d raw %d\n",ii,i,j,jj);
			if(box->lamin[ii][j]==0)
			fprintf(fpWrite,"%d %d %f %f %f %f %f %f %f %f %f %f %f\n",ii+1,0,box->patchposx[ii][j],box->patchposy[ii][j],box->patchposz[ii][j],box->actox[ii][j],box->actoy[ii][j],box->actoz[ii][j],mspam*box->radius[jj],box->mass[jj],p->k_spring_level,p->f_act_level,0.0);
			if(box->lamin[ii][j]==1)
			{
// 				printf("lamin patch for chrom %d raw %d -- patch %d of %d raw %d\n",ii,i,j,box->patchList[ii][0],jj);
				fprintf(fpWrite,"%d %d %f %f %f %f %f %f %f %f %f %f %f\n",ii+1,0,box->patchposx[ii][j],box->patchposy[ii][j],box->patchposz[ii][j],box->actox[ii][j],box->actoy[ii][j],box->actoz[ii][j],mspam*box->radius[jj],box->mass[jj],p->k_spring_level,0.0,p->lamin_scale);
			}
		}
	}
}


void writeSeparations(FILE *fpWrite,struct PARAINPUT *p, struct SYSTEM *box)
{
	int ii,jj,i,j;
	double rX,rY,rZ,rij;


	for(ii=0;ii<p->nDChrom;ii++)
	{
		i = box->chromList[ii];

		for(jj=0;jj<p->nDChrom;jj++)
		{
			j = box->chromList[jj];

			rX = box->X[j] - box->X[i];
			rY = box->Y[j] - box->Y[i];
			rZ = box->Z[j] - box->Z[i];

			rij = sqrt(rX*rX + rY*rY + rZ*rZ);

			fprintf(fpWrite,"%d %d %f %f %f\n",ii+1,jj+1,box->radius[i],box->radius[j],rij);
		}
	}
}
