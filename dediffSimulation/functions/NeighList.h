long int lround(double x);

int vnlist(struct PARAINPUT *p,struct SYSTEM *box)
{
	double delx,dely,delz,dr2;
	double xi,yi,zi;

// 	int *llist,head[10][10][10],nbin,**mycell;
    int nbin;
	double binsz,rskinmax,lxh;
	int ix,iy,iz;
	int ixn,iyn,izn;
	int i,j,neigcount_i,neigcount_j,lcount=0,k;
	int ii,jj,kk;
    int ch;
// 	llist=calloc(p->nAto,sizeof(int));
// 	mycell=calloc(p->nAto,sizeof(int *));

// 	for(i=0;i<p->nAto;i++)
// 	{
// 		mycell[i]=calloc(4,sizeof(int));
// 	}
    
    for(i=0;i<p->nAto;i++)
    {
        box->NEIGLIST[i][0]=0;
    }
	/* Prima definisco le linked lists. Poi le uso per la verlet neighborlist */

	lxh=box->lx/2.0;
// 	nbin=(int) (box->lx/(p->rc[0]+SKIN));
	nbin=(int) (box->lx/(p->l+SKIN));
	binsz = box->lx/nbin; 
//     printf("%d %f %f %d %f\n",nbin,binsz,lxh,ch,p->l);
  
	for(i=0;i<p->nAto;i++)
	{
		xi=box->X[i];
		yi=box->Y[i];
		zi=box->Z[i];

        //neigcount_i=0;
		neigcount_i=box->NEIGLIST[i][0]; 	// Trovo i vicini della particella i
// 		if(i==2 || i==213)
//         printf("\nneigcount_i %d %d \n",i,neigcount_i);
		for(j=0;j<p->nAto;j++)
        {
            if(j>i)
            {
                ch=box->type[i]+box->type[j];
                delx=xi-box->X[j];
//                 delx-=box->lx*lround(delx/box->lx);
                dely=yi-box->Y[j];
//                 dely-=box->ly*lround(dely/box->ly);
                delz=zi-box->Z[j];
//                 dr2=delx*delx+dely*dely+delz*delz;
                dr2=delx*delx+dely*dely+delz*delz;

    // 					if(dr2<=sqr(0.5*box->lx))
//                 if(dr2<=sqr(p->rc[ch]+SKIN))
				if(dr2<=sqr(2.0*p->l+2.0*SKIN))
                {
					
//                     if(i==2 || i==213)
//                     {
//                         printf("full ij %d %d \n",i,j);
//                     }
//                     if(j==2 || j==213)
//                     {
//                         printf("full ji %d %d \n",j,i);
//                     }
                    neigcount_i++; // Aggiorno i vicini di i
//                     printf("%d %d\n",neigcount_i,j);
					if(neigcount_i>box->nlen)
						printf("particle %d exceeded neighbour count %d\n",i,neigcount_i);
                    box->NEIGLIST[i][neigcount_i]=j;
                    neigcount_j=box->NEIGLIST[j][0]+1;	// Aggiorno i vicini di j
                    box->NEIGLIST[j][neigcount_j]=i;
                    // printf("%d %d\n",neigcount_j,i);
                    box->NEIGLIST[j][0]=neigcount_j;
                    //printf("atom %d has a neig %d\n",i,j);

                }	
            }
        }
        

		box->NEIGLIST[i][0]=neigcount_i;
//         if(i==2 || i==213)
       // printf("full calc %d %d %f\n",i,neigcount_i,box->ori[i]);
        //for(j=1;j<=neigcount_i;j++)
        //{
         //   printf(" %d",box->NEIGLIST[i][j]);
        //}
        //printf("\n");
		lcount=MAX(lcount,neigcount_i);

		box->dispX[i]=0.0;
		box->dispY[i]=0.0;
		box->dispZ[i]=0.0;
	}
// 	free(llist);
// 	for(i=0;i<p->nAto;i++)
// 	{
// 		free(mycell[i]);
// 	}
// 	free(mycell);
    //printf("finished nlist\n");
	return(nbin);
}

void vnlistpar(int par,struct SYSTEM *box,double delx,double dely,double delz,struct PARAINPUT *p)
{
	int j,lcount=0,neigcount;
	double rijsq,Xij,Yij,Zij;
	double displacesq, neipar;
	int ch;
	//displacesq=sqr(box->dispX[par]+delx)+sqr(box->dispY[par]+dely)+sqr(box->dispZ[par]+delz);
    displacesq=sqr(box->dispX[par]+delx)+sqr(box->dispY[par]+dely);
//     if(par==2 || par==213)
//     printf("displace %d %f\n",par,displacesq);
    box->NEIGPAR[0]=0;
	if(displacesq>sqr(0.5*SKIN))
	{
		for(j=0;j<p->nAto;j++)
		{
			if(j!=par)
			{
                ch=box->type[par]+box->type[j];
				Xij=box->X[j]-box->X[par];
// 				Xij-=box->lx*lround(Xij/box->lx);
				Yij=box->Y[j]-box->Y[par];
// 				Yij-=box->ly*lround(Yij/box->ly);
				Zij=box->Z[j]-box->Z[par];
// 				Zij-=box->lz*lround(Zij/box->lz);

				rijsq=Xij*Xij+Yij*Yij+Zij*Zij;
// 				rijsq=Xij*Xij+Yij*Yij;
                
				if(rijsq<=sqr(2.0*p->l+2.0*SKIN))
				{   

					lcount++;
					box->NEIGPAR[lcount]=j;
				}
			}
		}
		box->NEIGPAR[0]=lcount;
		// update full neiglist here??
	}
	else
	{
		neigcount=box->NEIGLIST[par][0];
		for(j=1;j<=neigcount;j++)
		{
			box->NEIGPAR[j]=box->NEIGLIST[par][j];
		}
		box->NEIGPAR[0]=neigcount;
	}
}

