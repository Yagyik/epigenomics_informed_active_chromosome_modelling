void countList(struct PARAINPUT *p, double **source, int size,int rank)
{
	int i,j,k,n;
	double eps,T;


	T=source[0][0];
    eps=source[0][1];

	p->nT=1;
	p->neps=1;
	printf("setting up counts\n");
	for(i=1;i<size;i++)
	{
		if(eps==source[i][1] && T!= source[i][0])
		{
			p->nT++;
		}

	}
	for(i=1;i<size;i++)
	{
		if(T==source[i][0] && eps!= source[i][1])
		{
			p->neps++;
		}
	}
	printf("accumu neps, nt\n");
	int spamcount = 1;
	for(i=1;i<size;i++)
	{
		if(T==source[i][0] && eps == source[i][1])
		{
			spamcount++;
		}
	}
	printf("%d %d %d\n",p->nT,p->neps,spamcount);


	printf("size spam  neps, nT %d %d %d %d\n",size,spamcount,p->neps,p->nT);
	if(p->neps*p->nT!=size && spamcount !=size)
	{
		printf("PT warning: grid of %d values, but using %d cores!!\n",p->neps*p->nT,size);
	}

	if(spamcount % 4 == 0)
	{
		p->neps = 4;
		p->nT = spamcount/4;
	}
	else if(spamcount % 3 == 0)
	{
		p->neps = 3;
		p->nT = spamcount/3;
	}
	else if(spamcount == 1)
	{
		p->neps = 1;
		p->nT = 1;
	}
	else
	{
		printf("number of procs needs to be a multiple of 3 or 4 -- curr value %d\n",size);
		exit(0);
	}
	printf("FINAL size neps, nT %d %d %d\n",size,p->neps,p->nT);
// 	exit(0);

	p->procList=calloc(p->neps,sizeof(int *));
	for(i=0;i<p->neps;i++)
	{
		p->procList[i]=calloc(p->nT,sizeof(int *));
	}


	n=0;
	printf("alloc procList\n");
	for(i=0;i<p->neps;i++)
    {
	for(j=0;j<p->nT;j++)
	{
		p->procList[i][j]=n;		// This order follows the makeSource.sh's one!!!
		if(rank==0)
		{
		printf("rank %d i %d %d - n %d\n",rank,i,j,n);
		printf(" source %d , %d i eps %f , %d j T %f\n",n,i,source[n][1],j,source[n][0]);
		}
		n++;


	}
	}

}

void TswapSend(struct SYSTEM *box, int proc)
{
                                        // tag 0 = Volume
	MPI_Send(&box->curren,1,MPI_DOUBLE,proc,0,MPI_COMM_WORLD);							// tag 0 = U
	MPI_Send(&box->currW,1,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);									// tag 1 = W
	MPI_Recv(&box->yestoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);			// tag 2 = esitoSwap
}

void TswapReceive(struct SYSTEM *box, int proc, double **source,struct PARAINPUT *p)
{
	double Uswap,Wswap,ran;
	MPI_Recv(&Uswap,1,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 0 = U
	MPI_Recv(&Wswap,1,MPI_DOUBLE,proc,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 1 = W


	box->deltaBeta=(1.0/p->mctemp-1.0/source[proc][0]); // find out what this does
	ran=drand48();	// Attento: questa formula e' vera se le finestre sono esattamente uguali come k,N0!


	if(ran<exp(box->deltaBeta*(box->curren-Uswap)))
	{
		box->yestoSwap=1;
		MPI_Send(&box->yestoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
	else
	{
		box->yestoSwap=0;
		MPI_Send(&box->yestoSwap,1,MPI_INT,proc,2,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
}

void epsSwapSend(struct SYSTEM *box, struct PARAINPUT *p, int proc)
{



	MPI_Send(&p->eps_hertz,1,MPI_DOUBLE,proc,0,MPI_COMM_WORLD);									// tag 0=eps


    // now send x,y,z,surfx,surfy,surfz
    MPI_Send(box->X,p->nAto,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
    MPI_Send(box->Y,p->nAto,MPI_DOUBLE,proc,2,MPI_COMM_WORLD);
    MPI_Send(box->Z,p->nAto,MPI_DOUBLE,proc,3,MPI_COMM_WORLD);

    MPI_Send(box->surfPointx,p->nSurfPoints,MPI_DOUBLE,proc,4,MPI_COMM_WORLD);
    MPI_Send(box->surfPointy,p->nSurfPoints,MPI_DOUBLE,proc,5,MPI_COMM_WORLD);
    MPI_Send(box->surfPointz,p->nSurfPoints,MPI_DOUBLE,proc,6,MPI_COMM_WORLD);

	MPI_Recv(&box->yestoSwap,1,MPI_INT,proc,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE);			// tag 2 = esitoSwap

}

void epsSwapRiceive(struct SYSTEM *box, int proc, double **source,struct PARAINPUT *p)
{
	double ran,ennew,enold;
    double epsrecv;
    double Xrecv[p->nAto];
    double Yrecv[p->nAto];
    double Zrecv[p->nAto];


    double surfXrecv[p->nSurfPoints];
    double surfYrecv[p->nSurfPoints];
    double surfZrecv[p->nSurfPoints];


	MPI_Recv(&epsrecv,1,MPI_DOUBLE,proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag0=eps
    // now receive the x,y,z,surfx,surfy,surfz
    MPI_Recv(&Xrecv,p->nAto,MPI_DOUBLE,proc,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&Yrecv,p->nAto,MPI_DOUBLE,proc,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&Zrecv,p->nAto,MPI_DOUBLE,proc,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    MPI_Recv(&surfXrecv,p->nSurfPoints,MPI_DOUBLE,proc,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&surfYrecv,p->nSurfPoints,MPI_DOUBLE,proc,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&surfZrecv,p->nSurfPoints,MPI_DOUBLE,proc,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE);


    // enold is curr eps with curr pos (enDiffEps with box var and p->eps) + enDiffEps with received box and received eps)
    // ennew is curr eps with diff pos (enDirffEps with box var and received eps) + enDiffEps with received box and p->eps)

    enold = enDiffEps(p,box,p->eps_hertz,box->X,box->Y,box->Z,box->surfPointx,box->surfPointy,box->surfPointz);
    enold += enDiffEps(p,box,epsrecv,&(*Xrecv),&(*Yrecv),&(*Zrecv),&(*surfXrecv),&(*surfYrecv),&(*surfZrecv));

    ennew = enDiffEps(p,box,epsrecv,box->X,box->Y,box->Z,box->surfPointx,box->surfPointy,box->surfPointz);
    ennew += enDiffEps(p,box,p->eps_hertz,&(*Xrecv),&(*Yrecv),&(*Zrecv),&(*surfXrecv),&(*surfYrecv),&(*surfZrecv));

	ran=drand48();

	if(ran<exp(-(ennew-enold)/p->mctemp))
	{
		box->yestoSwap=1;
		MPI_Send(&box->yestoSwap,1,MPI_INT,proc,7,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}
	else
	{
		box->yestoSwap=0;
		MPI_Send(&box->yestoSwap,1,MPI_INT,proc,7,MPI_COMM_WORLD);		// tag 2 = esitoSwap
	}

}

void SendTmpBox(struct SYSTEM *box,struct SYSTEM *tmpbox,struct PARAINPUT *p,int proc)
{
	int i;


	for(i=0;i<p->nAto;i++)
	{
		tmpbox->X[i]=box->X[i];
		tmpbox->Y[i]=box->Y[i];
		tmpbox->Z[i]=box->Z[i];


	}
	for(i=0;i<p->nSurfPoints;i++)
	{

		tmpbox->surfPointx[i]=box->surfPointx[i];
		tmpbox->surfPointy[i]=box->surfPointy[i];
		tmpbox->surfPointz[i]=box->surfPointz[i];

	}

	MPI_Send(tmpbox->X,p->nAto,MPI_DOUBLE,proc,13,MPI_COMM_WORLD);		// tag 13 = tmpbox->X
	MPI_Send(tmpbox->Y,p->nAto,MPI_DOUBLE,proc,14,MPI_COMM_WORLD);		// tag 14 = tmpbox->Y
	MPI_Send(tmpbox->Z,p->nAto,MPI_DOUBLE,proc,15,MPI_COMM_WORLD);		// tag 15 = tmpbox->Z


	MPI_Send(tmpbox->surfPointx,p->nSurfPoints,MPI_DOUBLE,proc,16,MPI_COMM_WORLD);		// tag 16 = tmpbox->dispX
	MPI_Send(tmpbox->surfPointy,p->nSurfPoints,MPI_DOUBLE,proc,17,MPI_COMM_WORLD);		// tag 17 = tmpbox->dispY
	MPI_Send(tmpbox->surfPointz,p->nSurfPoints,MPI_DOUBLE,proc,18,MPI_COMM_WORLD);		// tag 18 = tmpbox->dispZ

    // Sent all the info on tmpbox. Now I receive in box

	MPI_Recv(box->X,p->nAto,MPI_DOUBLE,proc,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 23 = box->X
	MPI_Recv(box->Y,p->nAto,MPI_DOUBLE,proc,24,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 24 = box->Y
	MPI_Recv(box->Z,p->nAto,MPI_DOUBLE,proc,25,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 25 = box->Z


	MPI_Recv(box->surfPointx,p->nSurfPoints,MPI_DOUBLE,proc,26,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 26 = box->dispX
	MPI_Recv(box->surfPointy,p->nSurfPoints,MPI_DOUBLE,proc,27,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 27 = box->dispY
	MPI_Recv(box->surfPointz,p->nSurfPoints,MPI_DOUBLE,proc,28,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 28 = box->dispZ

}

void ReceiveTmpBox(struct SYSTEM *box,struct SYSTEM *tmpbox,struct PARAINPUT *p,int proc)
{
	int i;


	for(i=0;i<p->nAto;i++)
	{
		tmpbox->X[i]=box->X[i];
		tmpbox->Y[i]=box->Y[i];
		tmpbox->Z[i]=box->Z[i];


	}
	for(i=0;i<p->nSurfPoints;i++)
	{

		tmpbox->surfPointx[i]=box->surfPointx[i];
		tmpbox->surfPointy[i]=box->surfPointy[i];
		tmpbox->surfPointz[i]=box->surfPointz[i];

	}

	MPI_Recv(box->X,p->nAto,MPI_DOUBLE,proc,13,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 13 = tmpbox->X
	MPI_Recv(box->Y,p->nAto,MPI_DOUBLE,proc,14,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 14 = tmpbox->Y
	MPI_Recv(box->Z,p->nAto,MPI_DOUBLE,proc,15,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 15 = tmpbox->Z


	MPI_Recv(box->surfPointx,p->nSurfPoints,MPI_DOUBLE,proc,16,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 10 = box->dispX
	MPI_Recv(box->surfPointy,p->nSurfPoints,MPI_DOUBLE,proc,17,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 11 = box->dispY
	MPI_Recv(box->surfPointz,p->nSurfPoints,MPI_DOUBLE,proc,18,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		// tag 12 = box->dispZ


	MPI_Send(tmpbox->X,p->nAto,MPI_DOUBLE,proc,23,MPI_COMM_WORLD);		// tag 13 = box->X
	MPI_Send(tmpbox->Y,p->nAto,MPI_DOUBLE,proc,24,MPI_COMM_WORLD);		// tag 14 = box->Y
	MPI_Send(tmpbox->Z,p->nAto,MPI_DOUBLE,proc,25,MPI_COMM_WORLD);		// tag 15 = box->Z


	MPI_Send(tmpbox->surfPointx,p->nSurfPoints,MPI_DOUBLE,proc,26,MPI_COMM_WORLD);		// tag 20 = box->dispX
	MPI_Send(tmpbox->surfPointy,p->nSurfPoints,MPI_DOUBLE,proc,27,MPI_COMM_WORLD);		// tag 21 = box->dispY
	MPI_Send(tmpbox->surfPointz,p->nSurfPoints,MPI_DOUBLE,proc,28,MPI_COMM_WORLD);		// tag 22 = box->dispZ


}
