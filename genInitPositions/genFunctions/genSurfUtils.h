// double gaussVar(double mu,double sig)
// {
//     double gauss_r,v1,v2,gauss_var;
//     gauss_r=2.0;
//     while(gauss_r>=1.0)
//     {
// //         v1=2.0*genrand_real1(&(*a))-1.0;
// //         v2=2.0*genrand_real1(&(*a))-1.0;
// 		v1 = 2.0*drand48()-1.0;
// 		v2 = 2.0*drand48()-1.0;
// //         printf("v1 v2 %f %f\n",v1,v2);
//         gauss_r = v1*v1 + v2*v2;
//     }
// //     printf("gr log gr %f %f\n",gauss_r,log(gauss_r));
//     gauss_var=v1*sqrt(-2.0*log(gauss_r)/gauss_r);
//     gauss_var=mu+sig*gauss_var;
//
//     return gauss_var;
// }


void genSurf(struct PARAINPUT *p, struct SYSTEM *box)
{
	// subroutine generates uniformly distributed points on surface of ellipsoid
	// main parameters are number of points, and packing fraction on surface
	// take the ellipsoid semi axes a>=b>=c (a=b=c implies sphere)

	int i,j;
	double curr[4];
	double rX,rY,rZ,r,pkeep,random;
	int atomdone,atomnotdone;

	double spam = 0.33333*(pow(p->ella*p->ellb,1.6) + pow(p->ella*p->ellc,1.6) + pow(p->ellb*p->ellc,1.6));
	double surf = 4.0*M_PI*pow(spam,1.0/1.6); // surface area of ellipsoid
	printf("fact %f, %f - fraction %f\n",spam,surf,p->surfPhi);
	double avesep = 2.0*pow(p->surfPhi*4.0*surf/(M_PI*p->nSurfPoints),0.5); // N*pi*r^2/surf = phi (sep = 2*r)
	printf("placing %d surface points %f distance apart\n",p->nSurfPoints,avesep);
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
				if(r<0.25*avesep*avesep) // too close
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
			box->surfPointx[i] = p->ella*curr[0]; // keeps getting replaced at iterations
			box->surfPointy[i] = p->ellb*curr[1]; // store transformed ellipsoidal positions
			box->surfPointz[i] = p->ellc*curr[2];



		}
// 		printf("surf point %d at %f %f %f\n",i,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i]);
// 		printf("%d %d surf %f %f %f -- fund %f %f %f\n",i,p->nSurfPoints,box->surfPointx[i],box->surfPointy[i],box->surfPointz[i],box->surfFundx[i],box->surfFundy[i],box->surfFundz[i]);
	}



	printf("done making surf\n");

}
