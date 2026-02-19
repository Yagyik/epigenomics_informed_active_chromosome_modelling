// this file contains subroutine related to rigid body manipulations
// first, easiest function is to find positions of all patches in the chromosome frame
// given ox,oy,oz which should be (1,0,0), find the normal vector (ox,oy,oz) x (1,0,0)
// find the angle cos(theta) = (ox,oy,oz) . (1,0,0)
// second subroutine rotates a vector about axis by the chosen angle

void findQDotG(double oldqx,double oldqy,double oldqz,double oldqw,double wx,double wy,double wz,double *qdot)
{
    // subroutine to find the time derivative of quaternion from angular velocities in global reference frame
//     double Wmat[16]; // used if we represent old q in a matrix for and call a quaternion multiplication function
    double ww=0.0;
//     qdot[0] = oldqw*wx - oldqx*wy - oldqy*wz - oldqz*0.0;
//     qdot[1] = oldqx*wx + oldqw*wy + oldqz*wz - oldqy*0.0;
//     qdot[2] = oldqy*wx - oldqz*wy + oldqw*wz + oldqx*0.0;
//     qdot[3] = oldqz*wx + oldqy*wy - oldqx*wz + oldqw*0.0;
	
	// NOTE: global frame angular velocity from global frame torque qdot = 1/2 w o q (where 'o' denotes quaternion product, q = (qw,qx,qy,qz) and w = (0,wx,wy,wz) is in global frame)
	// NOTE: because of the convention here (qw,qx,qy,qz) we have to re-order our qdot accordingly
	qdot[3] = 0.5*(ww*oldqw - wx*oldqx - wy*oldqy - wz*oldqz);
	qdot[0] = 0.5*(wx*oldqw + ww*oldqx - wz*oldqy + wz*oldqz);
	qdot[1] = 0.5*(wy*oldqw + wz*oldqx + ww*oldqy - wx*oldqz);
	qdot[2] = 0.5*(wz*oldqw - wy*oldqx + wx*oldqy + ww*oldqz);
//     printf("in calc qdot do we ruin oldq %f %f %f %f\n",oldqx,oldqy,oldqz,oldqw);
    // and we're done!
}

double findAxesQuat(double qx,double qy,double qz,double qw,double *axis)
{
    double angle,qnorm,norm,halftheta;
    halftheta = acos(qw);
// 	printf("qw %f cosinv %f and new cos angle %f\n",qw,acos(qw),cos(2.0*halftheta));
    angle = cos(2.0*halftheta);
    qnorm = qx*qx + qy*qy + qz*qz + qw*qw;
    axis[0] = qx/(sin(halftheta)*sqrt(qnorm));
    axis[1] = qy/(sin(halftheta)*sqrt(qnorm));
    axis[2] = qz/(sin(halftheta)*sqrt(qnorm));
    
    norm = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
    
    axis[0] /= sqrt(norm);
    axis[1] /= sqrt(norm);
    axis[2] /= sqrt(norm);
	
	return angle;
}

void findAxesInitOri(struct PARAINPUT *p,struct SYSTEM *box) // this subroutine only called in the beginning to rotate the chromosome by some amount
{
	// find the axis of each chromosome
	int ii,i,j;
	int ref_lamin;
	double norm1,norm2,norm3;
	double normprod;
	for(i=0;i<p->nChrom;i++)
	{
		ref_lamin = box->patchList[i][0]; // last patch lamin patch
		printf("\nchromosome %d of %d, real %d -- lamin %d raw %d\n",i,p->nChrom,box->chromList[i],box->patchList[i][0],box->patchList[i][box->patchList[i][0]]);
		// cross product between (-1,0,0),i.e., patch[i][0] and (ox,oy,oz)
		ii = box->chromList[i];

		printf("chrom %d is atom %d -- atom %d is reference at body position %f %f %f\n",i,ii,box->patchList[i][ref_lamin],box->patchposx[i][ref_lamin],box->patchposy[i][ref_lamin],box->patchposz[i][ref_lamin]);
		printf("chrom %d %d initial axes %f %f %f\n",i,ii,box->ox[ii],box->oy[ii],box->oz[ii]);
		norm1 = sqrt(box->patchposx[i][ref_lamin]*box->patchposx[i][ref_lamin] + box->patchposy[i][ref_lamin]*box->patchposy[i][ref_lamin] + box->patchposz[i][ref_lamin]*box->patchposz[i][ref_lamin]);


		norm2 = sqrt(box->ox[ii]*box->ox[ii] + box->oy[ii]*box->oy[ii] + box->oz[ii]*box->oz[ii]);

		if(norm2 == 0)
		{
			printf("initial orientation vector has zero norm -- in the default case should at least be negative of lamin vector, i.e., patch %d on chrom %d\n",ref_lamin,i);
			exit(0);
		}
		normprod = norm1*norm2;
		box->axisx[ii] = box->oy[ii]*box->patchposz[i][ref_lamin]/normprod - box->oz[ii]*box->patchposy[i][ref_lamin]/normprod;
		box->axisy[ii] = box->oz[ii]*box->patchposx[i][ref_lamin]/normprod - box->ox[ii]*box->patchposz[i][ref_lamin]/normprod;
		box->axisz[ii] = box->ox[ii]*box->patchposy[i][ref_lamin]/normprod - box->oy[ii]*box->patchposx[i][ref_lamin]/normprod;

		if(box->axisx[ii]*box->axisx[ii] + box->axisy[ii]*box->axisy[ii] + box->axisz[ii]*box->axisz[ii] <= 0.0)
		{
			printf("initial ori parallel to reference patch axis (cross product coming zero), rotating by delta for convenience\n");
			box->ox[ii] = box->ox[ii] + 0.1*(2.0*drand48()-1.0);
			box->oy[ii] = box->oy[ii] + 0.1*(2.0*drand48()-1.0);
			box->oz[ii] = box->oz[ii] + 0.1*(2.0*drand48()-1.0);
			printf("new orientations %d %d %f %f %f \n",i,ii,box->ox[ii],box->oy[ii],box->oz[ii]);
			norm2 = sqrt(box->ox[ii]*box->ox[ii] + box->oy[ii]*box->oy[ii] + box->oz[ii]*box->oz[ii]);

			normprod = norm1*norm2;
			printf("check norm chrom %d raw %d -- %f %f\n",i,ii,norm1,norm2);
			box->axisx[ii] = box->oy[ii]*box->patchposz[i][ref_lamin]/normprod - box->oz[ii]*box->patchposy[i][ref_lamin]/normprod;
			box->axisy[ii] = box->oz[ii]*box->patchposx[i][ref_lamin]/normprod - box->ox[ii]*box->patchposz[i][ref_lamin]/normprod;
			box->axisz[ii] = box->ox[ii]*box->patchposy[i][ref_lamin]/normprod - box->oy[ii]*box->patchposx[i][ref_lamin]/normprod;
		}



		norm3 = sqrt(box->axisx[ii]*box->axisx[ii] + box->axisy[ii]*box->axisy[ii] + box->axisz[ii]*box->axisz[ii]);


		printf("norm of ori is (%g,%g,%g) %f and patch is (%g,%g,%g) %f -- to axis %f\n",box->ox[ii],box->oy[ii],box->oz[ii],norm2,box->patchposx[i][ref_lamin],box->patchposy[i][ref_lamin],box->patchposz[i][ref_lamin],norm1,norm3);
		
		box->axisx[ii] /=norm3;
		box->axisy[ii] /=norm3;
		box->axisz[ii] /=norm3;

		
		// dot product divided by norm gives cos(theta)
		
		box->cosRotAngle[ii] = box->patchposx[i][ref_lamin]*box->ox[ii] + box->patchposy[i][ref_lamin]*box->oy[ii] + box->patchposz[i][ref_lamin]*box->oz[ii];
		
		// norm of both vectors should already be 1 but doesn't hurt to check
		box->cosRotAngle[ii] /= norm1*norm2;


		printf("final axes for chrom %d real ato %d -- axis %f %f %f, angle %f\n",i,ii,box->axisx[ii],box->axisy[ii],box->axisz[ii],box->cosRotAngle[ii]);
	}
	
}


void rotate(double rx,double ry, double rz, double cx,double cy,double cz,double costheta,double *rp)
{
	double dotprod,crossprodx,crossprody,crossprodz;
	double normr,normc;
	
	// use axis and angle to rotate a given vector
	// axis needs to be normalised, resultant vector will have same norm as original vector
	// rp = r*cos(theta) + (c.r)*c*(1-cos(theta)) + (cxr)*sin(theta)
	dotprod = cx*rx + cy*ry + cz*rz;
	crossprodx = cy*rz - cz*ry;
	crossprody = cz*rx - cx*rz;
	crossprodz = cx*ry - cy*rx;
	
	rp[0] = rx*costheta  + dotprod*cx*(1-costheta) + crossprodx*sqrt(1-costheta*costheta); // taking from body coords to global coords based on axis angle of chrom
	rp[1] = ry*costheta  + dotprod*cy*(1-costheta) + crossprody*sqrt(1-costheta*costheta); // same
	rp[2] = rz*costheta  + dotprod*cz*(1-costheta) + crossprodz*sqrt(1-costheta*costheta); // same
	

// 	rp[0] = rx;
// 	rp[1] = ry;
// 	rp[2] = rz;
// 	printf("norm of (%g,%g,%g) %f -- rotated thru (%g,%g,%g) %f is %f\n",rx,ry,rz,sqrt(rx*rx+ry*ry+rz*rz),cx,cy,cz,sqrt(cx*cx+cy*cy+cz*cz),sqrt(rp[0]*rp[0]+rp[1]*rp[1]+rp[2]*rp[2]));
}





// need to find subroutine for Euler angles for a given vector


// need to find subroutine to represent quaternions for given euler angle

// find axis from quaternions

// etc
