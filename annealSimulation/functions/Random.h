/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.						  

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

	 1. Redistributions of source code must retain the above copyright
		notice, this list of conditions and the following disclaimer.

	 2. Redistributions in binary form must reproduce the above copyright
		notice, this list of conditions and the following disclaimer in the
		documentation and/or other materials provided with the distribution.

	 3. The names of its contributors may not be used to endorse or promote 
		products derived from this software without specific prior written 
		permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/* Period parameters */  
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s,struct APPOGGIO *a)
{
	a->mt[0]= s & 0xffffffffUL;
	
	for(a->mti=1;a->mti<N;a->mti++)
	{	
// 		printf("%d is ok\n",a->mti);
		a->mt[a->mti]=(1812433253UL*(a->mt[a->mti-1]^(a->mt[a->mti-1]>>30))+a->mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].						*/
		/* 2002/01/09 modified by Makoto Matsumoto			 */
		a->mt[a->mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
	printf("completed genrand\n");
}

unsigned long genrand_int32(struct APPOGGIO *a)
{
	unsigned long y;
	static unsigned long mag01[2]={0x0UL, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if(a->mti>=N)
	{
		/* generate N words at one time */
		int kk;

		if(a->mti==N+1)		/* if init_genrand() has not been called, */
		{
			init_genrand(5489UL,&(*a));	/* a default initial seed is used */
		}

		for(kk=0;kk<N-M;kk++)
		{
			y=(a->mt[kk]&UPPER_MASK)|(a->mt[kk+1]&LOWER_MASK);
			a->mt[kk]=a->mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for(;kk<N-1;kk++)
		{
			y=(a->mt[kk]&UPPER_MASK)|(a->mt[kk+1]&LOWER_MASK);
			a->mt[kk]=a->mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y=(a->mt[N-1]&UPPER_MASK)|(a->mt[0]&LOWER_MASK);
		a->mt[N-1]=a->mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
		
		a->mti=0;
	}

	y=a->mt[a->mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}


/* generates a random number on [0,1]-real-interval */
double genrand_real1(struct APPOGGIO *a)
{
	return genrand_int32(&(*a))*(1.0/4294967295.0); 
	/* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(struct APPOGGIO *a)
{
	return genrand_int32(&(*a))*(1.0/4294967296.0); 
	/* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(struct APPOGGIO *a)
{
	return (((double)genrand_int32(&(*a))) + 0.5)*(1.0/4294967296.0); 
	/* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(struct APPOGGIO *app)
{ 
	unsigned long a=genrand_int32(&(*app))>>5, b=genrand_int32(&(*app))>>6; 
	return(a*67108864.0+b)*(1.0/9007199254740992.0); 
}
