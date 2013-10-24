/*
 *  ranmar.c
 *  kiko
 *
 *  Created by Wilson on 29/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "ranmar.h"

#define TRUE -1
#define FALSE 0
#define boolean int

static float u[98], c, cd, cm;
static int i97, j97;
static boolean test = FALSE;

void rmarin(int ij,int kl)
{
	/*
	 C This is the initialization routine for the random number generator RANMAR()
	 C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
	 C                                                      0 <= KL <= 30081
	 C The random number sequences created by these two seeds are of sufficient 
	 C length to complete an entire calculation with. For example, if sveral 
	 C different groups are working on different parts of the same calculation,
	 C each group could be assigned its own IJ seed. This would leave each group
	 C with 30000 choices for the second seed. That is to say, this random 
	 C number generator can create 900 million different subsequences -- with 
	 C each subsequence having a length of approximately 10^30.
	 C 
	 C Use IJ = 1802 & KL = 9373 to test the random number generator. The
	 C subroutine RANMAR should be used to generate 20000 random numbers.
	 C Then display the next six random numbers generated multiplied by 4096*4096
	 C If the random number generator is working properly, the random numbers
	 C should be:
	 C           6533892.0  14220222.0  7275067.0
	 C           6172232.0  8354498.0   10633180.0
	 */
	int i, j, k, l, ii, jj, m;
	float s, t;
	
	if (ij<0 || ij>31328 || kl<0 || kl>30081) {
		puts("The first random number seed must have a value between 0 and 31328.");
		puts("The second seed must have a value between 0 and 30081.");
		exit(1);
	}
	
	i = (ij/177)%177 + 2;
	j = ij%177 + 2;
	k = (kl/169)%178 + 1;
	l = kl%169;
	
	for (ii=1; ii<=97; ii++) {
		s = 0.0;
		t = 0.5;
		for (jj=1; jj<=24; jj++) {
			m = (((i*j)%179)*k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53*l + 1) % 169;
			if ((l*m)%64 >= 32) s += t;
			t *= 0.5;
		}
		u[ii] = s;
	}
	
	c = 362436.0 / 16777216.0;
	cd = 7654321.0 / 16777216.0;
	cm = 16777213.0 / 16777216.0;
	
	i97 = 97;
	j97 = 33;
	
	test = TRUE;
}

float ranmar(void)
/*
 C This is the random number generator proposed by George Marsaglia in 
 C Florida State University Report: FSU-SCRI-87-50
 */
{
	// int ivec;
	float uni;
	
	if (test==FALSE) {
		puts("Call the init routine rmarin() before calling ranmar().");
		exit(2);
	}
	
	uni = u[i97] - u[j97];
	if (uni < 0.0) uni += 1.0;
	u[i97] = uni;
	i97--;
	if (i97==0) i97 = 97;
	j97--;
	if (j97==0) j97 = 97;
	c -= cd;
	if (c<0.0) c += cm;
	uni -= c;
	if (uni<0.0) uni += 1.0;
	
	return uni;
}