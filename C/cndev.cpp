/*
 *  cndev.c
 *
 *  Created by Wilson on 29/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "cndev.h"
#include <math.h>

#define Pi 3.141592653589793238462643

static double A[4]={
	2.50662823884,
	-18.61500062529,
	41.39119773534,
	-25.44106049637
};
static double B[4]={
	-8.47351093090,
	23.08336743743,
	-21.06224101826,
	3.13082909833
};
static double C[9]={
	0.3374754822726147,
	0.9761690190917186,
	0.1607979714918209,
	0.0276438810333863,
	0.0038405729373609,
	0.0003951896511919,
	0.0000321767881768,
	0.0000002888167364,
	0.0000003960315187
};

double norminv(double U)
{
	/* Returns the inverse of the cumulative normal distribution
	 
	 Written by B. Moro, November 1994 */
	
	double X, R;
	
	X = U - 0.5;
	if (fabs(X)<0.42) {
		R = X*X;
		R = X*(((A[3]*R+A[2])*R+A[1])*R+A[0]) /
		((((B[3]*R+B[2])*R+B[1])*R+B[0])*R+1.0);
		return (R);
	}
	
	R = U;
	if (X > 0.)
		R=1.0-U;
	
	R = log(-log(R));
	R = C[0]+R*(C[1]+R*(C[2]+R*(C[3]+R*(C[4]+R*(C[5]+R*(C[6]+R*(C[7]+R*C[8])))))));
	if (X < 0.0)
		R = -R;
	return (R);
}

double normcdf( double X )
{
	
	double L, K, w ;
	
	double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
	double const a4 = -1.821255978, a5 = 1.330274429;
	
	L = fabs(X);
	K = 1.0 / (1.0 + 0.2316419 * L);
	w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));
	
	if (X < 0 ){
		w= 1.0 - w;
	}
	return w;
}
