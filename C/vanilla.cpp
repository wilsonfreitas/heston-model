/*
 *  vanilla.c
 *  bs-lewis
 *
 *  Created by Wilson on 30/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include "ranmar.h"
//#include "cndev.h"
#include "heston.h"

#define MAX0(a) (a > 0.0 ? a : 0.0)
#define randn (norminv(ranmar()))

int main (int argc, const char * argv[]) {
	
	const gsl_rng_type* rng_type;
	gsl_rng* rng;
	gsl_rng_env_setup();
	rng_type = gsl_rng_default;
	rng = gsl_rng_alloc (rng_type);
	
	printf ("generator type: %s\n", gsl_rng_name (rng));
	printf ("seed = %lu\n", gsl_rng_default_seed);
	printf ("first value = %lu\n", gsl_rng_get (rng));
	
	/* random seeds for the test case: */
//	int ij, kl, i;
//	int un = (int) time(NULL);
//	// ij = 1802;
//	ij = un % 31328;
//	kl = 9373;
//	rmarin(ij,kl);
	
	// vanilla parameters
	double S0 = 100;
	double K = 100;
	double r = 0.1;
	double q = 0.0;
	double T = 1;
	double sig = 0.2;
	// heston parameters
	// 1.80317 0.410092 -0.711811 0.0348168 0.203653
	// 2.96618 0.0890053 0.44904 0.0505632 0.141391
	// 2.47495 0.19612 -0.583731 0.0515928 0.127311
	// 3.98234 0.001 0.226036 0.0530761 0.00515548
	// 0.170288 0.249738 -0.764972 0.285047 0.124602
	// 2.94174 0.0588502 -0.870883 0.0447049 0.165728
	// 2.2822 0.00498456 0.375828 0.0590305 0.103704
	// 4.17006 0.499087 -0.875984 0.0259152 0.271258
	double lambda = 0.291935;
	double eta = 0.109045;
	double rho = -0.825172;
	double vbar = 0.245522;
	double v0 = 0.0792829;
	
	int scen = 100000;
	int steps = 252*T;
	double dt = 1.0/steps;
	
	double bs = 0.0;
	double hs = 0.0;
	double bs_sq = 0.0;
	double hs_sq = 0.0;
	double bs_std = 0.0;
	double hs_std = 0.0;
	double S, H, V;
	double z1, z2;
	
	for (int i=0; i<scen; i++) {
		z1 = gsl_ran_gaussian(rng, 1);
		
		S = S0*exp( (r - q - sig*sig/2.0)*T + sig*sqrt(T)*z1 );
		
		V = v0;
		H = S0;
		for (int j=0 ; j<steps ; j++) {
			z1 = gsl_ran_gaussian(rng, 1);
			z2 = rho*z1 + sqrt(1 - rho*rho)*gsl_ran_gaussian(rng, 1);
			
			H = H + (r - q)*H*dt + sqrt(V*dt)*H*z1;
			V = V - lambda*( MAX0(V) - vbar )*dt + eta*sqrt( MAX0(V)*dt )*z2;
//			V = fabs(V);
//			V = MAX0(V);
		}
		
		double bs_payoff = exp(-r*T)*MAX0(S - K);
		double hs_payoff = exp(-r*T)*MAX0(H - K);
		
		bs += bs_payoff;
		hs += hs_payoff;
		
		bs_sq += bs_payoff*bs_payoff;
		hs_sq += hs_payoff*hs_payoff;
	}
	
	gsl_rng_free(rng);
	
	bs /= scen;
	hs /= scen;
	bs_std = (bs_sq/scen - bs*bs)/sqrt(scen);
	hs_std = (hs_sq/scen - hs*hs)/sqrt(scen);
	
	double bs_f = bsm_call(S0, K, T, r, sig, q);
	double hs_f = heston_call(S0, K, T, r, v0, q, lambda, rho, eta, vbar);
	
	printf("\n * Simulation\n");
	printf("Black-Scholes-Merton %.10f %.10f\n", bs, bs_std);
	printf("Heston               %.10f %.10f\n", hs, hs_std);
	printf(" * Formula\n");
	printf("Black-Scholes-Merton %.10f\n", bs_f);
	printf("Heston               %.10f\n", hs_f);
	printf("\n");
	
    return 0;
}

