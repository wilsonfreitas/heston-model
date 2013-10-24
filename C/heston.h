/*
 *  heston.h
 *  bs-lewis
 *
 *  Created by Wilson on 02/10/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

//#include "gcomplex.h"

//static gcomplex heston_phi(double k, double v, double tau, void* p);
//double heston_phi_integrand(double k, void* p);
//double heston_phi_transform(double v, double tau, double X, void* p);

double heston_future_call(double F, double K, double tau, double sig, double lambda, double rho, double eta, double vbar);
double heston_call(double S, double K, double T, double r, double sig, double q, double lambda, double rho, double eta, double vbar);
double bsm_call(double S, double K, double T, double r, double sig, double q);
double black_call(double F, double K, double T, double sig);
