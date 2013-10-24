/*
 *  heston.cpp
 *  bs-lewis
 *
 *  Created by Wilson on 02/10/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "cndev.h"
#include "gcomplex.h"

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

typedef struct {
	double lambda;
	double rho;
	double eta;
	double vbar;
} heston_params;

typedef struct {
	heston_params* heston_p;
	double v;
	double x;
	double tau;
} heston_integrand_params;

using namespace gslcpp;

gcomplex heston_phi(double k_, double v, double tau, void* p) {
	
	gcomplex I = gcomplex::rect(0,1);
	gcomplex k = I*0.5 + k_;
	
	heston_params* params = (heston_params*)p;
	double lambda = params->lambda;
	double rho = params->rho;
	double eta = params->eta;
	double vbar = params->vbar;
	
	gcomplex b = I*(k*rho*eta) + lambda;
	gcomplex c = (I-k)*k*eta*eta*(-1);
	gcomplex d = gcomplex::sqrt( b*b + c );
	gcomplex g = (b - d)/(b + d);
	gcomplex T_m = (b - d)/(eta*eta);
	gcomplex T = T_m * ( gcomplex::exp(d*(-tau))*(-1) + 1 )/( g*gcomplex::exp(d*(-tau))*(-1) + 1 );
	gcomplex W = ( T_m*tau + gcomplex::log( ( g*gcomplex::exp(d*(-tau))*(-1) + 1 )/( g*(-1) + 1 ) )*(-2)/(eta*eta) )*lambda*vbar;
	
	return gcomplex::exp(W + T*v);
}

double heston_phi_integrand(double k, void* p) {
	heston_integrand_params* params = (heston_integrand_params*)p;
	
	double v = params->v;
	double tau = params->tau;
	double X = params->x;
	
	gcomplex I = gcomplex::rect(0,-1);
	gcomplex integrand = gcomplex::exp(I*k*X) * heston_phi(k, v, tau, params->heston_p);
	return 2 * integrand.real() / (k*k + 1.0/4.0);
}

double heston_phi_transform(double v, double tau, double X, void* p) {
	heston_params* hparams = (heston_params*)p;
	
	heston_integrand_params params;
	params.heston_p = hparams;
	params.v = v;
	params.x = X;
	params.tau = tau;
	
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &heston_phi_integrand;
	F.params = &params;
	double result, error;
	gsl_integration_qags (&F, 0, 100, 0, 1e-7, 1000, w, &result, &error);
	return result;
}

double heston_future_call(double F, double K, double tau, double sig, double lambda, double rho, double eta, double vbar) {
	double v = sig*sig;
	double X = log(F/K);
	
	heston_params hparams;
	hparams.lambda = lambda;
	hparams.eta = eta;
	hparams.rho = rho;
	hparams.vbar = vbar;
	
	double integral = heston_phi_transform(v, tau, X, &hparams);
	return F - (sqrt(K*F)/(2*M_PI)) * integral;
}

double heston_call(double S, double K, double tau, double r, double sig, double q, double lambda, double rho, double eta, double vbar) {
	double F = S*exp((r-q)*tau);
	return exp(-r*tau)*heston_future_call(F, K, tau, sig, lambda, rho, eta, vbar);
}

double bsm_call(double S, double K, double T, double r, double sig, double q) {
	double d1 = ( log(S/K) + (r - q + sig*sig/2.0)*T ) / ( sig*sqrt(T) );
	double d2 = d1 - sig*sqrt(T);
	return S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
}

double black_call(double F, double K, double T, double sig) {
	double d1 = ( log(F/K) + sig*sig*T/2.0 ) / ( sig*sqrt(T) );
	double d2 = d1 - sig*sqrt(T);
	return F*normcdf(d1) - K*normcdf(d2);	
}

