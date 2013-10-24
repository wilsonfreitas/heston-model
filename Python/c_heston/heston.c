
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

typedef struct {
	double lambda;
	double rho;
	double eta;
	double v;
	double vbar;
} heston_params;

typedef struct {
	heston_params* heston_p;
	double x;
	double tau;
} heston_integrand_params;

gsl_complex heston_phi(double k_, double tau, void* p) {
	
	// k = k_ + 0.5*I
	gsl_complex k = gsl_complex_add_real(gsl_complex_rect(0, 0.5), k_);
	
	heston_params* params = (heston_params*)p;
	double lambda = params->lambda;
	double rho = params->rho;
	double eta = params->eta;
	double vbar = params->vbar;
	double v = params->v;
	
	// b = lambda + I*rho*eta*k
	gsl_complex b = gsl_complex_add_real(gsl_complex_mul(k, gsl_complex_rect(0, rho*eta)), lambda);
	// c = eta**2 * k * (k - I)
	gsl_complex c = gsl_complex_mul(gsl_complex_sub(gsl_complex_rect(0, 1), k), gsl_complex_mul_real(k, -eta*eta));
	// d = sqrt(b**2 + c)
	gsl_complex d = gsl_complex_sqrt(gsl_complex_add(gsl_complex_pow_real(b, 2.0), c) );
	// g = (b - d)/(b + d)
	gsl_complex g = gsl_complex_div(gsl_complex_sub(b, d), gsl_complex_add(b, d));
	// T_m = (b - d)/(eta*eta)
	gsl_complex T_m = gsl_complex_div_real(gsl_complex_sub(b, d), eta*eta);
	// T_a = -exp(-d*tau)
	gsl_complex T_a = gsl_complex_mul_real(gsl_complex_exp(gsl_complex_mul_real(d, -tau)), -1);
	// T_1 = 1 + T_a
	gsl_complex T_1 = gsl_complex_add_real(T_a, 1);
	// T_2 = 1 + g*T_a
	gsl_complex T_2 = gsl_complex_add_real(gsl_complex_mul(T_a, g), 1);
	// T = T_m * T_1 / T_2
	// T = T_m * ( exp(d*(-tau))*(-1) + 1 )/( g*exp(d*(-tau))*(-1) + 1 );
	gsl_complex T = gsl_complex_mul(T_m, gsl_complex_div(T_1, T_2));
	// T_t = tau * T_m
	gsl_complex T_t = gsl_complex_mul_real(T_m, tau);
	// g_1 = 1 - g
	gsl_complex g_1 = gsl_complex_add_real(gsl_complex_mul_real(g, -1), 1);
	// T_l = ln(T_2/g_1)
	gsl_complex T_l = gsl_complex_log(gsl_complex_div(T_2, g_1));
	// T_e = -2/eta**2 * T_l
	gsl_complex T_e = gsl_complex_mul_real(T_l, -2.0/(eta*eta));
	// W = lambda*vbar*(T_t + T_e)
	gsl_complex W = gsl_complex_mul_real(gsl_complex_add(T_t, T_e), lambda*vbar);
	
	// exp(W + T*v)
	return gsl_complex_exp(gsl_complex_add(W, gsl_complex_mul_real(T, v)));
}

double heston_phi_integrand(double k, void* p) {
	heston_integrand_params* params = (heston_integrand_params*)p;
	
	double tau = params->tau;
	double X = params->x;
	
	gsl_complex phi = heston_phi(k, tau, params->heston_p);
	gsl_complex e = gsl_complex_mul_real(gsl_complex_rect(0,-1), k*X);
	gsl_complex integrand = gsl_complex_mul(gsl_complex_exp(e), phi);
	return 2 * GSL_REAL(integrand) / (k*k + 1.0/4.0);
}

double heston_phi_transform(double tau, double X, void* p) {
	heston_params* hparams = (heston_params*)p;
	
	heston_integrand_params params;
	params.heston_p = hparams;
	params.x = X;
	params.tau = tau;
	
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &heston_phi_integrand;
	F.params = &params;
	double result, error;
	gsl_integration_qags (&F, 0, 100, 0, 1e-7, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}

double heston_ucall(double F, double K, double tau, double v, double vbar, double lambda, double eta, double rho) {
	double X = log(F/K);
	
	heston_params hparams;
	hparams.lambda = lambda;
	hparams.eta = eta;
	hparams.rho = rho;
	hparams.vbar = vbar;
	hparams.v = v;
	
	double integral = heston_phi_transform(tau, X, &hparams);
	return F - (sqrt(K*F)/(2*M_PI)) * integral;
}

double heston_call(double S, double K, double tau, double r, double q, 
double v, double vbar, double lambda, double eta, double rho) {
	double F = S*exp((r-q)*tau);
	return exp(-r*tau)*heston_ucall(F, K, tau, v, vbar, lambda, eta, rho);
}

