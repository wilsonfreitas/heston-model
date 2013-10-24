
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

gsl_complex bsm_phi(double k, double v, double tau) {
	double m = -0.5*v*k*tau;
	gsl_complex z = gsl_complex_add_real(gsl_complex_rect(0, -1), k);
	z = gsl_complex_mul_real(z, m);
	z = gsl_complex_exp( z );
	printf("Real %f\n", GSL_REAL(z));
	printf("Imag %f\n", GSL_IMAG(z));
	return z;
}

int main() {
	gsl_complex z = bsm_phi(1, 1, 1);
	// printf("%f\n", GSL_REAL(bsm_phi(1, 1, 1)));
	return 0;
}