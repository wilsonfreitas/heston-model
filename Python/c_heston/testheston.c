
#include "heston.h"
#include <stdio.h>

int main() {
	double h = heston_u_call(100, 100, 0.5, 0.2, 0.1, 2, 0.2, -0.5);
	printf("heston = %.20f\n", h);
	return 0;
}