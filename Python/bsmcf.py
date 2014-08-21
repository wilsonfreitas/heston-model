#!/usr/bin/env python
# encoding: utf-8

from numpy import meshgrid, sqrt, diff
from scipy import inf, pi, exp, linspace, zeros, real, imag, array, log
from scipy.stats import norm
from scipy.integrate import quad

import bsm

def cf(k, v, tau):
	'''Black-Scholes-Merton characteristic function'''
	return exp( -0.5*v*k*(k - 1j)*tau )


def cf_transform(v, tau, b, x):
	'''Black-Scholes-Merton characteristic function transform'''
	def cf_integrand(k, v, tau, b):
		phi = cf(k + b*1j, v, tau)
		return ( 1./( (k + 1j*b)*(k + 1j*(b - 1)) ) ) * phi
	
	integrand = lambda k: real( exp(-1j*k*x) * cf_integrand(k, v, tau, b) )
	return (1.0/pi) * exp(b*x) * quad(integrand, 0, 100)[0]


def call(S, K, T, r, v, q, b=0.5):
	v = v**2.0
	F = S*exp((r-q)*T)
	x = log(F/K)
	integral = cf_transform(v, T, b, x)
	return S * exp(-q*T) - K * exp(-r*T) * integral


def bsm_call2(S, K, tau, r, v, q, b=2):
	v = v**2.0
	F = S*exp((r-q)*tau)
	x = log(F/K)
	integral = cf_transform(v, tau, b, x)
	return - K * exp(-r*tau) * integral


if __name__ == '__main__':
	print 'BSM                : ', bsm.call (100.00, 110.0, 1.0, 0.0, 0.2, 0.0)
	print 'BSM Characteristic : ', call (100.00, 110.0, 1.0, 0.0, 0.2, 0.0)
	print 'Formula II: ', bsm_call2(100.00, 110.0, 1.0, 0.0, 0.2, 0.0, 2.0)

