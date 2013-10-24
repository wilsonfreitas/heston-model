#!/usr/bin/env python
# encoding: utf-8
"""
bsm.py

Created by Wilson on 2010-08-30.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

from numpy import meshgrid, sqrt, diff
from scipy import inf, pi, exp, linspace, zeros, real, imag, array, log
from scipy.stats import norm
from scipy.integrate import quad
from pylab import figure, show, plot, axvline
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def bsmprice(S, K, T, r, v, q):
	d1 = ( log(S/K) + (r - q + 0.5*v**2)*T ) / (v*sqrt(T))
	d2 = d1 - v*sqrt(T)
	N1 = norm.cdf(d1)
	N2 = norm.cdf(d2)
	S0 = S*exp(-q*T)
	K0 = K*exp(-r*T)
	
	call = S0*N1 - K0*N2
	put = -S0*(1-N1) + K0*(1-N2)
	
	return call, put


def bsm_phi(k, v, tau):
	return exp( -0.5*v*k*(k - 1j)*tau )


def bsm_phi_integrand(k, v, tau, b):
	phi = bsm_phi(k + b*1j, v, tau)
	return ( 1./( (k + 1j*b)*(k + 1j*(b - 1)) ) ) * phi


def bsm_phi_transform(v, tau, b, x):
	integrand = lambda k: real( exp(-1j*k*x) * bsm_phi_integrand(k, v, tau, b) )
	return (1.0/pi) * exp(b*x) * quad(integrand, 0, 100)[0]


def bsm_call(S, K, tau, r, v, q, b):
	v = v**2.0
	F = S*exp((r-q)*tau)
	x = log(F/K)
	integral = bsm_phi_transform(v, tau, b, x)
	return S * exp(-q*tau) - K * exp(-r*tau) * integral


def bsm_call2(S, K, tau, r, v, q, b):
	v = v**2.0
	F = S*exp((r-q)*tau)
	x = log(F/K)
	integral = bsm_phi_transform(v, tau, b, x)
	return - K * exp(-r*tau) * integral


if __name__ == '__main__':
	print 'BSM       : ', bsmprice (100.00, 110.0, 1.0, 0.0, 0.2, 0.0)[0]
	print 'Formula I : ', bsm_call (100.00, 110.0, 1.0, 0.0, 0.2, 0.0, 0.5)
	print 'Formula II: ', bsm_call2(100.00, 110.0, 1.0, 0.0, 0.2, 0.0, 2.0)

