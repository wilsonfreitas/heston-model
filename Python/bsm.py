#!/usr/bin/env python
# encoding: utf-8

from numpy import sqrt
from scipy import exp, log
from scipy.stats import norm

def d1d2(S, K, T, r, v, q):
	d1 = ( log(S/K) + (r - q + 0.5*v**2)*T ) / (v*sqrt(T))
	d2 = d1 - v*sqrt(T)
	return d1, d2

def put(S, K, T, r, v, q):
	d1, d2 = d1d2(S, K, T, r, v, q)
	N1 = norm.cdf(d1)
	N2 = norm.cdf(d2)
	
	return -S*exp(-q*T)*(1-N1) + K*exp(-r*T)*(1-N2)


def call(S, K, T, r, v, q):
	d1, d2 = d1d2(S, K, T, r, v, q)
	N1 = norm.cdf(d1)
	N2 = norm.cdf(d2)
	
	return S*exp(-q*T)*N1 - K*exp(-r*T)*N2


if __name__ == '__main__':
	print 'BSM call: ', call(100.00, 110.0, 1.0, 0.0, 0.2, 0.0)
	print 'BSM put : ', put (100.00, 110.0, 1.0, 0.0, 0.2, 0.0)
