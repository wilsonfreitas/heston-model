#!/usr/bin/env python
# encoding: utf-8

import sys
import os
from scipy import loadtxt
import numpy as np
from numpy.linalg import norm
import heston
import _heston
import bsm
# from pyevolve import G1DList
# from pyevolve import GSimpleGA
# from pyevolve import Mutators
# from pyevolve import Initializators
# from pyevolve import GAllele
# from pyevolve import Consts
# from pyevolve import DBAdapters

fut_datafile = '20100526_USDBRL_FUT.txt'
vol_datafile = '20100526_USDBRL_IMPVOL.txt'

vol_data = loadtxt(vol_datafile)
terms = vol_data[0,1:].T
strikes = vol_data[1:,0]
vols = vol_data[1:,1:]

fut_data = loadtxt(fut_datafile)
futs = fut_data[:,1]
vanilla = np.zeros((len(strikes), 1))
term = 1
for i, k in enumerate(strikes):
	# for j, t in enumerate(terms):
	vanilla[i,0] = bsm.bsmprice(futs[term], k, terms[term], 0, vols[i,term], 0)[0]

def heston_evaluate(chromo):
	"""docstring for heston_evaluate"""
	v = chromo[0]
	vbar = chromo[1]
	lambd = chromo[2]
	eta = chromo[3]
	rho = chromo[4]
	diffs = np.zeros(len(strikes)*1)
	n = 0
	for i, k in enumerate(strikes):
		# for j, t in enumerate(terms):
		diffs[n] = vanilla[i,0] - _heston.ucall(futs[term], k, terms[term], v, vbar, lambd, eta, rho)
		n += 1
	score = norm(diffs)
	# print(score)
	return diffs

from scipy.optimize import leastsq
p0 = (0.03940606,  0.03301662,  2.48004912,  0.56642651,  0.63548062)
plsq = leastsq(heston_evaluate, p0)
print(plsq[0])
parms = plsq[0]

import matplotlib.pyplot as plt

hes = np.zeros(len(strikes))
for i, k in enumerate(strikes):
	hes[i] = _heston.ucall(futs[term], k, terms[term], parms[0], parms[1], parms[2], parms[3], parms[4])

plt.plot(strikes, vanilla)

# plt.title('Least-squares fit to noisy data')
# plt.legend(['Fit', 'Noisy', 'True'])
# plt.show()

# # Genome instance
# setOfAlleles = GAllele.GAlleles()
# setOfAlleles.add(GAllele.GAlleleRange(0.0, 5.0, True))
# setOfAlleles.add(GAllele.GAlleleRange(0.0, 5.0, True))
# setOfAlleles.add(GAllele.GAlleleRange(5.0, 15.0, True))
# setOfAlleles.add(GAllele.GAlleleRange(0.0, 5.0, True))
# setOfAlleles.add(GAllele.GAlleleRange(-1, 1, True))
# 
# genome = G1DList.G1DList(5)
# genome.setParams(allele=setOfAlleles)
# genome.evaluator.set(heston_evaluate)
# genome.mutator.set(Mutators.G1DListMutatorAllele)
# genome.initializator.set(Initializators.G1DListInitializatorAllele)
# 
# # Genetic Algorithm Instance
# ga = GSimpleGA.GSimpleGA(genome)
# # ga.setGenerations(100)
# # ga.setPopulationSize(100)
# # ga.setMutationRate(0.01)
# # ga.setCrossoverRate(0.90)
# ga.minimax = Consts.minimaxType["minimize"]
# # ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
# 
# id = 'default-11'
# 
# sqlite_adapter = DBAdapters.DBSQLite(identify="heston-"+id, resetIdentify=True,
# 	resetDB=False)
# ga.setDBAdapter(sqlite_adapter)
# 
# # Do the evolution, with stats dump
# # frequency of 10 generations
# ga.evolve(freq_stats=5)
# 
# # Best individual
# print ga.bestIndividual()
# f = file('output-' + id + '.txt', 'w')
# f.write(str(ga.bestIndividual()))
# f.close()

