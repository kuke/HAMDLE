#! /bin/python

import sys
import matplotlib.pyplot as  plt
import numpy as np

X = np.arange(-1, 1.05, 0.05)
order = 5

def generate_Legendre_matrix(order, X):
	m = X.size
	Legendre_mat = [[0 for col in range(m)] for row in range(order+1)]
	Legendre_mat[0] = [1 for i in range(0, m)]
	Legendre_mat[1] = X
	for i in range(2, order+1):
		Legendre_mat[i] = (2*i-1)*np.multiply(X,Legendre_mat[i-1])-(i-1)*Legendre_mat[i-2]
		Legendre_mat[i] = Legendre_mat[i]/i
	return Legendre_mat

Legendre_mat = generate_Legendre_matrix(order, X)
plt.figure("Legendre polynomials")

if 5==order:
	plt.plot(X, Legendre_mat[0], '-',  label="n=0",)
	plt.plot(X, Legendre_mat[1], '--', label="n=1")
	plt.plot(X, Legendre_mat[2], '-.', label="n=2")
	plt.plot(X, Legendre_mat[3], ':',  label="n=3")
	plt.plot(X, Legendre_mat[4], '-+', label="n=4")
	plt.plot(X, Legendre_mat[5], '-^', label="n=5")
else: 
	for i in range(0, order):
		plt.plot(X, Legendre_mat[i], label="n="+str(i))

plt.legend(loc='best', prop={'size':16})
plt.axis([-1, 1, -1, 2])	
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$P_n(x)$', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()
