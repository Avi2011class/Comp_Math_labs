#!/usr/bin/python3

import math
import numpy as np
import random
import matplotlib.pyplot as plt
import scipy.integrate

from pylab import rcParams
rcParams['figure.figsize'] = 7,5


N = 10

legendre_polynomials = []

legendre_polynomials.append(np.polynomial.Polynomial([1,]))
legendre_polynomials.append(np.polynomial.Polynomial([0,1]))

for k in range(2, N + 1):
    legendre_polynomials.append( 
        (2 * k - 1) / (k) * np.polynomial.Polynomial([0, 1]) * legendre_polynomials[k - 1] -
        (k - 1) / (k) * legendre_polynomials[k - 2]
       )

for k in range(N + 1):
    legendre_polynomials[k] /= math.sqrt(2 / (2 * k + 1))


af = lambda x: legendre_polynomials[9](x) * legendre_polynomials[9](x)


print(scipy.integrate.quad(af, -1, 1))

def f(x, *args):
    global legendre_polynomials
    return legendre_polynomials[args[0]].__call__(x) * np.sin(x)




fig = plt.figure()
ax = fig.add_subplot(111)

for k in [1, 2, 6]:
    legendre_coefficients = [ scipy.integrate.quad(f, -1, 1, args=i)[0] for i in range(k) ]
    result_polynomial = np.polynomial.Polynomial(legendre_coefficients)


    l = np.linspace( -1, 1, 1000 )
    ax.plot( l, result_polynomial(l), label=str(k) )
    
ax.plot( l, np.sin(l), label='sin' )
ax.grid()
ax.legend()

fig.savefig('plot.png')






