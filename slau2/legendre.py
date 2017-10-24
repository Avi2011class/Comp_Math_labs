#!/usr/bin/python3

import math
import numpy as np
import random
import matplotlib.pyplot as plt
import scipy.integrate

from pylab import rcParams
rcParams['figure.figsize'] = 6,4

# Приближаемая функция
def approx(x):
    return np.sin(x * 3) 

# Размер определяемого базиса
N = 15

# Построение базиса из полиномов Лежандра
legendre_polynomials = []

legendre_polynomials.append(np.polynomial.Polynomial([1,]))
legendre_polynomials.append(np.polynomial.Polynomial([0,1]))

for k in range(2, N + 1):
    legendre_polynomials.append( 
        (2 * k - 1) / (k) * np.polynomial.Polynomial([0, 1]) * legendre_polynomials[k - 1] -
        (k - 1) / (k) * legendre_polynomials[k - 2]
       )

# Нормировка полиномов к 1 по L_2[-1, 1]
for k in range(N + 1):
    legendre_polynomials[k] /= math.sqrt(2 / (2 * k + 1))

# Функция, интегрируемая по -1, 1 для определения args[0] - 
# компоненты по построенному базису
def f(x, *args):
    global legendre_polynomials
    return legendre_polynomials[args[0]].__call__(x) * approx(x)

# Скалярные произведения апроксимируемой функции на базисные векторы
legendre_coefficients = [ scipy.integrate.quad(f, -1, 1, args=i)[0] for i in range(N + 1) ]


# Построение графика для некоторых степеней
fig = plt.figure()
ax = fig.add_subplot(111)
l = np.linspace( -1, 1, 100 )
ax.plot( l, approx(l), label='sin(3x)', linewidth=8, color='green', alpha=0.5 )

color = {1:'red', 2:'blue', 4:'purple', 6:'black'}

# В частности для 0, 1, 3 и 5
for k in [1,2,4,6]:
    result_polynomial =  np.polynomial.Polynomial([0])
    for i in range(k):
        result_polynomial += legendre_polynomials[i] * legendre_coefficients[i]

    
    ax.plot( l, result_polynomial(l), label=str(k - 1)+' степень', linewidth=1, c=color[k] )
    

ax.grid()
ax.legend()

fig.savefig('plot.png')






