import numpy as np
import math
import os
import sys

import decimal
decimal.getcontext().prec = 80

N = 5
lmbd = 0.01
h = 2 * np.pi / N

K_generator =   lambda m, n:    h ** 2 * np.abs(n - m)
f =             lambda x:       (1 + 2 * lmbd) * np.cos(x / 2) ** 2 - lmbd * (x ** 2 + np.pi ** 2) / 2
F_generator =   lambda n:      f(-np.pi + (n - 0.5) * h)     

A = np.identity(N)
for i in range(N):
    for j in range(N):
        A[i][j] -= lmbd * K_generator(i + 1, j + 1)
        
F = np.matrix([F_generator(n) for n in range(1, N + 1)]).T

A = np.matrix([ [ decimal.Decimal(A[i, j]) for j in range(N) ] for i in range(N) ])
F = np.matrix([ [ decimal.Decimal(F[i, j]) for j in range(1) ] for i in range(N) ])


print('Матрица системы:')
print(A)
print('Столбец из правой части:')
print(F)
print('-' * 80)


L = np.matrix([ [ decimal.Decimal(0) for j in range(N) ] for i in range(N) ])
U = np.matrix([ [ decimal.Decimal(0) for j in range(N) ] for i in range(N) ])

# LU decompositio
for j in range(N):
    U[0, j] = A[0, j]
    L[j, 0] = A[j, 0] / U[0, 0]

for i in range(1, N):
    for j in range(i, N):
        
        su = decimal.Decimal(0)
        for k in range(0, i):
            su += L[i, k] * U[k, j]
        U[i, j] = A[i, j] - su

        su = decimal.Decimal(0)
        for k in range(0, i):
            su += L[j, k] * U[k, i]
        L[j, i] = decimal.Decimal(1) / U[i, i] * (A[j, i] - su)
        
print('Матрица L:')
print(L)
print('Матрица U:')
print(U)
print('-' * 80)
print('Отклонение LU от A по норме Фробениуса: ', np.linalg.norm(L * U - A, ord='fro'))
print('-' * 80)


solution = np.matrix([ [ decimal.Decimal(0) for j in range(1) ] for i in range(N) ])

F_tmp = F[:]

for i in range(0, N):
    for j in range(0, i):
        F_tmp[i, 0] -= solution[j, 0] * L[i, j]
    solution[i, 0] = F_tmp[i, 0] / L[i, i]
    
print(L * solution - F)



