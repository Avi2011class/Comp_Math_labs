import numpy as np
import math
import os
import sys
import scipy.linalg

N = 5
lmbd = 0.01
h = 2 * np.pi / N

K_generator =   lambda m, n:    h ** 2 * np.abs(n - m)
f =             lambda x:       (1 + 2 * lmbd) * np.cos(x / 2) ** 2 - lmbd * (x ** 2 + np.pi ** 2) / 2
F_generator =   lambda n:       f(-np.pi + (n - 0.5) * h)     

# Создание матрицы А
def generateA():
    global N, lmbd, h, K_generator, F_generator
    
    A = np.identity(N)
    for i in range(N):
        for j in range(N):
            A[i][j] -= lmbd * K_generator(i + 1, j + 1)
    
    return np.matrix(A)

# Создание столбца F
def generateF():
    return np.matrix([F_generator(n) for n in range(1, N + 1)]).T

# Решение системы Ax=B методом простых итераций, начиная с start_x и прекращая в тот момент, 
# когда бесконечная норма соседних станосится меньше eps
def iterations(A, B, start_x, eps):
    tau = 2 / np.linalg.norm(A, ord=np.inf)
    x = start_x
    n_iterations = 0
    while True:
        n_iterations += 1
        x_new = (np.identity(A.shape[0]) - tau * A)* x + tau * B
        err = np.linalg.norm(x_new - x)
        x = x_new
        if err < eps:
            return x_new, n_iterations
    
A = generateA()
F = generateF()
print('Матрица системы:')
print(A)
print('Столбец из правой части:')
print(F)
print('-' * 80)


solution, n = iterations(A, F, np.zeros(shape=(N, 1)), 1e-6)

print('Полученное решение:')
print(solution)
print('Число необходимых итераций:', n)

print('-' * 80)
print('||A.I * F - x||_1, для проверки и оценки ощибки:')
print(np.linalg.norm(A.I * F - solution, ord=np.inf))
