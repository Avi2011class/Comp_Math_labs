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

# LU разложение матрицы А
def LU_decomposition(A):
    global N
    L = np.matrix(np.zeros(shape=(N, N)))
    U = np.matrix(np.zeros(shape=(N, N)))
    for j in range(N):
        U[0, j] = A[0, j]
        L[j, 0] = A[j, 0] / U[0, 0]

    for i in range(1, N):
        for j in range(i, N):
            su = 0
            for k in range(0, i): su += L[i, k] * U[k, j]
            U[i, j] = A[i, j] - su

            su = 0
            for k in range(0, i): su += L[j, k] * U[k, i]
            L[j, i] = 1 / U[i, i] * (A[j, i] - su)

    return L, U

# Решение системы Lx = F
def L_solve(L, F):
    global N
    F_tmp = F.copy()
    solution = np.zeros(shape=N)

    for i in range(0, N):
        for j in range(0, i):
            F_tmp[i, 0] -= solution[j] * L[i, j]
        solution[i] = F_tmp[i, 0] / L[i, i]

    return np.matrix(solution).T

# Решение системы Ux =F
def U_solve(U, F):
    global N
    F_tmp = F.copy()
    solution = np.zeros(shape=N)

    for i in range(N - 1, -1, -1):
        for j in range(N - 1, i, -1):
             F_tmp[i, 0] -= solution[j] * U[i, j]
        solution[i] = F_tmp[i, 0] / U[i, i]

    return np.matrix(solution).T

#  Объединение предыдущих методов
def LU_system_solve(A, F):
    L, U = LU_decomposition(A)
    solution = L_solve(L, F)
    solution = U_solve(U, solution)
    return solution

# ------------------------------------------------
A = generateA()
F = generateF()

print('Матрица системы:')
print(A)
print('Столбец из правой части:')
print(F)
print('-' * 80)

L, U = LU_decomposition(A)

print('Матрица L:')
print(L)
print('Матрица U:')
print(U)
print('-' * 80)
print('Отклонение LU от A по норме Фробениуса: ', np.linalg.norm(L * U - A, ord='fro'))
print('-' * 80)


solution = L_solve(L, F)
print('Решение системы Lx = F:')
print(solution)

solution = U_solve(U, solution)
print('Решение системы Uy = x (совпадает с решением изначальной системы):')
print(solution)

err = solution - A.I * F
print('y - A^{-1}*F, вычисленное для проверки:')
print(err)
print('-' * 80)


# ----------------
print('Будем давать приращения вектору F и оценивать число обусловленности матрицы')
mu_arr = []
for i in range(100):
    el = np.linalg.norm(F, ord=np.inf)
    df = np.random.uniform(-el * 0.01, el * 0.01, size=(1, N))
    F_new = F + df
    solution_new = LU_system_solve(A, F_new)
    du = solution_new - solution

    mu = np.linalg.norm(du, ord='fro') / np.linalg.norm(solution_new, ord='fro') / (np.linalg.norm(df, ord='fro') / np.linalg.norm(F_new, ord='fro'))
    print('mu =', mu)
    mu_arr.append(mu)

print('Максимальное из найденных mu:', np.max(mu_arr))




