import numpy as np
import math
import random

def LDM_decomposition(A):
    N = len(A)
    L = np.matrix(np.zeros(shape=(N, N)))
    D = np.matrix(np.zeros(shape=(N, N)))
    M = np.matrix(np.zeros(shape=(N, N)))
    for j in range(N):
        M[0, j] = A[0, j]
        L[j, 0] = A[j, 0] / M[0, 0]

    for i in range(1, N):
        for j in range(i, N):
            su = 0
            for k in range(0, i): su += L[i, k] * M[k, j]
            M[i, j] = A[i, j] - su

            su = 0
            for k in range(0, i): su += L[j, k] * M[k, i]
            L[j, i] = 1 / M[i, i] * (A[j, i] - su)

    for i in range(N):
        D[i, i] = M[i, i]
        M[i] /= M[i, i]

    return L, D, M


Mat = np.matrix(np.random.randint(-10, 10, (3, 3)).astype(float))
print(Mat)
if np.abs(np.linalg.det(Mat)) < 0.001:
    print('Матрица вырождена, LDM разложение невозможно')
    exit(0)
else:
    print('Матрица невырожена, LDM разложение ниже')



print('-' * 50)
L, D, M = LDM_decomposition(Mat)

print(L, end='\n\n')
print(D, end='\n\n')
print(M.T, end='\n\n')

print('-' * 50)
print(L * D * M, end='\n\n')
print(L * D * M - Mat, end='\n\n')
