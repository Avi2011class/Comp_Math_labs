import numpy as np
import math

# Параметры задачи
N = 20
h = np.pi / N

# Зададим коэффициенты системы
alpha_n = lambda n: -1 if n != 0 and n != N else 0
beta_n = lambda n: -2 - h ** 2 if n != 0 and n != N else 1
gamma_n = lambda n: -1 if n != 0 and n != N else 0
delta_n = lambda n: 2 * (h ** 2) * np.sin( n * h ) if n != 0 and n != N else 0
    
# Прямой ход
P = [0 for i in range(N + 1)]
Q = [0 for i in range(N + 1)]

P[0] = gamma_n(0) / beta_n(0)
Q[0] = -delta_n(0) / beta_n(0)

for i in range(1, N + 1):
    P[i] = gamma_n(i) / (beta_n(i) - alpha_n(i) * P[i - 1])
    Q[i] = (alpha_n(i) * Q[i - 1] - delta_n(i)) / (beta_n(i) - alpha_n(i) * P[i - 1])

# Обратный ход
x = [0 for i in range(N + 1)]
x[N] = (alpha_n(N) * Q[N - 1] - delta_n(N)) / (beta_n(N) - alpha_n(N) * P[N - 1])
for i in range(N - 1, -1, -1):
    x[i] = P[i] * x[i + 1] + Q[i]
print('u =', x, end='\n\n')
error = [abs(x[i] - np.sin(i * h)) for i in range(0, N + 1)]
print('Error = |x_i - sin(i * h)| =', error, end='\n\n')
print('Max(error) =', max(error), end='\n\n')

