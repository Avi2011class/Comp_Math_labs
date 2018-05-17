import numpy as np
import math
import functools
import time
import tqdm

import matplotlib
from matplotlib.animation import FuncAnimation

import matplotlib.pyplot as plt
matplotlib.use("Qt5Agg")

from scipy.linalg import solve_banded

N = 50     # Число точек на стороне квадрата

# Начальное положение частицы
x_0 = 1/4
y_0 = 1/2
sigma_0 = 1/10

# Волновое число
k = 200

# Характеристики "пика" потенциала
x_v = 1/2
y_v = 1/2
sigma_v = 1/40

# -----------------------------------------------------------

# Технические константы, необходимые для вычислений
h = 1 / (N - 1)
tau = h / k / 2
rho = 1j * tau / 4 / h ** 2

frame_size = 10 # Количество итераций метода между двумя кадрами анимации
n_frames = 100 # Число кадров
delay = 50 # Задержка кадров

# ---------------------------------------------------------

# Генерация поля с начальным распределением температуры
# Значения на границе обязательно должны быть равны нулям

def generate_field():
    arr = np.zeros(shape=(N, N), dtype=np.complex128)

    for i in range(N):
        for j in range(N):
            x = i / (N - 1)
            y = j / (N - 1)

            arr[i, j] = np.exp( -1.0j * k * x ) * \
                np.exp( -( (x - x_0) ** 2 + (y - y_0) ** 2) / sigma_0 ** 2 )

    return arr


# Вспомогательная функция V, вычисление потенциала
def V(x, y):
    return k * k * np.exp( -( (x - x_v) ** 2 + (y - y_v) ** 2 ) / sigma_v ** 2 )

v = np.zeros(shape=(N, N), dtype=np.complex128)
for i in range(N):
    for j in range(N):
        x = i / (N - 1)
        y = j / (N - 1)
        v[i, j] = 1j * tau * V(x, y) / 4


banded_matrix = np.zeros(shape=(3, N), dtype=np.complex128)
b_matrix = np.zeros(N, np.complex128)

# Сама схема расщепления, неявная, один шаг
def computation_step(arr, buff):
    global banded_matrix, b_matrix

    for i in range(1, N - 1):
        for j in range(N):
            banded_matrix[0, j] = -rho
            banded_matrix[1, j] = 1 + 2 * rho + v[i, j]
            banded_matrix[2, j] = -rho

        banded_matrix[0, 0] = 0
        banded_matrix[0, 1] = 0

        banded_matrix[-1, -1] = 0
        banded_matrix[-1, -2] = 0

        banded_matrix[1, 0] = 1
        banded_matrix[1, -1] = 1

        b_matrix[0] = 0
        b_matrix[-1] = 0

        for j in range(1, N - 1):
            b_matrix[j] = rho * arr[i - 1, j] + \
                (1 - 2 * rho - v[i, j]) * arr[i, j] + \
                rho * arr[i + 1, j]

        solution = solve_banded((1, 1), banded_matrix, b_matrix)
        buff[i, :] = solution

    for j in range(1, N - 1):
        banded_matrix = np.zeros(shape=(3, N), dtype=np.complex128)
        b_matrix = np.zeros(N, np.complex128)

        for i in range(N):
            banded_matrix[0, i] = -rho
            banded_matrix[1, i] = 1 + 2 * rho + v[i, j]
            banded_matrix[2, i] = -rho

        banded_matrix[0, 0] = 0
        banded_matrix[0, 1] = 0

        banded_matrix[-1, -1] = 0
        banded_matrix[-1, -2] = 0

        banded_matrix[1, 0] = 1
        banded_matrix[1, -1] = 1

        b_matrix[0] = 0
        b_matrix[-1] = 0

        for i in range(1, N - 1):
            b_matrix[i] = rho * buff[i, j - 1] + \
                (1 - 2 * rho - v[i, j]) * buff[i, j] + \
                rho * arr[i, j + 1]

        solution = solve_banded((1, 1), banded_matrix, b_matrix)
        arr[:, j] = solution

#
# ----------------------------------------------------
#

arr = generate_field()

buff = np.zeros(shape=(N, N), dtype=np.complex128)


#
# -----------------------------------------------------
#
#fig = plt.figure(figsize=(6, 10))
#ax = fig.subplots(2)

#l = np.linspace(0, 1, N)
#xx, yy = np.meshgrid(l, l)

#ax[0].grid()
#ax[0].contourf(xx, yy, np.absolute(arr), cmap='hot')
#ax[0].set_title('Модуль волновой функции')

#ax[1].grid()
#ax[1].contourf(xx, yy, np.angle(arr), cmap='hot')
#ax[1].set_title('Аргумент волновой функции')

#fig.savefig('initial.png')
#
# ------------------------------------------------------
#

fig = plt.figure(figsize=(6, 6))

ax = fig.add_subplot(111)

l = np.linspace(0, 1, N)
xx, yy = np.meshgrid(l, l)

ax.grid()
a0 = ax.contourf(xx, yy, np.absolute(arr), cmap='hot')
ax.set_title('Модуль волновой функции')


def update(frame):
    ax.clear()

    for i in range(frame_size):
        computation_step(arr, buff)

    a0 = ax.contourf(xx, yy, np.absolute(arr), 100,  cmap='jet')

    return a0,

anim = FuncAnimation(fig, update, tqdm.tqdm(list(range(1, n_frames))), interval=delay)
anim.save('animation.mp4', writer='ffmpeg')
#plt.show()
#for i in tqdm.tqdm(range(1000)):
#    computation_step(arr, buff)


#
# ------------------------------------------------------
#
#







