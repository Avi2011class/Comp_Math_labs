import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 6,6

df = pd.read_csv("bessel_log.csv")
df['Error'] = np.abs(df['custom'] - df['boost'])
print(" Принимая во внимание тот факт, что эталонные значения \
функции были вычислены функцией библиотеки Boost с точностью 10^(-20), \
а для новой кастомной функции задана предельная погрешность 10^(-6), \
видим, что максимальное расхождение вычисления в точках отрезка -10, 10 с шагом 0.01 равно:")
print('{:.12f}={:g}'.format(np.max(df[df['x'] <= 10]['Error'].values), np.max(df[df['x'] <= 10]['Error'].values)))
print("А на отрезке 20 - 220 с шагом в 10:")
print('{:.12f}={:g}'.format(np.max(df[df['x'] > 10]['Error'].values), np.max(df[df['x'] > 10]['Error'].values)))
print("Дальнейшие вычисления по данной формуле являются бессмысленными в силу необходимовсти в слишком большом количестве членов ряда Тейлора")

fig, axarr = plt.subplots(2)
axarr[0].set_title('График вычисленной части функции $J_0(x)$, $J_0\'(x)$')
axarr[0].plot(df[df['x'] <= 10]['x'], df[df['x'] <= 10]['custom'], label='$J_0(X)$')
axarr[0].plot(df[df['x'] <= 10]['x'], df[df['x'] <= 10]['derivation'], label='$J_0\'(X)$')
axarr[0].grid()
axarr[0].legend()

axarr[1].set_title('Распределение ошибки вычисления')
df['Error'].plot.hist()
axarr[1].grid()

fig.savefig("bessel.png")
