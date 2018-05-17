import numpy as np
import math
import os
import matplotlib.pyplot as plt

def f(x):
    return np.array([x[1], -np.sin(x[0])])

def solve(N_steps):
    print(N_steps, 'точек')
    x_current = np.array([np.pi / 2, 0])
    T = 1.854074677301372 * 4
    step = T / N_steps
    l = list( range(N_steps + 1) )
    r = None
    result = []

    try:
        import tqdm
        r = tqdm.tqdm( l )
    except:
        r = l[:]

    for i in r:
        result.append([step * i] + list(x_current))

        x_1 = f( x_current )
        x_2 = f( x_1 * step / 2 + x_current )
        x_3 = f( x_2 * step / 2 + x_current )
        x_4 = f( x_current + step * x_3 )

        x_current = x_current + step / 6 * (x_1 + 2 * x_2 + 2 * x_3 + x_4)


    result = np.array(result)
    return result

if __name__ == "__main__":
    # сначала просто график решения
    result = solve(100)

    fig, ax = plt.subplots(2, figsize=(8, 9))
    ax[0].plot(result[:, 0],  result[:, 1], label='$\\varphi(t)$')
    ax[0].plot(result[:, 0], result[:, 2], label='$\\omega(t)$')
    ax[0].grid()
    ax[0].legend()

    # теперь зависимость ошибки от размера шага

    err = []

    for n in [10, 20, 30, 100, 500, 1000, 10000, 5000, 2000, 20000, 70, 700]:
        result = solve(n)
        err_phi = abs(result[0, 1] - result[-1, 1])
        err_omega = abs(result[0, 2] - result[-1, 2])
        err.append( [n, err_phi + err_omega] )

    err = np.array(err)

    ax[1].scatter( -np.log(err[:, 0]), np.log(err[:, 1]), label='Зависимость логарифма ошибки от логарифма размера шага' )
    ax[1].grid()

    from sklearn.linear_model import LinearRegression
    lr = LinearRegression()
    lr.fit( -np.log(err[:, 0]).reshape(-1, 1), np.log(err[:, 1]) )
    s = np.linspace(np.min(-np.log(err[:, 0])), np.max(-np.log(err[:, 0])), 100)
    r = lr.predict(s.reshape(-1, 1))
    ax[1].plot(s, r, label='$y~{:03f}x ~ {}x$'.format(lr.coef_[0], round(lr.coef_[0])))
    ax[1].legend()


    fig.savefig('plot.png')
