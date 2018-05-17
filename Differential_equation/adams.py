import numpy as np
from math import *
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 10, 8

def f(t, u):
    ret = -2 * t * u + 2 * t
    return ret

def solve(tau, N):
    arr = [ [0, 2], [tau, 2 + tau * f(0, 2)] ]
    for i in range(N):
        t1, u1 = arr[-1]
        t2, u2 = arr[-2]
        
        u_new = ( 3 * f(t1, u1) - f(t2, u2) ) / 2 * tau + u1
        t_new = t1 + tau
        
        arr.append((t_new, u_new))
    
    arr = np.array(arr)
    return arr[-1][1]

if __name__ == "__main__":
    tau = 0.1
    arr = [ [0, 2], [tau, 2 + tau * f(0, 2)] ]
    N = 20
    for i in range(N):
        t1, u1 = arr[-1]
        t2, u2 = arr[-2]
        
        u_new = ( 3 * f(t1, u1) - f(t2, u2) ) / 2 * tau + u1
        t_new = t1 + tau
        
        arr.append((t_new, u_new))
    
    arr = np.array(arr)
    
    fig, ax = plt.subplots(2)
    
    l = np.linspace(0, N * tau, 1000)
    y = 1 + np.exp(-l**2)
    ax[0].plot(l, y, linewidth=15, alpha=0.5)
    ax[0].plot(arr[:, 0], arr[:, 1], linewidth=3, c='red')
    ax[0].grid()
    
   
    
    err = []
    x_fin = 100
    
    N = 200
    
    ll = np.linspace(0.0001, 0.001, 10)
    for tau in ll:
        N = int(x_fin / tau)
        err.append(solve(tau, N) - 1 - np.exp( -(tau * N) ** 2 ))
        
    ax[1].scatter(ll, ll)
    ax[1].grid()    
    fig.savefig('solution.png')
    print(err)

