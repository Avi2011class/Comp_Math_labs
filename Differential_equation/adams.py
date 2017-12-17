import numpy as np
from math import *
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 10, 8

def f(t, u):
    ret = -2 * t * u + 2 * t
    print(t, u, ret)
    return ret

if __name__ == "__main__":
    tau = 0.01
    arr = [ [0, 1.00001], [tau, 2] ]
    
    for i in range(10):
        t1, u1 = arr[-1]
        t2, u2 = arr[-2]
        
        u_new = ( 3 * f(t1, u1) - f(t2, u2) ) / 2 * tau + u1
        t_new = t1 + tau
        
        arr.append((t_new, u_new))
    
    arr = np.array(arr)
    print(arr)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(arr[:, 0], arr[:, 1])
    ax.grid()
    fig.savefig('solution.png')
        

