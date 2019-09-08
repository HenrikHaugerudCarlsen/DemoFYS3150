from LU import pylib
import numpy as np
import time

def f_func(x):
    '''
    calculating the source term of the poisson equation
    '''
    sol = 100*np.exp(-10*x)
    return sol

n_max = 4
i = np.linspace(1, n_max, n_max)
N = 10**i
#N_test = 10**5
h = 1/(N+1)

for i in range(len(N)):
    n = int(N[i])
    x = np.linspace(0,1,n+2)
    A = np.zeros((n, n))
    A[0,0] = 2; A[0,1] = -1; A[-1,-1] = 2; A[-1,-2] = -1
    b = f_func(x[1:-1])*(h[i]**2)

    for i in range(1, n-1): # generating the tridiagonal matrix A
        for j in range(1, n-1):
            if i == j:
                A[i, j] = 2
                A[i, j+1] = -1
                A[i, j-1] = -1

    start = time.time()
    pl = pylib() # using the library functions from Github to calculate the LU decomposistion and solve the equation
    LU, index, d = pl.luDecomp(A)
    u = pl.luBackSubst(LU, index, b)
    end = time.time()

    print('time taken to compute the LU-decomposition and the solution is %.4f seconds with n = %g' % ((end-start), n))