import numpy as np
import matplotlib.pyplot as plt
from Project2_jacobi import max_nondiag, rotate

N = 100
p0 = 0
p_max = 12
h = float(p_max)/(N+1)
j = np.linspace(1,N,N)
w_r = np.array([0.01, 0.5, 1, 5])

p = p0 + j*h
def plot(w_r, ind):
    e = -1./(h**2)
    d = 2/(h**2) + w_r**2 * p**2 #+ (1/p) repulsion term

    A = np.zeros((N,N))
    A[-1, -1] = d[-1]

    for i in range(N - 1):  # generating the tridiagonal matrix A
        A[i, i] = d[i]
        A[i, i + 1] = e
        A[i+1, i] = e

    U = np.zeros((N, N))
    for i in range(N):
        U[i, i] = 1

    tol = 1e-8
    iteration = 0
    max_iteration = N**3

    max = max_nondiag(A, N)[0]

    while max > tol and iteration < max_iteration:
        max, k, l = max_nondiag(A, N)
        A, U = rotate(A, U, k, l, N)
        iteration += 1

    eig_A = np.zeros(N)
    for i in range(N):
        eig_A[i] = A[i, i]

    idx = np.where(eig_A == np.min(eig_A))[0]

    plt.subplot(4,1,ind+1)
    plt.plot(np.hstack([0,p,p_max]), np.hstack([0,U[:,idx[0]],0]), label=r'groundstate for $\omega_r = %g$' % w_r)
    plt.xlabel(r'$\rho$')
    plt.ylabel(r'$u(\rho)$')


def analytical_eigv(x):
    lamb = 4*x-1
    return lamb

for i in range(len(w_r)):
    plot(w_r[i], i)
    plt.legend()

plt.suptitle(r'the groundstate for different frequencies $\omega_r$')
plt.savefig('groundstates.pdf')
plt.show()