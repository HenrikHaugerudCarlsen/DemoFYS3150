import numpy as np
import matplotlib.pyplot as plt
from Project2_jacobi import max_nondiag, rotate

N = 99
p0 = 0
p_max = 10
h = float(p_max)/(N+1)
j = np.linspace(1,N,N)

p = p0 + j*h

e = -1./(h**2)
d = 2/(h**2) + p**2

A = np.zeros((N,N))
A[0, 0] = d[0]; A[0, 1] = e; A[-1, -1] = d[-1]; A[-1, -2] = e

for i in range(1, N - 1):  # generating the tridiagonal matrix A
    A[i, i] = d[i]
    A[i, i + 1] = e
    A[i, i - 1] = e

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
print(idx)
plt.plot(np.hstack([0,p,p_max]), np.hstack([0,U[:,idx[0]],0]))
plt.show()

eig_A = np.sort(eig_A)


def analytical_eigv(x):
    lamb = 4*x-1
    return lamb

plt.plot(j, eig_A, 'r.')
plt.plot(j, analytical_eigv(j), 'g*')
plt.show()
plt.plot(j, abs(eig_A - analytical_eigv(j)))
plt.show()
