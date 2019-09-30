import numpy as np, matplotlib.pyplot as plt
from Project2_jacobi import max_nondiag, rotate

# plot formatting
plt.style.use('ggplot')

font = {'family' : 'serif',
        'weight' : 'normal',
        'style'  : 'normal',
        'size'   : 12}
plt.rc('font', **font)

N = np.array([10,50,100,200])

p_m = 10

def error_eigA(N, p):
    n = int(N)
    p0 = 0
    p_max = p
    h = float(p_max) / (n + 1)
    j = np.linspace(1, n, n)

    p = p0 + j * h

    e = -1. / (h ** 2)
    d = 2 / (h ** 2) + p ** 2

    A = np.zeros((n, n))
    A[-1, -1] = d[-1]

    for i in range(n - 1):  # generating the tridiagonal matrix A
        A[i, i] = d[i]
        A[i, i + 1] = e
        A[i + 1, i] = e

    U = np.zeros((n, n))
    for i in range(n):
        U[i, i] = 1

    tol = 1e-8
    iteration = 0
    max_iteration = n ** 3

    max = max_nondiag(A, n)[0]

    while max > tol and iteration < max_iteration:
        max, k, l = max_nondiag(A, n)
        A, U = rotate(A, U, k, l, n)
        iteration += 1

    eig_A = np.zeros(n)
    for i in range(n):
        eig_A[i] = A[i, i]

    analytical_eig = 3
    err = abs(np.min(eig_A)-analytical_eig)

    return err, eig_A

def analytical_eigv(x):
    lamb = 4*x-1
    return lamb

number_N = len(N)
error = np.zeros(number_N)

for i in range(number_N):
    n = N[i]
    error[i], eig_A = error_eigA(n, p_m)
    eig_A = np.sort(eig_A)
    j = np.linspace(1, n, n)
    plt.plot(j,abs(eig_A - analytical_eigv(j)), label = 'abs difference, N = %g' % n)

plt.title('|Analytical value - computed value|')
plt.xlabel('Diagonal element i')
plt.ylabel('Absolute difference')
plt.savefig('Erroroflambda.pdf')
plt.legend()
plt.close()

error_tol = np.full(number_N, 1e-3)

plt.title(r'error as a function of N, $\rho_{max} = 10$')
plt.plot(N, error, 'r--', label=r'$\epsilon(N)$')
plt.plot(N, error_tol, 'b-', label=r'error toleration of $\epsilon_i < 10^{-3}$')
plt.yscale('log')
plt.xlabel('integration points, N')
plt.ylabel(r'error of $\lambda_0$')
plt.legend()
plt.savefig('errorfunc.pdf')
plt.show()
