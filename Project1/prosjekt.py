import numpy as np
import matplotlib.pyplot as plt
import time

def f_func(x):
    '''
    calculating the source term of the poisson equation
    '''
    sol = 100*np.exp(-10*x)
    return sol

def closed_form(x):
    '''
    analytical solution of the one dimensional poisson equation
    '''
    u = 1-(1-np.exp(-10))*x-np.exp(-10*x)
    return u

# generating the size n for the arrays and steplength h
n_max = 7
i = np.linspace(1,n_max,n_max)
n = 10**i
h = 1/(n+1)

def general(n, h):
    '''
    The general Thomas Algorithm for calculating the general solution of a tri-diagonal matrix equation
    '''
    n = int(n)

    #initial arrays b, a, c, x and f
    x = np.linspace(0, 1, n + 2)
    f = f_func(x[1:-1])*(h**2)
    b = np.ones(n)*2
    a = np.ones(n-1)*-1
    c = np.ones(n-1)*-1

    #initializing the arrays used in the forward and backwards substitution
    b_g = np.zeros(n)
    f_g = np.zeros(n)
    u_g = np.zeros(n)

    b_g[0] = 2;f_g[0] = f[0]

    start = time.time()
    for i in range(n-1):  # algorithm for the forward substitution
        b_g[i+1] = b[i+1] - a[i]*c[i]/b_g[i]
        f_g[i+1] = f[i+1] - a[i]*f_g[i]/b_g[i]

    u_g[-1] = f_g[-1]/b_g[-1]

    for i in range(n-2, -1, -1):  # algorithm for the backwards substitution
        u_g[i] = (f_g[i] - c[i]*u_g[i+1])/b_g[i]
    end = time.time()

    print('time spent for n = %g points with general algorithm was %.2e s' % (n, end-start))
    u_g = np.hstack([0, u_g, 0])  # adding the boundary conditions on the general solution
    return x, u_g

def optimized(n, h):
    '''
    The optimized Algorithm for calculating the specific solution for our tri-diagonal matrix equation
    '''
    N = int(n)
    i = np.linspace(1, N, N)

    # initialize the known arrays x and f
    x = np.linspace(0, 1, N+2)
    f = f_func(x[1:-1])*(h**2)

    # generating the arrays b, u and f for the special case
    u_s = np.zeros(N)
    b_s = (i+1)/i
    f_s = np.zeros(N)
    f_s[0] = f[0]

    start = time.time()
    for i in range(N-1):  # optimized algorithm for the forward substitution
        f_s[i+1] = f[i+1] + f_s[i]/b_s[i]

    u_s[-1] = f_s[-1]/b_s[-1]

    for i in range(N-2, -1, -1):  # optimized algorithm for the backwards substitution
        u_s[i] = (f_s[i] + u_s[i+1])/b_s[i]
    end = time.time()

    print('time spent for n = %g points with optimized algorithm was %.7e s' % (n, end - start))

    u_s = np.hstack([0, u_s, 0]) # adding the boundary conditions on the optimized solution

    # finding the relative error of the algorithm
    u_a = closed_form(x[1:-1])
    error_a = np.log10(np.abs((u_s-u_a)/u_a))

    return x, u_s, error_a, max(error_a)

# comparing the time spent with the general and optimized algorithm with n=10^j points
j = 4
x_g, ug = general(n[j], h[j])
x_s, us, error_a, error_max1 = optimized(n[j], h[j])

plt.title('numerical and analytical solution with n=%g' % n[j])
plt.plot(x_g, closed_form(x_g), label='analytic')
plt.plot(x_g, ug, 'ro', label='general numeric')
plt.plot(x_s, us, label='special numeric')
plt.legend()
plt.show()

# finding the max error of the optimized algorithm and plots it as a function of the briggsian logarithm of the
# steplength h
error_max = np.zeros(len(n))
fig, ax = plt.subplots(1)
for l in range(len(n)):
    x_s1, u_s, error, error_max[l] = optimized(n[l], h[l])
    ax.plot(x_s1[1:-1], error)
plt.show()

plt.plot(np.log10(h), error_max)
plt.xlabel('$log_{10}(h)$')
plt.ylabel('max relative error $\epsilon_i$')
plt.show()