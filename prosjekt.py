import numpy as np
import matplotlib.pyplot as plt
import time

def f_func(x):
    sol = 100*np.exp(-10*x)
    return sol

def closed_form(x):
    u = 1-(1-np.exp(-10))*x-np.exp(-10*x)
    return u

i = np.linspace(1,6,6)

n = 10**i
h = 1/(n+1)

def general(n, h):
    n = int(n)

    x = np.linspace(0,1,n+2)

    f = f_func(x)*(h**2)
    b = np.ones(n)*2
    a = np.ones(n-1)*-1
    c = np.ones(n-1)*-1

    b_g = np.zeros(n)
    f_g = np.zeros(n)
    u_g = np.zeros(n+2)

    b_g[0] = 2;f_g[0] = f[1]

    start = time.time()
    for i in range(n-1):
        b_g[i+1] = b[i+1] - a[i]*c[i]/b_g[i]
        f_g[i+1] = f[i+2] - a[i]*f_g[i]/b_g[i]

    u_g[n] = f_g[-1]/b_g[-1]

    for i in range(n-1):
        u_g[n-1-i] = f_g[-2-i] - c[-1-i]*u_g[n-i]/(b_g[-2-i])
    end = time.time()

    print('time spent for n = %g points with general algorithm was %.2e s' % (n, end-start))
    return x, u_g

def optimized(n, h):
    n = int(n)
    i = np.linspace(1,n,n)

    x = np.linspace(0,1,n+2)
    f = f_func(x)*(h**2)

    u_s = np.zeros(n+2)
    b_s = (i+1)/i
    f_s = np.zeros(n)
    f_s[0] = f[1]

    start = time.time()
    for i in range(n-1):
        f_s[i+1] = f[i+2] + f_s[i]/b_s[i]

    u_s[n] = f_s[-1]/b_s[-1]

    for i in range(n-1):
        u_s[n-1-i] = f_s[-2-i] + u_s[n-i]/b_s[-2-i]
    end = time.time()

    error_a = np.zeros(n)
    u = closed_form(x)

    for j in range(n):
        error_a[j] = np.log10(abs((u_s[j+1]-u[j+1])/u[j+1]))

    print('time spent for n = %g points with optimized algorithm was %.2e s' % (n, end - start))

    return x, u_s, error_a, max(error_a)

x, u_g = general(n[1], h[1])
error_max = np.zeros(len(n))
for l in range(len(n)):
    x_s, u_s, error, error_max[l] = optimized(n[l], h[l])

plt.plot(x_s, closed_form(x_s), label='analytic')
plt.plot(x, u_g, label='general numeric')
plt.plot(x_s, u_s, label='special numeric')
plt.legend()
plt.show()

plt.plot(x_s[1:-1], error)
plt.show()

plt.plot(np.log10(h), error_max)
plt.show()