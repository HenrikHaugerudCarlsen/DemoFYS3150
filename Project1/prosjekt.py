import numpy as np
import matplotlib.pyplot as plt
import time

def f_func(x):
    sol = 100*np.exp(-10*x)
    return sol

def closed_form(x):
    u = 1-(1-np.exp(-10))*x-np.exp(-10*x)
    return u

n_max = 7
i = np.linspace(1,n_max,n_max)

n = 10**i
h = 1/(n+1)

def general(n, h):
    n = int(n)

    x = np.linspace(0,1,n+2)

    f = f_func(x[1:-1])*(h**2)
    b = np.ones(n)*2
    a = np.ones(n-1)*-1
    c = np.ones(n-1)*-1

    b_g = np.zeros(n)
    f_g = np.zeros(n)
    u_g = np.zeros(n)

    b_g[0] = 2;f_g[0] = f[0]

    start = time.time()
    for i in range(n-1):
        b_g[i+1] = b[i+1] - a[i]*c[i]/b_g[i]
        f_g[i+1] = f[i+1] - a[i]*f_g[i]/b_g[i]

    u_g[-1] = f_g[-1]/b_g[-1]

    for i in range(n-2, -1, -1):
        u_g[i] = (f_g[i] - c[i]*u_g[i+1])/b_g[i]
    end = time.time()

    print('time spent for n = %g points with general algorithm was %.2e s' % (n, end-start))
    u_g = np.hstack([0, u_g, 0])
    return x, u_g

def optimized(n, h):
    n = int(n)
    i = np.linspace(1,n,n)

    x = np.linspace(0,1,n+2)
    f = f_func(x[1:-1])*(h**2)

    u_s = np.zeros(n)
    b_s = (i+1)/i
    f_s = np.zeros(n)
    f_s[0] = f[0]

    start = time.time()
    for i in range(n-1):
        f_s[i+1] = f[i+1] + f_s[i]/b_s[i]

    u_s[-1] = f_s[-1]/b_s[-1]

    for i in range(n-2,-1,-1):
        u_s[i] = (f_s[i] + u_s[i+1])/b_s[i]

    end = time.time()

    u_a = closed_form(x[1:-1])

    error_a = np.log10(np.abs((u_s-u_a)/u_a))

    print('time spent for n = %g points with optimized algorithm was %.7e s' % (n, end - start))
    u_s = np.hstack([0,u_s,0])
    return x, u_s, error_a, max(error_a)

j = 4
x, u_g = general(n[j], h[j])
x_s, u_s, error_a, error_max1 = optimized(n[j], h[j])

plt.title('numerical and analytical solution with n=%g' % n[j])
plt.plot(x, closed_form(x), label='analytic')
plt.plot(x, u_g, 'ro', label='general numeric')
plt.plot(x, u_s, label='special numeric')
plt.legend()
plt.show()

error_max = np.zeros(len(n))
fig,ax = plt.subplots(1)
for l in range(len(n)):
    x_s1, u_s, error, error_max[l] = optimized(n[l], h[l])
    ax.plot(x_s1[1:-1], error)
plt.show()

plt.plot(np.log10(h), error_max)
plt.xlabel('$log_{10}(h)$')
plt.ylabel('max relative error $\epsilon_i$')
plt.show()