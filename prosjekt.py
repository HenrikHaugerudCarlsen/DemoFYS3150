import numpy as np
import matplotlib.pyplot as plt
import time

def f_func(x):
    sol = 100*np.exp(-10*x)
    return sol

def closed_form(x):
    u = 1-(1-np.exp(-10))*x-np.exp(-10*x)
    return u

n = 10
h = 1/(n+1)

u_g = np.zeros(n+2)
x = np.linspace(0,1, n+2)
f = f_func(x)*(h**2)
b = np.ones(n)*2
a = np.ones(n-1)*-1
c = np.ones(n-1)*-1

b_g = np.zeros(n)
f_g = np.zeros(n)

b_g[0] = 2;f_g[0] = f[1]

start1 = time.time()
for i in range(n-1):
    b_g[i+1] = b[i+1] - a[i]*c[i]/b_g[i]
    f_g[i+1] = f[i+2] - a[i]*f_g[i]/b_g[i]

u_g[n] = f_g[-1]

for i in range(n-1):
    u_g[n-1-i] = f_g[-2-i] - c[-1-i]*u_g[n-i]/(b_g[-2-i])
end1 = time.time()

i = np.linspace(1,n,n)

u_s = np.zeros(n+2)
b_s = (i+1)/i
f_s = np.zeros(n)
f_s[0] = f[1]

start2 = time.time()
for i in range(n-1):
    f_s[i+1] = f[i+2] + f_s[i]/b_s[i]

u_s[n] = f_s[-1]

for i in range(n-1):
    u_s[n-1-i] = f_s[-2-i] + u_s[n-i]/b_g[-2-i]
end2 = time.time()

print(end1-start1, end2-start2)
print(end1, start1)

plt.plot(x, closed_form(x))
plt.plot(x, u_g)
plt.plot(x, u_s)
plt.show()
