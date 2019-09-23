import numpy as np, matplotlib.pyplot as plt

#defining constants
tol = 1e-7
n = 10
h = 1./(n+1)
d = 2./(h**2)
a = -1./(h**2)
max_iteration = n**2
iteration = 0

#setting matrix elements
A = np.zeros((n, n))
A[0,0] = d; A[0,1] = a; A[-1,-1] = d; A[-1,-2] = a


for i in range(1, n-1): # generating the tridiagonal matrix A
    for j in range(1, n-1):
        if i == j:
            A[i, j] = d
            A[i, j+1] = a
            A[i, j-1] = a

#generating eigenvector matrix
U = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if i == j:
            U[i,j] = 1

print(A)
eig = np.linalg.eigh(A)
print(eig[0])

j = np.linspace(1,n,n)
print(j)
lam = d + 2*a*np.cos((j*np.pi)/(n+1))
print(lam)

#Jacobis method

#tau = np.zeros(n-1)
def max_nondiag(A, n):
    max = 0
    k = 0
    l = 0
    for i in range(0, n):
        for j in range(i+1, n):
            if abs(A[i,j]) > max:
                max = abs(A[i,j])
                k = i
                l = j
    return max, k, l
max, k, l = max_nondiag(A,n)
print(k,l)

def rotate(A, U, k, l, n):
    if A[k,l] != 0:
        tau = (A[l,l]-A[k,k])/(2*A[k,l])

        if tau > 0:
            t = 1. / (-tau + np.sqrt(1. + tau*tau))
        else:
            t = -1. / (-tau + np.sqrt(1. + tau*tau))

        c = 1 / np.sqrt(1 + t**2)
        s = t*c

    else :
        c = 1; s = 0

    a_kk = A[k,k]
    a_ll = A[l,l]
    a_kl = A[k,l]
    A[k,k] = a_kk*c**2 - 2*a_kl*c*s + a_ll*s**2
    A[l,l] = a_ll*c**2 + 2*a_kl*c*s + a_kk*s**2
    A[k,l] = 0
    A[l,k] = 0

    for i in range(n):
        if i != k and i != l:
            a_ik = A[i,k]
            a_il = A[i,l]

            A[i,k] = a_ik*c - a_il*s
            A[k,i] = A[i,k]
            A[i,l] = a_il*c + a_ik*s
            A[l,i] = A[i,l]

        u_ik = U[i,k]
        u_il = U[i,l]
        U[i,k] = u_ik*c - u_il*s
        U[i,l] = u_il*c + u_ik*s

    return A, U

max = max_nondiag(A, n)[0]

while max > tol and iteration < max_iteration:
    max, k, l = max_nondiag(A,n)
    A, U = rotate(A, U, k, l, n)
    iteration += 1
    #print(iteration, k, max)

for i in range(n):
    print(A[i,i])
