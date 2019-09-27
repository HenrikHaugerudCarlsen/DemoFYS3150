import numpy as np, matplotlib.pyplot as plt
plt.style.use('ggplot')

font = {'family' : 'serif',
        'weight' : 'normal',
        'style'  : 'normal',
        'size'   : 10}
plt.rc('font', **font)

#Jacobis method

def max_nondiag(A, n): #defining function for finding the max nondiagonal element
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

def rotate(A, U, k, l, n): # defining function for unitary transforms of the matrix elements
    if A[k,l] != 0: #checking if elements are non-zero
        #defining tau, t, c and s
        tau = (A[l,l]-A[k,k])/(2*A[k,l])

        if tau > 0:
            t = 1. / (tau + np.sqrt(1. + tau**2))
        else:
            t = 1. / (tau - np.sqrt(1. + tau*tau))

        c = 1 / np.sqrt(1 + t**2)
        s = t*c

    else: # if elements are not non-zero:
        c = 1; s = 0
#Setting element-values as constants to reduce memory reads
    a_kk = A[k,k]
    a_ll = A[l,l]
    a_kl = A[k,l]
#Updating elements with algorithm
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

if __name__ == '__main__':
    # defining constants
    tol = 1e-8
    n = 5
    h = 1. / (n + 1)
    d = 2. / (h ** 2)
    a = -1. / (h ** 2)
    max_iteration = n ** 3
    iteration = 0

    # setting matrix elements
    A = np.zeros((n, n))
    A[0, 0] = d; A[0, 1] = a; A[-1, -1] = d; A[-1, -2] = a

    for i in range(1, n - 1):  # generating the tridiagonal matrix A
        A[i, i] = d
        A[i, i + 1] = a
        A[i, i - 1] = a

    # generating eigenvector matrix
    U = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                U[i, j] = 1

    max = max_nondiag(A, n)[0]

    while max > tol and iteration < max_iteration:
        max, k, l = max_nondiag(A, n)
        A, U = rotate(A, U, k, l, n)
        iteration += 1

    eig_A = np.zeros(n)
    for i in range(n):
        eig_A[i] = A[i, i]
    eig_A = np.sort(eig_A)

    f_n = open('eigvalues_%g.txt' % n, 'w+')

    for i in range(n):
        f_n.write('%.12f\n' % eig_A[i])

    p = np.linspace(0,1,n+2)

    plt.plot(p, np.hstack([0,U[:,0],0]), label='solution with eigenvalue %.2f' % A[0, 0])

    plt.legend()
    plt.show()
