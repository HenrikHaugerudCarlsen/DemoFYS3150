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

def solution(N):
    tol = 1e-8
    n = int(N)
    h = 1. / (n + 1)
    d = 2. / (h ** 2)
    a = -1. / (h ** 2)
    max_iteration = n ** 3
    iteration = 0

    # setting matrix elements
    A = np.zeros((n, n))
    A[-1, -1] = d;

    for i in range(n - 1):  # generating the tridiagonal matrix A
        A[i, i] = d
        A[i, i + 1] = a
        A[i + 1, i] = a

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

    return U, eig_A, iteration

if __name__ == '__main__':
    n_min = 10
    n_max = 100
    j = 10
    n = np.linspace(n_min, n_max, j)

    for i in range(j):
        iteration = solution(n[i])[-1]
        print('number of iterations is %g for n = %g' % (iteration, n[i]))

    U, eig_A, iteration = solution(n_max)

    f_n = open('eigvalues_%g.txt' % n_max, 'w+')

    for i in range(n_max):
        f_n.write('%.12f\n' % eig_A[i])

    p = np.linspace(0,1,n_max+2)

    idx1 = np.where(eig_A == np.min(eig_A))[0]
    max2 = np.max(eig_A)
    idx2 = 0
    for i in range(len(eig_A)):
        if eig_A[i] != np.min(eig_A):
            if eig_A[i] < max2:
                max2 = eig_A[i]
                idx2 = i

    max3 = np.max(eig_A)
    idx3 = 0
    for i in range(len(eig_A)):
        if i != idx2 and i != idx1[0]:
            if eig_A[i] < max3:
                max3 = eig_A[i]
                idx3 = i

    plt.title(r'Displacement of beam for $\lambda_0$ = %.2f, $\lambda_1$ = %.2f and $\lambda_2$ = %.2f'\
    % (eig_A[idx1[0]], eig_A[idx2], eig_A[idx3]))
    figure = plt.plot(p, np.hstack([0,U[:,idx2],0]), label=r'u($\rho$), $\lambda_1$ = %.2f' % eig_A[idx2])
    plt.plot(p, np.hstack([0,U[:,idx1[0]],0]), label=r'u($\rho$), $\lambda_0$ = %.2f' % eig_A[idx1[0]])
    plt.plot(p, np.hstack([0,U[:,idx3],0]), label=r'u($\rho$), $\lambda_2$ = %.2f' % eig_A[idx3])
    plt.xlabel(r'$\rho$')
    plt.ylabel(r'u($\rho$)')
    plt.legend()
    plt.savefig('BucklingBeamlambda.pdf')
    plt.close()
    plt.show()
