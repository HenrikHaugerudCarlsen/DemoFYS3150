import numpy as np, matplotlib.pyplot as plt
from Project2_jacobi import max_nondiag, rotate



#plot formatting
plt.style.use('ggplot')

font = {'family' : 'serif',
        'weight' : 'normal',
        'style'  : 'normal',
        'size'   : 12}
plt.rc('font', **font)

list = [10,50,100,200]

p0 = 0
p_max = 10

for N in list:
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

    tol = 1e-9
    iteration = 0
    max_iteration = N**3

    max = max_nondiag(A, N)[0]

    while max > tol and iteration < max_iteration:
        max, k, l = max_nondiag(A, N)
        A, U = rotate(A, U, k, l, N)
        iteration += 1
        #print(iteration, max_iteration)

    eig_A = np.zeros(N)
    for i in range(N):
        eig_A[i] = A[i, i]

    idx = np.where(eig_A == np.min(eig_A))[0]


    #print(idx)


    #plt.title('Wavefunction')

    #yolo = plt.plot(np.hstack([0,p,p_max]), np.hstack([0,U[:,idx[0]],0]), label = 'psi, N = %g' %(N))
    #plt.savefig('yolo.pdf')
    #plt.close()
#plt.legend()
#plt.show()

    eig_A = np.sort(eig_A)


    def analytical_eigv(x):
        lamb = 4*x-1
        return lamb

    #print(abs(eig_A - analytical_eigv(j)))



    plt.plot(j,abs(eig_A - analytical_eigv(j)), label = 'abs difference, N = %g' %(N))
    plt.legend()
    plt.title('|Analytical value - computed value|')
    plt.xlabel('Diagonal element i')
    plt.ylabel('Absolute difference')
plt.savefig('Erroroflambda.pdf')
plt.close()


#plt.show()



    #plt.title('Eigenvalues')
    #plt.xlabel('Diagonal element i')
    #plt.ylabel('Eigenvalue')

    #plt.plot(j, eig_A, 'r.', label = 'Computed value')
    #plt.legend()
    #plt.plot(j, analytical_eigv(j), 'g*', label = 'Analytical value')
    #plt.legend()

    #plt.show()
