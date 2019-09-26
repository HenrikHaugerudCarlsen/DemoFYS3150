import numpy as np
from Project2_jacobi import max_nondiag, rotate #importing functions from other project files



def Test_eigvalues(): # making testfunction for testing computed eigenvalues against analytial ones

    tol_eig = 1e-12
    n = 5
    #opening and reading the txt file containing eigenvalues and putting eigenvalues into array
    eig_A = np.zeros(n)
    o_file = open('eigvalues_5.txt', 'r')
    o1 = o_file.readlines()
    for i in range(len(o1)):
        eig_A[i] = float(o1[i])

    h = 1./(n+1)
    d = 2. / (h**2)
    a = -1. /(h**2)
    j = np.linspace(1,n,n)
    eig = d + 2*a*np.cos((j*np.pi)/(n+1)) # analytical solution

    for i in range(n): #success statement
        success = abs(eig[i] - eig_A[i]) < tol_eig
        #print(abs(eig[i]-eig_A[i]))
        msg = 'Analytical eigenvalue = %.2f , computed eigenvalue = %.2f' %(eig[i], eig_A[i])

        assert success, msg

Test_eigvalues() # calling testfuntion

def Test_maxnondiag(): #making testfunction for testing that the right max non_diagonal elements are found
    tol_diag = 1e-12
    n = 5
    A = np.zeros((n,n))
    for i in range(n-1):
        for j in range(i+1,n):
            A[i,j] = np.random.randint(1,25) #generation random upper diagonal matrix
            A[j,i] = A[i,j] # generating random symmetrical matrix
    #print(A)
    max = np.max(A)
    k,l = np.where(A == max
    max_c, k_c, l_c = max_nondiag(A, n)
    #print(max, max_c, k, k_c, l, l_c)

    success1 = abs(max-max_c) < tol_diag
    success2 = k_c < l_c and k_c == k[0] and l_c == l[0] #checking that we are in the elements above the diagonal,
    #and that the analytical position matches the computed position of the max value element

    msg1 = 'max nondiag element of A = %g and computed max = %g' % (max, max_c)
    msg2 = 'analytical: [k,l] = [%g , %g] , computed: [k,l] = [%g, %g] ' %(k[0],l[0],k_c,l_c)

    assert success1, msg1
    assert success2, msg2


Test_maxnondiag()
