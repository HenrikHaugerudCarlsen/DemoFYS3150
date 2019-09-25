import numpy as np
from Project2_jacobi import max_nondiag, rotate



def Test_eigvalues():

    tol_eig = 1e-10
    n = 5
    #opening and reading the txt file containing eigenvalues
    o_file = open(eigvalues_n.txt)
    if o_file.mode == 'r':
        contents = o_file.read()

    h = 1./(n+1)
    d = 2. / (h**2)
    a = -1. /(h**2)
    j = np.linspace(1,n,n)
    eig = d + 2*a*np.cos((j*np.pi)/(n+1))

    for i in range(n):
        success = abs(eig[i] - eig_A[i]) < tol_eig
        msg = 'Analytical eigenvalue = %.2f , computed eigenvalue = %.2f' %(eig[i], eig_A[i])

        assert success, msg

Test_eigvalues()
