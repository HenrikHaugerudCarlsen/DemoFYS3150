import matplotlib.pyplot as plt, numpy as np, random as random
"""
A = np.array(([1.+0.j , 0.-2.j],[0.+2.j , 5.+0.j]))

eig = np.linalg.eig(A) # displays the complex part of eigenvalue
eigh = np.linalg.eigh(A) #displays only the real part of eigenvalue
#print(eig)
#print(eigh)

diag = np.diag(A)
print(diag)
"""
#trying to generate an nxn matrix
N = 3
nxn = np.random.rand(N,N)
#print(nxn)
#print(np.diag(nxn))
print(np.linalg.eig(nxn))
print(np.linalg.eigh(nxn))
