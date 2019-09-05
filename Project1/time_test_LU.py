from LU import pylib
import numpy as np

N = 5
A = np.zeros((N,N))
A[0,0] = 2; A[0,1] = -1; A[-1,-1] = 2; A[-1,-2] = -1

for i in range(1, N-1):
    for j in range(1, N-1):
        if i == j:
            A[i, j] = 2
            A[i, j+1] = -1
            A[i, j-1] = -1

pl = pylib()
print(pl.luDecomp(A))
