# -*- coding: utf-8 -*-

import numpy, math
import os, sys


class pylib(object):
    """
    Implements many of the functions found in M.H.Jensens
    c++-library, used in Computational Physics. This is again heavily based
    on what is found in Numerical recipes.

    Ported to Python by
    Kyrre Ness Sjøbæk      (k DÅTT n DÅTT sjobak ÆTT fys DÅTT uio DÅTT no),
    Magnar Kopangen Bugge  (magnarkb ÆTT student DÅTT matnat DÅTT uio DÅTT no),
    Marit Sandstad         (marit DÅTT sandstad ÆTT fys DÅTT uio DÅTT no)
    """

    ZERO = 1.0E-10
    """Used as a meassure of machine precision in some algos"""

    yTol = None
    """Vector or number giving the tolerance for the stepsize control functions"""
    guard = 0.95
    """Used for the adaptive stepsize control functions. Standard value (0.95) is just a guess..."""

    def __init__(self, inputcheck=False, cpp=False):
        """
        Constructor,
        Set inputcheck = True in order to do input checking (debug, slower)
        Set cpp = True to use compiled versions of the functions where aviable
        """
        self.inputcheck = inputcheck
        self.cpp = cpp
        if cpp:
            sys.path.append(os.path.abspath(".") + "/cpp")
            global pylib_cpp
            import pylib_cpp

    def luDecomp(self, A):
        """
        LU-decomposes a matrix A, and returns the LU-decomposition of a
        rowwise permutation of A. Used in combination with luBackSubst function
        to solve an equation-set Ax=B.
        Returns: A tuple containing
        - LU-decomposed matrix (upper part and diagonal = matrix U, lower part = matrix L)
        - Array which records the row permutation from partial pivoting
        - Number which depends on the number of row interchanges was even (+1) or odd (-1)
        BIG FAT WARNING: Destroys input A in calling program!
        (A is set equal to the returned LU-decomposed matrix)
        Send it a copy if this is bad.
        This function has the ability to switch between Python and C++ backends, see __init__()
        """

        if self.inputcheck:
            self.checkSquare(A)

        d = 1.0  # Records row interchange (parity)
        N = A.shape[0]
        index = numpy.zeros(N, numpy.int32)

        if self.cpp:
            return self.luDecomp_cpp(A, N, index, d)
        else:
            return self.luDecomp_python(A, N, index, d)

    def luDecomp_cpp(self, A, N, index, d):
        """
        C++ backend for luDecomp, using routine
        in the library
        """
        pylib_cpp.ludcmp(A, index, d)
        return (A, index, d)

    def luDecomp_python(self, A, N, index, d):
        """
        Python backend for luDecomp
        """

        # Loop over rows to get scaling info, check if matrix is singular
        vv = numpy.zeros(N)
        for i in range(N):
            big = numpy.fabs(A[i]).max()
            if big < self.ZERO:
                raise SingularError(A)
            vv[i] = 1.0 / big

        # Loop over columns, Crout's method
        for j in range(N):
            # i < j
            for i in range(j):
                sum = A.item(i, j)
                for k in range(i):
                    sum = sum - A.item(i, k) * A.item(k, j)
                A.itemset(i, j, sum)

            # i >= j
            imax = None
            big = self.ZERO
            for i in range(j, N):
                sum = A.item(i, j)
                for k in range(j):
                    sum = sum - A.item(i, k) * A.item(k, j)
                A.itemset(i, j, sum)
                # Find biggest entry in vv (scaled by THIS sum)
                dum = vv[i] * math.fabs(sum)
                if dum >= big:
                    big = dum
                    imax = i
            # Do we need to interchange rows?
            if j != imax:
                # Oui!
                dum2 = A[imax].copy()
                A[imax] = A[j]
                A[j] = dum2
                d = d * (-1)
                vv[imax] = vv[j]
            index[j] = imax

            # If the pivot element is to small, the matrix is singular
            # to working order. For some applications of singular matrices,
            # it is desirable to substitute in self.ZERO
            if math.fabs(A.item(j, j) < self.ZERO):
                A.itemset(j, j, self.ZERO)

            # Divide by pivot element
            if j < (N - 1):
                A[(j + 1):N, j] = A[(j + 1):N, j] / A.item(j, j)

        return (A, index, d)

    def luBackSubst(self, A, index, b):
        """
        Back-substitution of LU-decomposed matrix
        Solves the set of linear equations A x = b of dimension n
        The input A is the LU-decomposed version of A obtained from
        pylib.luDecomp(),
        index is the pivoting permutation vector obtained
        from the same function, and
        b is the right-hand-side vector b as a numpy array.
        Returns the solution x as an numpy array.
        BIG FAT WARNING: Destroys input b in calling program!
        (b is set equal to x after calculation has finished)
        Send it a copy if this is bad.
        """

        if self.inputcheck:
            self.checkSquare(A)
            if A.shape[0] != b.shape[0]:
                raise SizeError("A:(%d,%d), b:(%d)" % (A.shape[0], A.shape[1], b.shape[0]))

        N = A.shape[0]

        if self.cpp:
            return self.luBackSubst_cpp(A, N, index, b)
        else:
            return self.luBackSubst_python(A, N, index, b)

    def luBackSubst_cpp(self, A, N, index, b):
        """
        C++ backend for luDecomp, using routine
        in the library
        """
        pylib_cpp.luBackSubst(A, index, b)
        return b

    def luBackSubst_python(self, A, N, index, b):
        """
        Python backend for luBackSubst
        """

        ii = -1
        for i in range(N):
            ip = index[i]
            sum = b[ip]
            b[ip] = b[i]
            if ii > -1:
                for j in range(ii, i):
                    sum = sum - A[i, j] * b[j]
            elif sum != 0:
                ii = i
            b[i] = sum
        for i in range((N - 1), -1, -1):
            sum = b[i]
            for j in range(i + 1, N):
                sum = sum - A[i, j] * b[j]
            b[i] = sum / A[i, i]

        return b