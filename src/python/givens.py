# This script is for testing of calculating the determinant

from math import copysign,sqrt,log
import numpy as np
from scipy.linalg import eigvals

def prettyprint(M):
    dim = len(M)
    for i in range(dim):
        for j in range(dim):
            print("%+.3g " % M[i][j], end="")
        print()


def swap(M,i,j):
    dim = len(M)

    # swap columns
    for k in range(dim):
        M[k][i],M[k][j] = M[k][j],M[k][i]

    # swap rows
    for k in range(dim):
        M[i][k],M[j][k] = M[j][k],M[i][k]


def pivot(M):
    dim = len(M)
    for k in range(dim):
        # find maximum
        index = k 
        for z in range(k+1,dim):
            if abs(M[z][z]) < abs(M[index][index]):
                index = z 

        if k != index:
            swap(M,k,index)


def QR(M,dim):
    """
    Perform QR decomposition of real matrix M using Givens rotations. The
    matrix M will be decomposed to M = QR, where Q is an orthogonal matrix and
    R is an upper triangular matrix. This method will not compute Q.
    On exit: M = R
    See: https://en.wikipedia.org/wiki/Givens_rotation
    """

    for j in range(dim-1):
        for i in range(j+1,dim):
            a = M[j][j]
            b = M[i][j]

            if abs(b) > abs(a):
                t = a/b
                u = copysign(sqrt(1+t**2),b)
                s = -1/u
                c = -s*t
                r = b*u
            else:
                t = b/a
                u = copysign(sqrt(1+t**2),a)
                c = 1/u
                s = -c*t
                r = a*u

                
            for n in range(j,dim):
                Min = M[i][n]
                Mjn = M[j][n]

                M[i][n] = c*Min + s*Mjn
                M[j][n] = c*Mjn - s*Min

            M[i][j] = 0


def logdet(M):
    dim = len(M)
    QR(M,dim)
    
    det = 0
    for i in range(dim):
        det += log(abs(M[i][i]))
    return det


def load(filename,dim):
    M = []

    for i in range(dim):
        M.append([0]*dim)

    with open(filename) as f:
        for line in f:
            i,j,Mij = line.split(",")
            i,j = int(i), int(j)
            M[i][j] = float(Mij)

    return M


if __name__ == "__main__":
    M = [
        [100,2,3,4],
        [5,5,7,8],
        [9,10,11,12],
        [13,14,15,16],
    ]

    pivot(M)
    prettyprint(M)
    print(logdet(M))
    print()


    #M = load("testmatrix", dim)
    #pivot(M)

    #print(logdet(M,dim))
