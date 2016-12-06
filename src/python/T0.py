import numpy as np
import libcasimir
from multiprocessing import Pool
from multiprocessing.context import TimeoutError
from time import sleep
from math import fsum

def logdetD(LbyR,xi,m,lmax):
    """Calculate log(det(1-M))"""
    casimir = libcasimir.Casimir(LbyR,xi,lmax=lmax)
    return casimir.logdetD(1,m)

class Queue:
    def __init__(self,LbyR,lmax,processes=4):
        self.LbyR = LbyR
        self.lmax = lmax
        self.queue = {}
        self.pool = Pool(processes=processes)

    def submit(self,xi,m):
        LbyR, lmax = self.LbyR, self.lmax
        self.queue[(xi,m)] = self.pool.apply_async(logdetD, (LbyR,xi,m,lmax))

    def join(self):
        d = []
        for key in list(self.queue.keys()):
            try:
                v = self.queue[key].get(False)
                del self.queue[key]
                xi,m = key
                d.append((xi,m,v))
            except TimeoutError:
                pass

        return d

    def __len__(self):
        return len(self.queue)


def accurate(l, eps):
    """Determine if we have to calculate more terms for higher values of m"""
    l0 = l[0]
    
    if l0 > 0:
        return False

    if l0 == 0:
        return True

    for j in range(1,len(l)):
        lj = l[j]
        if lj > 0:
            return False
        if abs(lj/l0) < eps:
            return True

    return True


if __name__ == "__main__":
    prec = 1e-12
    LbyR = 0.1
    processes = 6
    deg = 80
    lmax = 80
    idle = 2 # in ms

    alpha = 2*LbyR/(1+LbyR);
    xk,wk = np.polynomial.laguerre.laggauss(deg)

    d = {}
    for i,x in enumerate(xk):
        d[x/alpha] = i

    queue = Queue(LbyR, lmax)

    # initialize results
    M = np.empty((deg,lmax))
    M.fill(1)

    for m in range(lmax):
        stop = True

        for j,xj in enumerate(xk):
            if not accurate(M[j,:],prec):
                stop = False
                xi = xj/alpha
                queue.submit(xi,m)

            while len(queue) >= processes:
                for xi_,m_,v in queue.join():
                    M[d[xi_],m_] = v
                sleep(idle/1000)

        if stop:
            break

    # get all remaining processes
    while len(queue) > processes:
        for xi_,m_,v in queue.join():
            M[d[xi_],m_] = v
        sleep(idle/1000)


    # set all matrix 0 which are > 0
    mask = np.where(M > 0)
    M[mask] = 0

    # m=0 term factor 1/2
    M[:,0] /= 2

    terms_integral = []
    for i,xi in enumerate(xk):
        value = fsum(M[i,:])
        terms_integral.append( wk[i]*np.exp(xk[i])*value )
        print("# k=%d, x=%g, logdetD(xi=x/alpha)=%.15g" % (i, xi, value))

    integral = fsum(terms_integral)/alpha/np.pi
    print("#")
    print("# L/R, lmax, order, alpha, F(T=0)*(L+R)/(ħc)")
    print("%g, %d, %d, %.15g, %.15g" % (LbyR, lmax, deg, alpha, integral))
