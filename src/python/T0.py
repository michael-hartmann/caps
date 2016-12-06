import numpy as np
import libcasimir
from multiprocessing import Pool
from multiprocessing.context import TimeoutError
from random import random
from time import sleep

class Queue:
    def __init__(self,processes=4):
        self.queue = {}
        self.pool = Pool(processes=processes)

    def submit(self,xi,m,lmax):
        self.queue[(xi,m,lmax)] = self.pool.apply_async(logdetD, (xi,m,lmax))

    def join(self):
        d = []
        for key in list(self.queue.keys()):
            try:
                v = self.queue[key].get(False)
                del self.queue[key]
                xi,m,lmax = key
                d.append((xi,m,lmax,v))
            except TimeoutError:
                pass

        return d

    def __len__(self):
        return len(self.queue)


def logdetD(xi,m,lmax):
    casimir = libcasimir.Casimir(1,xi,lmax=lmax)
    return casimir.logdetD(1,m)


def accurate(l, eps):
    l0 = l[0]
    if l0 == 0:
        return False

    for j in range(len(l)):
        if l[j] == 0:
            if abs(l[j-1]/l0 ) < prec:
                return True
            else:
                return False


if __name__ == "__main__":
    prec = 1e-6
    LbyR = 0.1
    processes = 4
    deg = 80
    lmax = 40
    idle = 2 # in ms

    alpha = 2*LbyR/(1+LbyR);
    xk,wk = np.polynomial.laguerre.laggauss(deg)
    xk /= alpha

    d = {}
    for i,x in enumerate(xk):
        d[x] = i

    queue = Queue()

    # initialize results
    M = np.zeros((deg,lmax))

    for m in range(lmax):
        stop = True
        #print("m=%d" % m)
        for i,xi in enumerate(xk):
            if not accurate(M[i,:],prec):
                stop = False
                queue.submit(xi,m,lmax)

            while len(queue) >= processes:
                for xi_,m_,lmax_,v in queue.join():
                    print("# xi=%g, m=%d, v=%g" % (xi_,m_,v))
                    M[d[xi_],m_] = v

                sleep(idle/1000)

        if stop:
            break


    print("get all remaining processes")
    while len(queue) > processes:
        for xi_,m_,v in queue.join():
            print("# xi=%g, m=%d, v=%g" % (xi_,m_,v))
            M[d[xi_],m_] = v


    for i,xi in enumerate(xk):
        print(xi, np.sum(M[i,:]))
