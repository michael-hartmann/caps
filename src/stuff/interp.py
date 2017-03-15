import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt


def slurp(filename):
    data = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if "# k=" in line:
                line = line[2:]
                line = line.replace("k=", "")
                line = line.replace("logdetD=", "")
                line = line.replace("t=", "")
                line = line.replace("xi=", "")
                data.append(list(map(float, line.split(","))))
    return data


def interp(x,y, eps=1e-7):
    f_linear = interpolate.interp1d(x, y, kind="linear")
    f_cubic  = interpolate.interp1d(x, y, kind="cubic")

    def g(x):
        lin = f_linear(x)
        if abs(lin) < eps:
            return lin

        cub = f_cubic(x)
        if cub > 0:
            return lin
        else:
            return cub

    return np.vectorize(g)


# k, xi, logdetD, t
data = slurp("foo.out")

x = []
y = []
for k,xi,logdetD,t in data:
    x.append(xi)
    y.append(logdetD)

f = interp(x,y)

xnew = np.arange(0, 0.5, 0.00001)
ynew = f(xnew)
plt.plot(x, y, 'o', xnew, ynew, '-')
plt.xlim([0,0.5])
plt.show()
