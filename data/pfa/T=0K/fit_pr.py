# get fitting parameters for PR

import numpy as np
from scipy.optimize import curve_fit
from sys import argv

theta1 = 1/3 - 20/np.pi**2

# PFA formula
def pfa(LbyR):
    return -np.pi**3/720.*(LbyR+1)/LbyR**2

left  = 0.002
right = 0.008

try:
    left  = float(argv[1])
    right = float(argv[2])
except:
    pass

data = np.loadtxt("pr/data_eta10.csv", delimiter=",", usecols=(0,3))

xdata = []
ydata = []

for i,(LbyR,corr) in enumerate(data):
    if left <= LbyR <= right:
        y = corr/pfa(LbyR)-1-theta1*LbyR

        xdata.append(LbyR)
        ydata.append(y)

xdata = np.array(xdata)
ydata = np.array(ydata)

y = ydata/xdata**1.5
x = np.sqrt(xdata)
theta3_polyfit, theta2_polyfit = np.polyfit(x,y,1)

func = lambda x,theta2, theta3: theta2*x**1.5+theta3*x**2

popt, pcov = curve_fit(func, xdata, ydata, p0=(3, 3))
theta2_curve_fit = popt[0]
theta3_curve_fit = popt[1]

print("E ≈ E_PFA ( 1 + θ_1·x + θ_2·x^1.5 + θ_3·x² )")
print("θ_1 = 1/3-20/π² ≈ %+.12g" % theta1)
print()
print("parameters obtained using linear regression (%g ≤ L/R ≤ %g):" % (left,right))
print("\tθ_2 = %+.8g" % theta2_polyfit)
print("\tθ_3 = %+.8g" % theta3_polyfit)
print()
print("parameters obtained using curve fit (%g ≤ L/R ≤ %g):" % (left,right))
print("\tθ_2 = %+.8g" % theta2_curve_fit)
print("\tθ_3 = %+.8g" % theta3_curve_fit)
