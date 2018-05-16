import numpy as np
from data import teo_theta1

# python fit.py */data.csv

def theta_fit(filename, fit_min=0, fit_max=0.002001):
    LbyR = []
    corr = []
    with open(filename) as f:
        inva = float("nan")
        for line in f:
            line = line.strip()
            if "# omega_p*L/c" in line:
                inva = float(line.split("=")[1])
            elif len(line) and line[0] != "#":
                d = list(map(float, line.split(",")))
                x,c = d[0], d[-1]

                if fit_min <= x <= fit_max:
                    LbyR.append(x)
                    corr.append(c)

    LbyR = np.array(LbyR)
    corr = np.array(corr)

    x = np.sqrt(LbyR)
    y = corr/LbyR**1.5
    theta3,theta2 = np.polyfit(x,y,1)
    theta1 = teo_theta1(1/inva)
    return inva, theta1, theta2, theta3


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fit fit parameters θ_1, θ_2, θ_3")
    parser.add_argument("filenames", nargs="+", help="data.csv")
    parser.add_argument("--fit_min", help="points used for fit: FIT_MIN ≤ L/R ≤ FIT_MAX", type=float, default=0)
    parser.add_argument("--fit_max", help="points used for fit: FIT_MIN ≤ L/R ≤ FIT_MAX", type=float, default=0.0011)
    args = parser.parse_args()

    l = []
    for filename in args.filenames:
        l.append(theta_fit(filename, fit_min=args.fit_min, fit_max=args.fit_max))

    print("# L*omegap/c, theta1, theta2, theta3")
    for inva, theta1, theta2, theta3 in sorted(l):
        print("%.10g, %.10g, %.10g, %.10g" % (inva, theta1, theta2, theta3))
