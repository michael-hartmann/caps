from glob import glob
from sys import path,stderr
path.append("/cfs/home/h/a/hartmmic/git/libcasimir-dev/src/python")
from PFA import PFA, hbarc, hbar_eV

def slurp(filename):
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if len(line) and line[0] != "#":
                LbyR, L, R, T, ldim, F_ = map(float, line.split(","))
                return L, R, T, F_

    raise LookupError("File does not contain results")

if __name__ == "__main__":
    #theta1 = -1.8325715+/-0.0000026
    theta1 = -1.8325715

    T = 0
    L = 150e-9
    omegap_eV = 9
    gamma_eV  = 0.035

    omegap = omegap_eV/hbar_eV # in rad/s
    gamma  = gamma_eV/hbar_eV  # in rad/s
    epsm1 = lambda xi: omegap**2/(xi*(xi+gamma))

    data = []
    for filename in glob("raw/slurm-*.out"):
        try:
            L, R, T, E = slurp(filename)
            data.append((L,R,T,E))
        except LookupError:
            print("%s did not contain result" % filename, file=stderr)

    print("# Drude model, T=%gK, L=%gnm, omegap=%geV, gamma=%gmeV, theta1=%.8g" % (T, L*1e9, omegap_eV, gamma_eV*1000, theta1))
    print("# L/R, E/E_PFA-1, E/E_PFA-1-theta1*(L/R)")
    for L,R,T,E in sorted(data):
        pfa = PFA(R, T, epsm1)
        E_PFA = (L+R)/hbarc*pfa.E(L)
        print(L/R, E/E_PFA-1, E/E_PFA-1-theta1*L/R)