from math import pi

hbar    = 1.0545718e-34   # m² kg / s
hbar_eV = 6.582119514e-16 # eV s
kB      = 1.38064852e-23  # m² kg / ( K s² )
c       = 299792458       # m/s

def scaled(R, LbyR, T, ωp, γ):
    T_scaled  = 2*pi*kB*R*(1+LbyR)*T/(hbar*c)
    ωp_scaled = ωp/(hbar_eV*c)*R*(1+LbyR)
    γ_scaled  =  γ/(hbar_eV*c)*R*(1+LbyR)

    return T_scaled, ωp_scaled, γ_scaled


if __name__ == "__main__":
    from sys import argv

    if len(argv) < 6:
        print("%s R L/R T ωp γ" % argv[0])
        exit(1)

    R    = eval(argv[1])
    LbyR = eval(argv[2])
    T    = eval(argv[3])
    ωp   = eval(argv[4])
    γ    = eval(argv[5])

    T_scaled, ωp_scaled, γ_scaled = scaled(R, LbyR, T, ωp, γ)

    print("R   = %.15g" % R)
    print("L/R = %.15g" % LbyR)
    print("T   = %.15g" % T_scaled)
    print("ωp  = %.15g" % ωp_scaled)
    print("γ   = %.15g" % γ_scaled)
