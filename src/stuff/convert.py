from math import pi

hbar    = 1.0545718e-34   # m² kg / s
hbar_eV = 6.582119514e-16 # eV s
kB      = 1.38064852e-23  # m² kg / ( K s² )
c       = 299792458       # m/s

def scaled(L, R, T, ωp, γ):
    T_scaled  = 2*pi*kB*(L+R)*T/(hbar*c)
    ωp_scaled = ωp/(hbar_eV*c)*(L+R)
    γ_scaled  =  γ/(hbar_eV*c)*(L+R)

    return T_scaled, ωp_scaled, γ_scaled


if __name__ == "__main__":
    from sys import argv

    try:
        L  = float(eval(argv[1]))
        R  = float(eval(argv[2]))
        T  = float(eval(argv[3]))
        ωp = float(eval(argv[4]))
        γ  = float(eval(argv[5]))
    except:
        print("%s L R T ωp γ" % argv[0])
        exit(1)

    T_scaled, ωp_scaled, γ_scaled = scaled(L, R, T, ωp, γ)

    print("L  = %.15g" % L)
    print("R  = %.15g" % R)
    print("T  = %.15g" % T_scaled)
    print("ωp = %.15g" % ωp_scaled)
    print("γ  = %.15g" % γ_scaled)
