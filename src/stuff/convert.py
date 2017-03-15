def eV2scaled(R, LbyR, ωp, γ):
    ħ = 6.582119514e-16 # reduced Planck constant [eV/s]
    c = 299792458 # speed of light

    ωp_scaled = ωp/(ħ*c)*R*(1+LbyR)
    γ_scaled  =  γ/(ħ*c)*R*(1+LbyR)

    return ωp_scaled, γ_scaled


if __name__ == "__main__":
    from sys import argv

    if len(argv) < 5:
        print("%s R L/R ωp γ" % argv[0])
        exit(1)

    R    = float(argv[1])
    LbyR = float(argv[2])
    ωp   = float(argv[3])
    γ    = float(argv[4])

    print("%.15e %.15e" % eV2scaled(R, LbyR, ωp, γ))
