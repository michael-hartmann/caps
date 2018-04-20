from math import pi

# coefficients β_ij, Table I, [Teo, PRD 88, 045019 (2013)]
d_β = {
    (0,0): +1,
    (1,0): -4/3,
    (2,0): +9/5,
    (1,1): +18/5,
    (3,0): -16/7+32/735*pi**2,
    (2,1): -48/7,
    (4,0): +25/9-326/1323*pi**2,
    (3,1): +100/9-326/1323*pi**2,
    (2,2): +50/3,
    (5,0): -36/11+1220/1617*pi**2-379/32340*pi**4,
    (4,1): -180/11+2440/1617*pi**2,
    (3,2): -360/11+1220/1617*pi**2
}

# coefficients λ_ij, Table II, [Teo, PRD 88, 045019 (2013)]
d_λ = {
    (0,0): -20/pi**2+1/3,
    (1,0): +56/(3*pi**2)-32/45,
    (0,1): +56/(3*pi**2)-14/45,
    (2,0): -398/21/pi**2+401/315,
    (1,1): -796/21/pi**2+454/315,
    (0,2): -398/21/pi**2+113/315,
    (3,0): +410/21/pi**2-37/18+286/6615*pi**2,
    (2,1): +410/7/pi**2-26/7,
    (1,2): +410/7/pi**2-16/7,
    (0,3): +410/21/pi**2-79/126+pi**2/6615,
    (4,0): -69824/3465/pi**2+35141/10395-28022/99225*pi**2,
    (3,1): -279296/3465/pi**2+84176/10395-2774/14175*pi**2,
    (2,2): -139648/1155/pi**2+742/99+32/11025*pi**2,
    (1,3): -279296/3465/pi**2+43856/10395-46558/1091475*pi**2,
    (0,4): -69824/3465/pi**2+14981/10395-11962/1091475*pi**2,
    (5,0): +26732/1287/pi**2-150368/27027+4937399/5675670*pi**2-1142/63063*pi**4,
    (4,1): +133660/1287/pi**2-35026/2079+773884/567567*pi**2,
    (3,2): +267320/1287/pi**2-548024/27027+26212/51597*pi**2,
    (2,3): +267320/1287/pi**2-415724/27027+16826/81081*pi**2,
    (1,4): 133660/1287/pi**2-256888/27027+19984/81081*pi**2,
    (0,5): 26732/1287/pi**2-84218/27027+3329/62370*pi**2+8059/2522520*pi**4,
    (5,0): 26732/1287/pi**2-150368/27027+4937399/5675670*pi**2-1142/63063*pi**4
}

def β(i,j):
    """Return coefficient β or 0 if i+j>5"""
    if i >= j:
        key = (i,j)
    else:
        key = (j,i)
    
    return d_β.get(key,0)


def λ(i,j):
    """Return coefficient λ or 0 if i+j>5"""
    key = (i,j)
    return d_λ.get(key,0)


def E_plasma(a):
    """Compute the free energy in the plasma model.

    Parameters: a=c/(ω_P*L)

    E ≈ -π³ħcR/(720L²)*(s_1 + s_2 L/R)

    References: Teo, PRD 88, 045019 (2013)

    Returns: s_1, s_2
    print(E_plasma(float(argv[1])))
    """
    s1 = 0
    s2 = 0

    maximum = 5
    for i in range(maximum+1):
        for j in range(maximum+1):
            s1 += β(i,j)*a**(i+j)
            s2 += λ(i,j)*a**(i+j)

    return s1,s2


if __name__ == "__main__":
    from sys import argv, stderr
    from math import pi
    from PFA import PFA, hbar, hbar_eV, c

    try:
        a = float(argv[1])
    except:
        print("Usage: %s c/(ω_P*L)" % argv[0], file=stderr)
        exit(1)

    R = 0.11
    L = 0.01
    omegap = c/(a*L)
    E_PR = pi**3*hbar*c*R/(720*L**2)

    epsm1 = lambda xi: (omegap/xi)**2
    pfa = PFA(R, 0, epsm1)
    E = pfa.E(L)
    theta0_besser = E/E_PR

    theta0,theta1 = E_plasma(a)
    print("# theta0, theta1, theta1/theta0")
    print(theta0, theta1, theta1/theta0, theta1/theta0_besser)
